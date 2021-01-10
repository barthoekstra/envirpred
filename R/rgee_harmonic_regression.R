# Most of the harmonic regression procedure is based on the excellent time
# series analysis lab in GEE by the Google team (see below). Modifications
# include dynamic naming of output bands, explicit detrending, and tweaks to
# suit this specific use-case.
# See the following url for the GEE time series analyis lab:
# https://docs.google.com/document/d/1mNIRB90jwLuASO1JYas1kuOXCLbOoy1Z4NlV1qIXM10/edit

rgee_harmonic_regression <- function(dataset, dependent, harmonics = 1L, aoi,
                                     resolution = 10000, dsn = NULL, 
                                     mask = NULL, refdataset = NULL) {
  aoi <- sf_as_ee(aoi)
  
  dataset <- ee$ImageCollection(dataset)$
    select(dependent)
  
  if (!is.null(refdataset)) {
    ref_daterange <- ee$List(ee$ImageCollection(refdataset)$get('date_range'))
    dataset <- dataset$
      filterDate(ref_daterange$get(0L), ref_daterange$get(1L))
  }
  
  # The dependent variable we're modelling
  dependent <- dependent
  
  # The number of cycles per year to model
  harmonics <- harmonics
  
  # Make a list of harmonic frequencies to model.
  # These also serve as band name suffixes.
  harmonicFrequencies <- ee$List$sequence(1L, harmonics)
  
  # Function to get a sequence of band names for harmonic terms.
  getNames <- function(base, list) {
    return(ee$List(list)$map(
      ee_utils_pyfunc(
        function(i) {
          return(ee$String(base)$cat(ee$Number(i)$format('%d')))
        }
      )
    ))
  }
  
  # Construct lists of names for the harmonic terms
  cosNames <- getNames('cos_', harmonicFrequencies)
  sinNames <- getNames('sin_', harmonicFrequencies)
  
  # Independent variables.
  independents = ee$List(c('constant', 't'))$
    cat(cosNames)$
    cat(sinNames)
  
  # Function to add a constant band.
  addConstant <- function(image) {
    return(image$addBands(ee$Image$constant(1)))
  }
  
  # Function to add a time band.
  addTime <- function(image) {
    # Compute time in fractional years since the epoch.
    date <- ee$Date(image$get('system:time_start'))
    years <- date$difference(ee$Date('1970-01-01'), 'year')
    timeRadians <- ee$Image(years$multiply(2L * pi))
    return(image$addBands(timeRadians$rename('t')$float()))
  }
  
  # Function to compute the specified number of harmonics and add them as bands.
  # Assumes the time band is present.
  addHarmonics <- function(freqs) {
    return(function(image) {
      # Make an image of frequencies.
      frequencies <- ee$Image$constant(freqs)
      # This band should represent time in radians.
      time <- ee$Image(image)$select('t')
      # Get the cosine terms.
      cosines <- time$multiply(frequencies)$cos()$
        rename(cosNames)
      # Get the sine terms.
      sines <- time$multiply(frequencies)$sin()$
        rename(sinNames)
      return(image$addBands(cosines)$addBands(sines))
    })
  }
  
  # Datasets are global, so we need to clip.
  clipToAoi <- function(image) {
    return(image$clip(aoi))
  }
  
  # Mask out certain values (e.g. water bodies where NDVI is -1).
  maskValues <- function(image) {
    mask <- image$select(dependent)$gt(mask)
    return(image$updateMask(mask))
  }
  
  # Filter to the area of interest, add variables.
  harmonicDep <- dataset$
    map(addConstant)$
    map(addTime)$
    map(addHarmonics(harmonicFrequencies))$
    map(clipToAoi)
  
  if (!is.null(mask)) {
    harmonicDep <- harmonicDep$map(maskValues)
  }
  
  # Detrend data using linear regression.
  detrendIndependents <- ee$List(c('constant', 't'))
  
  linearTrend <- harmonicDep$select(detrendIndependents$add(dependent))$
    reduce(ee$Reducer$linearRegression(detrendIndependents$length(), 1L))
  
  linearTrendCoefficients <- linearTrend$select('coefficients')$
    arrayProject(list(0))$
    arrayFlatten(c(detrendIndependents))$
    rename(ee$List(c('intercept', 'slope')))$
    float()
  
  detrendedHarmonicDep <- harmonicDep$map(
    function(image) {
      return(image$select(dependent)$subtract(
        image$select(detrendIndependents)$multiply(linearTrendCoefficients)$reduce('sum')
      ))$
        rename(dependent)$
        copyProperties(image, ee$List(c('system:time_start')))
    }
  )
  
  # Combine detrended and harmonic image collections
  harmonicDep <- harmonicDep$combine(detrendedHarmonicDep)
  bandNames <- harmonicDep$first()$bandNames()
  lastBandID <- ee$Number(bandNames$length()$subtract(1L))
  origDependentName <- ee$String(dependent)$cat(ee$String('_original'))
  newBandNames <- bandNames$set(lastBandID, dependent)$set(0L, origDependentName)
  
  # Rename the bands to keep original along with detrended dependent
  harmonicDep <- harmonicDep$map(function(image) {
    return(image$rename(newBandNames))
  })
  
  # The output of the regression reduction is a 4x1 array image.
  harmonicTrend <- harmonicDep$
    select(independents$add(dependent))$
    reduce(ee$Reducer$linearRegression(independents$length(), 1L))
  
  harmonicTrendResiduals <- harmonicTrend$select('residuals')$
    arrayFlatten(list(list('residuals')))$
    float()
  
  # Turn the array into a multi-band image of coefficients.
  harmonicTrendCoefficients <- harmonicTrend$select('coefficients')$
    arrayProject(list(0))$
    arrayFlatten(c(independents))
  
  # Compute fitted values.
  fittedHarmonic <- harmonicDep$map(function(image) {
    return(image$addBands(
      image$select(independents)$
        multiply(harmonicTrendCoefficients)$
        reduce('sum')$
        rename('fitted')
    ))
  })
  
  # Calculate phase
  calculatePhase <- function(harmonic) {
    phase <- harmonicTrendCoefficients$
      select(ee$String('sin_')$cat(ee$Number(harmonic)$format('%d')))$
      atan2(harmonicTrendCoefficients$select(ee$String('cos_')$cat(ee$Number(harmonic)$format('%d'))))$
      rename(ee$String('phase')$cat(ee$Number(harmonic)$format('%d')))$
      float()
    
    phaseCorrector <- phase$add(2L * pi)
    phaseMask <- phase$lt(0L)
    phaseCorrected <- phase$
      where(phaseMask, phaseCorrector)$
      rename(ee$String('phase')$cat(ee$Number(harmonic)$format('%d'))$cat(ee$String('_cor')))$
      float()
    
    return(ee$Image$cat(c(phase, phaseCorrected)))
  }
  
  phase <- harmonicFrequencies$map(
    ee_utils_pyfunc(
      function(x) calculatePhase(x)
    )
  )
  
  phaseBandNames <- harmonicFrequencies$map(
    ee_utils_pyfunc(
      function(freq) {
        base <- ee$String('phase')$cat(ee$Number(freq)$format('%d'))
        base_corr <- base$cat(ee$String('_cor'))
        return(ee$List(c(base, base_corr)))
        # return(ee$List(c(base, base_corr)))
      }
    )
  )
  
  phases <- ee$ImageCollection$fromImages(phase)$
    toBands()$
    rename(phaseBandNames$flatten())
  
  # Calculate amplitude
  calculateAmplitude <- function(harmonic) {
    return(harmonicTrendCoefficients$
             select(ee$String('sin_')$cat(ee$Number(harmonic)$format('%d')))$
             hypot(harmonicTrendCoefficients$select(ee$String('cos_')$cat(ee$Number(harmonic)$format('%d'))))$
             rename(ee$String('amplitude')$cat(ee$Number(harmonic)$format('%d')))$
             float()
    )
  }
  
  amplitude <- harmonicFrequencies$map(
    ee_utils_pyfunc(
      function(harmonic) calculateAmplitude(harmonic)
    )
  )
  amplitudeBandNames <- harmonicFrequencies$map(
    ee_utils_pyfunc(
      function(freq) ee$String('amplitude')$cat(ee$Number(freq)$format('%d'))
    )
  )
  
  amplitudes <- ee$ImageCollection$fromImages(amplitude)$
    toBands()$
    rename(amplitudeBandNames)
  
  # Combine bands in a single image
  origMeanDep <- harmonicDep$select(origDependentName)$mean()$float()
  
  envirpred <- origMeanDep$
    addBands(harmonicTrendResiduals)$
    addBands(phases)$
    addBands(amplitudes)$
    addBands(linearTrendCoefficients)
  
  out <- ee_as_stars(envirpred, region = aoi$geometry(), via = "drive", scale = resolution,
                     dsn = dsn)
  
  return(out)
}