library(bookdown)
rm(list = ls(all = TRUE))
# repro modes
full_repro <- FALSE
gee_repro <- FALSE
gdal_repro <- FALSE
dataprep_repro <- TRUE
model_repro <- TRUE

if (full_repro) {
  gee_repro <- TRUE
  gdal_repro <- TRUE
  dataprep_repro <- TRUE
  model_repro <- TRUE
}
