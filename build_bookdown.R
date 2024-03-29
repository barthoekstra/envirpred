library(bookdown)
rm(list = ls(all = TRUE))
ggplot2::set_last_plot(NULL)
gc()
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
