library(tidyverse)
library(mlr3)
library(mlr3learners)
library(mlr3extralearners)
library(mlr3spatiotempcv)
library(mlr3pipelines)
library(mlr3tuning)
library(mlr3viz)
library(paradox)

data <- readRDS("data/processed/data.RDS")
lut <- readRDS("data/processed/ecoregions_lut.RDS")

unused_predictors <- c("trend_long", "trend_short", "count_long", "count_short",
                       "def", "tmmx", "tmmn", "aet", "pet", "pr", "ndvi", "x", "y",
                       lut$ecoregions, lut$biomes, lut$realms)
predictors <- colnames(data)[!colnames(data) %in% unused_predictors]

task <- TaskRegrST$new(id = "data", backend = data, target = "trend_long", extra_args = list(coordinate_names = c("x", "y")))
task$select(all_of(predictors))

learners <- list(
  lrn("regr.gamboost", id = "gamboost"),
  lrn("regr.gbm", id = "gbm"),
  lrn("regr.xgboost", id = "xgboost")
)
learners_ids <- sapply(learners, function(x) x$id)

gr_subsample <- po("subsample") %>>%
  po("branch", options = learners_ids) %>>%
  gunion(lapply(learners, po)) %>>%
  po("unbranch")

grlrn_subsample <- GraphLearner$new(gr_subsample)

min_nrows <- 100
max_nrows <- 50000

min_nrows_s <- log10(min_nrows) / log10(2)
max_nrow_s <- log10(max_nrows) / log10(2)

ps_subsample <- ParamSet$new(list(
  ParamDbl$new("subsample.frac", lower = min_nrows_s, upper = max_nrow_s),
  ParamFct$new("branch.selection", levels = learners_ids)
))

# Distribution of subsample.frac in grid can be visualised as follows
# x <- seq(from = min_nrows_s, to = max_nrow_s, length.out = 20)
# y <- 2^x
# plot(x, y)
# Beware: division by nrow(data) has to be hard-coded in the function below, somehow...

ps_subsample$trafo <- function(x, param_set) {
  x$subsample.frac <- (2^x$subsample.frac) / 314868
  return(x)
}

resample_inner <- rsmp("cv", folds = 5)
resample_outer <- rsmp("cv", folds = 5)
measure <- msr("regr.rmse")
terminator <- trm("none")
tuner <- tnr("grid_search", resolution = 20)

at <- AutoTuner$new(
  learner = grlrn_subsample,
  resampling = resample_inner,
  measure = measure,
  search_space = ps_subsample,
  terminator = terminator,
  tuner = tuner,
  store_tuning_instance = TRUE,
  store_benchmark_result = TRUE
)

future::plan("multisession", workers = 3)
rr <- resample(task = task, learner = at, resampling = resample_outer, store_models = TRUE)
saveRDS(as.data.table(rr), file = "data/processed/rr_dataset_size.RDS")

rr_data <- bind_rows(lapply(rr$learners, function(x) x$model$tuning_instance$archive$data()))
saveRDS(rr_data, file = "data/processed/rr_dataset_size_df.RDS")
