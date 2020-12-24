library(mboost)
library(tidyverse)

model <- readRDS("data/processed/models/mod_resids.RDS")
outname <- "cvm_resids"

# Determine and save cross-validation folds
cv10f <- cv(model.weights(model), type = "kfold")
saveRDS(cv10f, file = "data/processed/models/cv10f_resids.RDS")

# Bootstrap optimal boosting iteration mstop
cvm <- cvrisk(model, folds = cv10f, grid = seq(from = 250, to = 30000, by = 250))
saveRDS(cvm, file = paste0("data/processed/models/", outname, ".RDS"))