library(mboost)
library(tidyverse)

model_name <- "trends"

# Load model
model <- readRDS(paste0("data/processed/models/mod_", model_name, ".RDS"))

# Determine and save cross-validation folds
cv10f <- cv(model.weights(model), type = "kfold")
saveRDS(cv10f, file = paste0("data/processed/models/cv10f_", model_name, ".RDS"))

# Bootstrap optimal boosting iteration mstop
cvm <- cvrisk(model, folds = cv10f, grid = seq(from = 250, to = 30000, by = 250))
saveRDS(cvm, file = paste0("data/processed/models/cvm_", model_name, ".RDS"))