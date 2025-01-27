# 2024.07.25.
# Emese Bata & Anna Lovrics

library(readxl) # another package to read xls files
library(dplyr)
library(tidyr) # to transform data for ggplot2 plots
library(tibble)
library(ggplot2) # for visualization
library(ggpmisc) # for visualization
library(EnvStats) # for statistic tests
library(ggcorrplot) # for nice correlation plots
library(cowplot) # to plot q-q plots in a grid
library(caret)

phome <- "/media/anna/bigdisk/Anna/Enzimhome_local/Projects/StatisticsCore/Ras"
datafolder <- file.path(phome, "Data")
figfolder <- file.path(phome, "Figures")
if (!file.exists(figfolder)) dir.create(figfolder)
calcfolder <- file.path(phome, "Data_calculated")
if (!file.exists(calcfolder)) dir.create(calcfolder)
# source content of utility file
rscriptfolder <- file.path(phome, "Rscript")
util_file <- file.path(rscriptfolder, "anal_utils.R")
source(util_file)

# read in ras data
ras_file <- file.path(calcfolder, "ras_all.rds")
ras_data <- readRDS(ras_data, file=ras_file)

set.seed(74)
train_rows <- createDataPartition(
  as.factor(paste(ras_data$sample_type2, ras_data$cancer_type, sep="_")), 
  p = 0.8, list = FALSE)
indata_train <- ras_data[train_rows,]
indata_test <- ras_data[-train_rows,]

indep_column_indices <- 3:(ncol(ras_data)-2)
depvar="sample_type2"
input_x <- as.matrix(indata_train[,indep_column_indices])
input_y <- factor(indata_train[,depvar, drop=T])
holdout_x <- as.matrix(indata_test[,indep_column_indices])
holdout_y <- factor(indata_test[,depvar, drop=T])

# define train_control
train_control <- caret::trainControl(
  method = "cv", # cross-validation
  number = 3, # with n folds 
  verboseIter = FALSE, # no training log
  allowParallel = FALSE # FALSE for reproducible results 
)

# define default hyperparameters
grid_default <- expand.grid(
  nrounds = 1000,
  max_depth = 3,
  eta = 0.3,
  gamma = 0.9,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

xgb_base <- caret::train(
  x = input_x,
  y = input_y,
  trControl = train_control,
  tuneGrid = grid_default,
  method = "xgbTree",
  verbose = TRUE
)

# check model performance
base_train_fit <- predict(xgb_base, newdata=input_x)
conf_mx <- confusionMatrix(data=base_train_fit, # a factor of predicted classes
                           reference=input_y, # a factor of true classes
                           positive="Tumor")


base_test_fit <- predict(xgb_base, newdata=holdout_x)
conf_mx <- confusionMatrix(data=base_test_fit, 
                           reference=holdout_y,
                           positive="Tumor")

