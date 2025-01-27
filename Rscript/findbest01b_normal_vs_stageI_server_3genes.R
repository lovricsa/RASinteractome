# 2024.10.01.
# Anna Lovrics

library(readxl) # another package to read xls fileslibrary(dplyr)
library(tidyr) # to transform data for ggplot2 plots
library(data.table) # while I use tibble mostly, fread is best input option
library(tibble)
library(ggplot2) # for visualization
library(caret) # for Classification And REgression Training
library(limma) # for the strsplit2 function
library(gridExtra) # to plot a list of figures in a grid
library(ggpubr) # for the same purpose to arrange figures in a grid
library(optparse)

option_list <- list( 
  make_option(c("--inID"), type="integer", 
              help="index of cancer type", default=2),
  make_option(c("--num_genes"), type="integer", default=3)
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

# now select cancer type
rtype_vector <- c("BRCA", "COAD", "LUAD", "LUSC")
rtype <- rtype_vector[opt$inID]
print(rtype)
figfolder <- file.path(figfolder, rtype)
calcfolder <- file.path(calcfolder, rtype)

# source content of utility file
rscriptfolder <- file.path(phome, "Rscript")
util_file <- file.path(rscriptfolder, "anal_utils.R")
source(util_file)

# use opt parse to parse input number
suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
option_list <- list( 
  make_option(c("--inID"), type="integer", 
              help="index of cancer type")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

# now select cancer type
rtype_vector <- c("BRCA", "COAD", "LUAD", "LUSC")
rtype <- rtype_vector[opt$inID]
print(rtype)
figfolder <- file.path(figfolder, rtype)
calcfolder <- file.path(calcfolder, rtype)

# read data
ras_list <- readRDS(file=file.path(calcfolder,
                                   paste("TCGA_", rtype, ".rds", sep="")))
ras_data <- ras_list[["data"]]
indep_vars <- ras_list[["indep_vars"]]

# select normal and stage samples
ras_data <- ras_data[ras_data$category2 %in% c("Normal", "Tumor_Stage_I"),]
# also delete unused factors
ras_data$category_factor <- droplevels(ras_data$category_factor)

# calculate aic of possible three gene models and 
# also the balanced accuracy of a confusion matrix obtained by cross-validation
set.seed(7004)

n <- length(indep_vars)
k <- opt$num_genes

num_models <- 100
mod_list_aic <- as.list(rep(NA, times=num_models))
aic_vals <- 1001:(1000+num_models)
mod_list_bic <- as.list(rep(NA, times=num_models))
bic_vals <- 1001:(1000+num_models)
mod_list_acc <- as.list(rep(NA, times=num_models))
acc_vals <- rep(0, times=num_models)
for (i in 1:choose(n, k)){
  #for (i in 1:1000){
  if (i == 1){
    cbn <- 1:k
  }else{
    cbn <- gen.next.cbn(cbn, n)
  }
  it_model <-  glm(as.formula(paste0("category_factor ~ ",
                          paste(indep_vars[cbn],collapse="+"))),
              data = ras_data, family = "binomial")
  # it_model <- glm(as.formula(sample_type_factor ~ .), 
  #                 data = ras_data[, c("sample_type_factor", indep_vars[cbn])], 
  #                 family = "binomial")
  it_aic <- it_model$aic
  #it_aic <- BIC(it_model)
  # keep best models based on aic
  if (it_aic < max(aic_vals)){
    change_ind <- which(aic_vals == max(aic_vals))[1]
    aic_vals[change_ind] <- it_aic
    mod_list_aic[[change_ind]] <- it_model
  }
  it_bic <- BIC(it_model)
  # keep best models based on bic
  if (it_bic < max(bic_vals)){
    change_ind <- which(bic_vals == max(bic_vals))[1]
    bic_vals[change_ind] <- it_bic
    mod_list_bic[[change_ind]] <- it_model
  }
  it_acc <- runCV_LOOCV(ras_data, it_model,
                        depvar="category_factor",
                        depvar_levels=c("Normal", "Tumor_Stage_I"),
                        posclass="Tumor_Stage_I",
                        #numCV_outer=5,
                        numCV_inner=3,
                        verbose=F,
                        return_conf_mx=F)
  
  # keep best models based on acc
  if (it_acc > min(acc_vals)){
    change_ind <- which(acc_vals == min(acc_vals))[1]
    acc_vals[change_ind] <- it_acc
    mod_list_acc[[change_ind]] <- it_model
  }
  if (i %% 100 == 0){
    print(paste(as.character(i), as.character(choose(n, k)), sep="/"))
  }
}

# save obtained results
resfile <- file.path(resfolder, 
            paste(rtype, "best_three_gene_models_normal_vs_stageI.rds", sep=""))
saveRDS(list("mod_list_aic"=mod_list_aic,
             "aic_vals"=aic_vals,
             "mod_list_bic"=mod_list_bic,
             "bic_vals"=bic_vals,
             "mod_list_acc"=mod_list_acc,
             "acc_vals"=acc_vals),
        file=resfile)

