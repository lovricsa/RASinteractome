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


#phome <- "/media/anna/bigdisk/Anna/Enzimhome_local/Projects/StatisticsCore/Ras"
phome <- "/bigdisk/users/alovrics/Projects/StatisticsCore/Ras"
datafolder <- file.path(phome, "Data")
figfolder <- file.path(phome, "Figures")
if (!file.exists(figfolder)) dir.create(figfolder)
calcfolder <- file.path(phome, "Data_calculated")
if (!file.exists(calcfolder)) dir.create(calcfolder)
resfolder <- file.path(phome, "Results")
if (!file.exists(resfolder)) dir.create(resfolder)

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
ras_data <- ras_data[ras_data$category2 %in% c("Tumor_Stage_I", "Tumor_Stage_IV"),]
# also delete unused factors
ras_data$category_factor <- droplevels(ras_data$category_factor)

# calculate aic of possible three gene models and 
# also the balanced accuracy of a confusion matrix obtained by cross-validation
set.seed(7004)

n   <- length(indep_vars)
k   <- 2

num_models <- choose(n,k)
mod_list <- list()
aic_vals <- 1001:(1000+num_models)
bic_vals <- 1001:(1000+num_models)
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
  mod_list[[i]] <- it_model
  aic_vals[i] <- it_model$aic
  bic_vals[i] <- BIC(it_model)
  acc_vals[i] <-  runCV_LOOCV(ras_data, it_model,
                        depvar="category_factor",
                        depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                        posclass="Tumor_Stage_IV")
  if (i %% 100 == 0){
    print(paste(as.character(i), as.character(choose(n, k)), sep="/"))
  }
}
# save obtained results
resfile <- file.path(resfolder, 
            paste(rtype, "best_two_gene_models_stageI_vs_stageIV.rds", sep=""))
saveRDS(list("mod_list"=mod_list,
             "aic_vals"=aic_vals,
             "bic_vals"=bic_vals,
             "acc_vals"=acc_vals),
        file=resfile)

