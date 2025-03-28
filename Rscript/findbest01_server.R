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

#phome <- "/media/anna/bigdisk/Anna/Enzimhome_local/Projects/StatisticsCore/Ras"
phome <- "/bigdisk/users/alovrics/Projects/StatisticsCore/Ras"
datafolder <- file.path(phome, "Data")
figfolder <- file.path(phome, "Figures")
calcfolder <- file.path(phome, "Data_calculated")
resfolder <- file.path(phome, "Results")

# source content of utility file
rscriptfolder <- file.path(phome, "Rscript")
util_file <- file.path(rscriptfolder, "anal_utils.R")
source(util_file)


option_list <- list( 
  make_option(c("--inID"), type="integer", 
              help="index of cancer type", default=2),
  make_option(c("--comparison_index"), type="integer", 
              help="index of comparison type", default=1),
  make_option(c("--num_genes"), type="integer", default=3)
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# # temporary to test
# opt <- list("inID"=2,
#             "comparison_index"=1,
#             "num_genes"=3)

# now select cancer type
rtype_vector <- c("BRCA", "COAD", "LUAD", "LUSC")
rtype <- rtype_vector[opt$inID]
print(rtype)
figfolder <- file.path(figfolder, rtype)
calcfolder <- file.path(calcfolder, rtype)

# and also select dependent variable plus levels
comptype_vector <- c("normal_vs_tumor", "normal_vs_stageI", 
                     "stageI_vs_stageIV")
comptype <- comptype_vector[opt$comparison_index]
print(comptype)
depvar_vector <- c("sample_type_factor", "category_factor", "category_factor")
depvar_levels_list <- list(c("Normal", "Tumor"),
                           c("Normal", "Tumor_Stage_I"),
                           c("Tumor_Stage_I", "Tumor_Stage_IV"))

# read data
ras_list <- readRDS(file=file.path(calcfolder,
                                   paste("TCGA_", rtype, ".rds", sep="")))
ras_data <- ras_list[["data"]]
indep_vars <- ras_list[["indep_vars"]]


# if needed select samples
if (!comptype == "normal_vs_tumor"){
  ras_data <- ras_data[ras_data$category2 %in% 
                        depvar_levels_list[[opt$comparison_index]],]
# also delete unused factors
ras_data$category_factor <- droplevels(ras_data$category_factor)
}

# calculate aic of possible three gene models and 
# also the balanced accuracy of a confusion matrix obtained by cross-validation
set.seed(7004)

n <- length(indep_vars)
k  <- opt$num_genes

num_models <- 100
mod_list <- as.list(rep(NA, times=num_models))
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
  it_model <-  glm(as.formula(paste0(depvar_vector[opt$comparison_index], " ~ ",
                          paste(indep_vars[cbn],collapse="+"))),
              data = ras_data, family = "binomial")
  # it_model <- glm(as.formula(sample_type_factor ~ .), 
  #                 data = ras_data[, c("sample_type_factor", indep_vars[cbn])], 
  #                 family = "binomial")
  it_aic <- it_model$aic
  it_bic <- BIC(it_model)
  it_acc <- runCV_LOOCV(ras_data, it_model,
                         depvar=depvar_vector[opt$comparison_index],
                         depvar_levels=depvar_levels_list[[opt$comparison_index]],
                         posclass=depvar_levels_list[[opt$comparison_index]][2],
                         numCV_inner=3,
                         verbose=F,
                         return_conf_mx=F)
  # keep best models based on acc
  if (it_acc > min(acc_vals)){
    change_ind <- which(acc_vals == min(acc_vals))[1]
    acc_vals[change_ind] <- it_acc
    aic_vals[[change_ind]] <- it_aic
    bic_vals[[change_ind]] <- it_bic
    mod_list[[change_ind]] <- it_model
  }
  if (i %% 100 == 0){
    print(paste(as.character(i), as.character(choose(n, k)), sep="/"))
  }
}

# save obtained results
resfile <- file.path(resfolder, 
                     paste(rtype, "best_", as.character(opt$num_genes), 
            "gene_models_", comptype,".rds", sep=""))
saveRDS(list("mod_list"=mod_list,
             "aic_vals"=aic_vals,
             "bic_vals"=bic_vals,
             "acc_vals"=acc_vals),
        file=resfile)

