# 2024.10.01.
# Anna Lovrics

library(readxl) # another package to read xls files
library(dplyr)
library(tidyr) # to transform data for ggplot2 plots
library(data.table) # while I use tibble mostly, fread is best input option
library(tibble)
library(ggplot2) # for visualization
library(caret) # for Classification And REgression Training
library(limma) # for the strsplit2 function
library(gridExtra) # to plot a list of figures in a grid
library(ggpubr) # for the same purpose to arrange figures in a grid
library(dplyr)

library(genalg) # for the genetic algorithm


phome <- "/media/anna/bigdisk/Anna/Enzimhome_local/Projects/StatisticsCore/Ras"
#phome <- "/bigdisk/users/alovrics/Projects/StatisticsCore/Ras"
datafolder <- file.path(phome, "Data")
figfolder <- file.path(phome, "Figures")
calcfolder <- file.path(phome, "Data_calculated")
resfolder <- file.path(phome, "Results")

# source content of utility file
rscriptfolder <- file.path(phome, "Rscript")
util_file <- file.path(rscriptfolder, "anal_utils.R")
source(util_file)

# test some models
compname <- "stageI_vs_stageIV"
rtype <- "COAD"
ras_list <- readRDS(file=file.path(calcfolder, rtype,
                                   paste("TCGA_", rtype, ".rds", sep="")))
ras_data <- ras_list[["data"]]
indep_vars <- ras_list[["indep_vars"]]
ras_data <- ras_data[ras_data$category2 %in% c("Tumor_Stage_I", "Tumor_Stage_IV"),]
# also delete unused factors
ras_data$category_factor <- droplevels(ras_data$category_factor)
indep_vars <- colnames(ras_data)[4:49]
suggest_genes <- c("RGL2", "RASSF1", "RADIL", "RASSF7", "RAPGEF4", "KRAS")
suggest_genes <- c("RGL2", "RASSF1", "RADIL", "RASSF7", "RAPGEF4", "KRAS", "RIN3")
suggest_binary <- rep(0, times=length(indep_vars))
suggest_binary[match(suggest_genes, indep_vars)] <- 1
suggest_binary <- matrix(suggest_binary, nrow=1, ncol=length(indep_vars))

calc_accurary <- function(indices){
  test_genes <- indep_vars[which(indices==1)]
  it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                      paste(test_genes,collapse="+"))),
                    data = ras_data, family = "binomial")
  balanced_accuracy <- runCV_LOOCV(ras_data, it_model,
                         depvar="category_factor",
                         depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                         posclass="Tumor_Stage_IV",
                         return_conf_mx=F)
  return(1-balanced_accuracy) # objective function is minimized
}
# 
# # also write a monitoring function
# monitor <- function(obj){
#   minEval = min(obj$evaluations);
#   # plot(obj, type="hist");
# }
# 
# rgba_res <- rbga.bin(size=length(indep_vars), # number of genes in 'chromosome'
#                   suggestions = suggest_binary,
#                   popSize=50, # default is 200 
#                   iters=50, # default is 100
#                   mutationChance=0.05, 
#                   elitism = 1,
#                   zeroToOneRatio=10,
#                   evalFunc=calc_accurary, 
#                   showSettings = TRUE,
#                   verbose=TRUE, 
#                   monitorFunc=monitor)
# 
# # save result
# resfile <- file.path(resfolder, 
#                 paste("genetic_best_", rtype, "_", compname, ".rds", sep=""))
# saveRDS(rgba_res, file=resfile)
# print(rgba_res)

# read result

resfile <- file.path(resfolder,
                paste("genetic_best_nosuggestion_", rtype, "_", compname, ".rds", sep=""))
rgba_res <- readRDS(file=resfile)
print(rgba_res)
best_res <- suggest_binary[1,]
calc_accurary(best_res)

bi <- which(rgba_res$evaluations==min(rgba_res$evaluations))
#for (ii in bi){
  ii <- bi[2]
  best_res <- rgba_res$population[ii,]
  test_genes <- indep_vars[which(best_res==1)]
  print(test_genes)
  it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                      paste(test_genes,collapse="+"))),
                    data = ras_data, family = "binomial")
  balanced_accuracy <- runCV_LOOCV(ras_data, it_model,
                                   depvar="category_factor",
                                   depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                                   posclass="Tumor_Stage_IV",
                                   return_conf_mx=F)
  print(balanced_accuracy)
  print("########")
#}

suggest_genes <- indep_vars[which(suggest_binary[1,]==1)]
print(suggest_genes)

# ####################################################
# compname <- "normal_vs_stageI"
# rtype <- "COAD"
# ras_list <- readRDS(file=file.path(calcfolder, rtype,
#                                    paste("TCGA_", rtype, ".rds", sep="")))
# ras_data <- ras_list[["data"]]
# ras_data <- ras_data[ras_data$category2 %in% c("Normal", "Tumor_Stage_I"),]
# # also delete unused factors
# ras_data$category_factor <- droplevels(ras_data$category_factor)
# 
# conf_mx <- runCV(ras_data, it_model,
#                  depvar="category_factor",
#                  depvar_levels=c("Normal", "Tumor_Stage_I"),
#                  posclass="Tumor_Stage_I",
#                  return_conf_mx=T)
# print(conf_mx)
# 
# 
# test_genes <- c("RIN1", "RAF1")
# it_model  <-  glm(as.formula(paste0("category_factor ~ ",
#                                     paste(test_genes,collapse="+"))),
#                   data = ras_data, family = "binomial")
# conf_mx <- runCV(ras_data, it_model,
#                  depvar="category_factor",
#                  depvar_levels=c("Normal", "Tumor_Stage_I"),
#                  posclass="Tumor_Stage_I",
#                  return_conf_mx=T)
# print(conf_mx)

