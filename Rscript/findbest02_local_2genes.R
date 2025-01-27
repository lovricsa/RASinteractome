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
library(dplyr)


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

######################
thr_list <- list("normal_vs_stageI"=0.95,
                 "stageI_vs_stageIV"=0.7)
best_genes_list <- list()
for (compname in c("normal_vs_stageI", "stageI_vs_stageIV")){
  summ_df_list <- list()
  best_genes_list[[compname]] <- list()
  for (rtype in c("BRCA", "COAD", "LUAD", "LUSC")){
    # original data
    # read data
    ras_list <- readRDS(file=file.path(calcfolder, rtype,
                                       paste("TCGA_", rtype, ".rds", sep="")))
    ras_data <- ras_list[["data"]]
    indep_vars <- ras_list[["indep_vars"]]
    
    # select normal and stage samples
    if (compname=="normal_vs_stageI"){
      ras_data <- ras_data[ras_data$category2 %in% c("Normal", "Tumor_Stage_I"),]
    }
    if (compname=="stageI_vs_stageIV"){
      ras_data <- ras_data[ras_data$category2 %in% c("Tumor_Stage_I", "Tumor_Stage_IV"),]
    }
    # also delete unused factors
    ras_data$category_factor <- droplevels(ras_data$category_factor)
    
    ######################x
    # read obtained results
    resfile <- file.path(resfolder, 
                         paste(rtype,"best_two_gene_models_",
                         compname, ".rds", sep=""))
    reslist <- readRDS(file=resfile)
    
    # plot histogram of results
    aic_vals_df <- as.data.frame(reslist[["aic_vals"]])
    colnames(aic_vals_df) <- "aic"
    p1 <- ggplot(data=aic_vals_df, mapping=aes(x=aic))+
      geom_histogram(bins=50)
    #plot(p1)
    
    bic_vals_df <- as.data.frame(reslist[["bic_vals"]])
    colnames(bic_vals_df) <- "bic"
    p2 <- ggplot(data=bic_vals_df, mapping=aes(x=bic))+
      geom_histogram(bins=50)
    #plot(p2)
    
    # plot histogram of results
    acc_vals_df <- as.data.frame(reslist[["acc_vals"]])
    colnames(acc_vals_df) <- "acc"
    p3 <- ggplot(data=acc_vals_df, mapping=aes(x=acc))+
      geom_histogram(bins=50)+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
    #plot(p3)
    
    plot_all <- cowplot::plot_grid(p1, p2, p3)
    figfile <- file.path(figfolder, rtype, paste("best_model_2genes_",
                         compname, "_aic_bic_acc_counts.pdf", sep=""))
    ggsave(figfile, plot=plot_all, width=8, height=8)
    
    # create a summary dataframe for the best models
    summ_df <- as_tibble(cbind(aic_vals_df, bic_vals_df, acc_vals_df))
    summ_df$gene1 <- NA
    summ_df$gene2 <- NA
    # add corresponding genes
    for (ii in 1:length(reslist[["mod_list"]])){
      coefs <- reslist[["mod_list"]][[ii]]$coefficients
      genes <- sort(names(coefs)[2:length(coefs)])
      summ_df$gene1[ii] <- genes[1]
      summ_df$gene2[ii] <- genes[2]
    }
    # sort by accuracy values
    summ_df <- summ_df[order(summ_df$acc, decreasing=T),]
    best_genes <- unique(c(summ_df$gene1[which(summ_df$acc>thr_list[[compname]])],
                     summ_df$gene2[which(summ_df$acc>thr_list[[compname]])]))
    # print(rtype)
    # print(length(best_genes))
    # print(setdiff(indep_vars, best_genes))
    best_genes_list[[compname]][[rtype]] <- best_genes
    summ_df_list[[rtype]] <- summ_df
  }
  # write best models two file
  outfile <- file.path(resfolder, 
      paste("best_model_2genes_", compname, ".xlsx", sep=""))
  openxlsx2::write_xlsx(summ_df_list, file=outfile)
}

# test some models
# compname <- "stageI_vs_stageIV"
compname <- "stageI_vs_stageIV"
rtype <- "COAD"
ras_list <- readRDS(file=file.path(calcfolder, rtype,
                                   paste("TCGA_", rtype, ".rds", sep="")))
ras_data <- ras_list[["data"]]
indep_vars <- ras_list[["indep_vars"]]
ras_data <- ras_data[ras_data$category2 %in% c("Tumor_Stage_I", "Tumor_Stage_IV"),]
# also delete unused factors
ras_data$category_factor <- droplevels(ras_data$category_factor)
# resfile <- file.path(resfolder, 
#                      paste(rtype,"best_two_gene_models_",
#                            compname, ".rds", sep=""))
# reslist <- readRDS(file=resfile)
# it_model <- reslist[["mod_list"]][[which(
#            reslist[["acc_vals"]]==max(reslist[["acc_vals"]]))[1]]]
# it_aic <- it_model$aic
# it_acc_vec <- c(1:10)
# for (jj in 1:10){
#   it_acc <- runCV(ras_data, it_model,
#                   depvar="category_factor",
#                   depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
#                   posclass="Tumor_Stage_IV",
#                   return_conf_mx=F)
#   it_acc_vec[jj] <- it_acc
# }
# print(it_acc_vec)
# 
# print(c(it_aic, it_acc, BIC(it_model)))
# conf_mx <- runCV(ras_data, it_model,
#                  depvar="category_factor",
#                  depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
#                  posclass="Tumor_Stage_IV",
#                  return_conf_mx=T)
# print(conf_mx)

test_genes <- c("RASIP1", "RASSF1", "ARAP2", "GRB7", "RGL2")
test_genes <- c("RGL2", "ARAP2", "RASSF1", "RAPGEF4", "GRB7")
test_genes <- c("KRAS", "MLLT4", "SNX27", "NRAS", "HRAS")
test_genes <- c("RASSF7", "RAPGEF6", "KRAS", "RASSF1", "NRAS", "RAPGEF2")
test_genes <- c("RGL2", "RASSF1", "RADIL", "RASSF7", "RAPGEF4", "KRAS")
test_genes <- c("ARAP2", "RADIL")
it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                    paste(test_genes,collapse="+"))),
                  data = ras_tmp2, family = "binomial")
conf_mx <- runCV_LOOCV(ras_data, it_model,
                 depvar="category_factor",
                 depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                 posclass="Tumor_Stage_IV",
                 return_conf_mx=T,
                 verbose=T)
print(conf_mx)
############xxxx


compname <- "normal_vs_stageI"
rtype <- "COAD"
ras_list <- readRDS(file=file.path(calcfolder, rtype,
                                   paste("TCGA_", rtype, ".rds", sep="")))
ras_data <- ras_list[["data"]]
indep_vars <- ras_list[["indep_vars"]]
ras_data <- ras_data[ras_data$category2 %in% c("Normal", "Tumor_Stage_I"),]
# also delete unused factors
ras_data$category_factor <- droplevels(ras_data$category_factor)
resfile <- file.path(resfolder, 
                     paste(rtype,"best_two_gene_models_",
                           compname, ".rds", sep=""))
reslist <- readRDS(file=resfile)
it_model <- reslist[["mod_list"]][[which(
  reslist[["acc_vals"]]==max(reslist[["acc_vals"]]))[1]]]
# it_model_bests <- reslist[["mod_list_acc"]][which(
#   reslist[["acc_vals"]]==min(reslist[["acc_vals"]]))]
it_aic <- it_model$aic
it_acc_vec <- c(1:10)
for (jj in 1:10){
  it_acc <- runCV(ras_data, it_model,
                  depvar="category_factor",
                  depvar_levels=c("Normal", "Tumor_Stage_I"),
                  posclass="Tumor_Stage_I",
                  return_conf_mx=F)
  it_acc_vec[jj] <- it_acc
}
print(it_acc_vec)
print(c(it_aic, it_acc, BIC(it_model)))
conf_mx <- runCV(ras_data, it_model,
                 depvar="category_factor",
                 depvar_levels=c("Normal", "Tumor_Stage_I"),
                 posclass="Tumor_Stage_I",
                 return_conf_mx=T)
print(conf_mx)


test_genes <- c("RIN1", "RAF1")
it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                    paste(test_genes,collapse="+"))),
                  data = ras_data, family = "binomial")
conf_mx <- runCV(ras_data, it_model,
                 depvar="category_factor",
                 depvar_levels=c("Normal", "Tumor_Stage_I"),
                 posclass="Tumor_Stage_I",
                 return_conf_mx=T)
print(conf_mx)

