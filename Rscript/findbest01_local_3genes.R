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
for (compname in c("normal_vs_stageI", "stageI_vs_stageIV")){
  summ_df_list <- list()
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
                         paste(rtype,"best_three_gene_models_",
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
    figfile <- file.path(figfolder, rtype, paste("best_model_",
                         compname, "_aic_bic_acc_counts.pdf", sep=""))
    ggsave(figfile, plot=plot_all, width=8, height=8)
    
    # # now find the corresponding models
    # aic_best <- which(reslist[["aic_vals"]] < 80)
    # aic_best_models <- reslist[["mod_list_aic"]][aic_best]
    
    # create a summary dataframe for the best models
    summ_df <- as.data.frame(rep(0, length(indep_vars)))
    rownames(summ_df) <- indep_vars
    colnames(summ_df) <- "AIC"
    summ_df <- as_tibble(summ_df, rownames="gene")
    summ_df$BIC <- summ_df$AIC
    summ_df$ACC <- summ_df$AIC
    
    for (ii in 1:length(reslist[["mod_list_aic"]])){
      curr_mod <- reslist[["mod_list_aic"]][[ii]]
      sel_genes <- names(curr_mod$coefficients)[2:length(curr_mod$coefficients)]
      for (gene in sel_genes){
        gi <- which(summ_df$gene==gene)
        summ_df$AIC[gi] <- summ_df$AIC[gi]+1
      }
    }
    
    for (ii in 1:length(reslist[["mod_list_bic"]])){
      curr_mod <- reslist[["mod_list_bic"]][[ii]]
      sel_genes <- names(curr_mod$coefficients)[2:length(curr_mod$coefficients)]
      for (gene in sel_genes){
        gi <- which(summ_df$gene==gene)
        summ_df$BIC[gi] <- summ_df$BIC[gi]+1
      }
    }
    
    for (ii in 1:length(reslist[["mod_list_acc"]])){
      curr_mod <- reslist[["mod_list_acc"]][[ii]]
      sel_genes <- names(curr_mod$coefficients)[2:length(curr_mod$coefficients)]
      for (gene in sel_genes){
        gi <- which(summ_df$gene==gene)
        summ_df$ACC[gi] <- summ_df$ACC[gi]+1
      }
    }
    # sort by best models
    sortval <- summ_df$ACC + summ_df$AIC + summ_df$BIC
    summ_df <- summ_df[order(sortval, decreasing = T), ]
    
    # print(rtype)
    # print(summ_df)
    colnames(summ_df)[2:ncol(summ_df)] <- paste(rtype, 
                                                colnames(summ_df)[2:ncol(summ_df)], sep="_")
    summ_df_list[[rtype]] <- summ_df
  }
  # write results to file:
  # 1) try merging the results
  summ_df_all <- purrr::reduce(summ_df_list, full_join, by = "gene")
  sortval <- rowSums(select(summ_df_all, !gene))
  summ_df_all <- summ_df_all[order(sortval, decreasing = T),]
  # write to file
  outfile <- file.path(resfolder, paste("best_model_",
                       compname, "_aic_bic_acc_counts.csv", sep=""))
  write.csv(summ_df_all, file=outfile, row.names = F)
}

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
resfile <- file.path(resfolder, 
                     paste(rtype,"best_two_gene_models_",
                           compname, ".rds", sep=""))
reslist <- readRDS(file=resfile)
ri <- which( reslist[["acc_vals"]]==max(reslist[["acc_vals"]]))
it_model <- reslist[["mod_list"]][[ri[1]]]
it_aic <- it_model$aic
it_acc <- runCV(ras_data, it_model,
                depvar="category_factor",
                depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                posclass="Tumor_Stage_IV",
                return_conf_mx=F)
it_acc2 <- runCV_LOOCV(ras_data, it_model,
                       depvar="category_factor",
                       depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                       posclass="Tumor_Stage_IV",
                       return_conf_mx=F)
print(c(it_aic, it_acc, BIC(it_model)))
conf_mx <- runCV(ras_data, it_model,
                 depvar="category_factor",
                 depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                 posclass="Tumor_Stage_IV",
                 return_conf_mx=T)
print(conf_mx)

test_genes <- c("RASSF1", "RASSF7", "RASSF8")
test_genes <- c("ARAP2", "NRAS")
test_genes <- c("RASSF1", "RASSF7")
it_model  <-  glm(as.formula(paste0("sample_type_factor ~ ",
                                    paste(test_genes,collapse="+"))),
                  data = ras_data, family = "binomial")
conf_mx <- runCV(ras_data, it_model,
                 depvar="category_factor",
                 depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                 posclass="Tumor_Stage_IV",
                 return_conf_mx=T)
print(conf_mx)


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
                     paste(rtype,"best_three_gene_models_",
                           compname, ".rds", sep=""))
reslist <- readRDS(file=resfile)
it_model <- reslist[["mod_list_acc"]][[which(
  reslist[["acc_vals"]]==min(reslist[["acc_vals"]]))[1]]]
it_model_bests <- reslist[["mod_list_acc"]][which(
  reslist[["acc_vals"]]==min(reslist[["acc_vals"]]))]
it_aic <- it_model$aic
it_acc <- runCV(ras_data, it_model,
                depvar="category_factor",
                depvar_levels=c("Normal", "Tumor_Stage_I"),
                posclass="Tumor_Stage_I",
                return_conf_mx=F)
print(c(it_aic, it_acc, BIC(it_model)))
conf_mx <- runCV(ras_data, it_model,
                 depvar="category_factor",
                 depvar_levels=c("Normal", "Tumor_Stage_I"),
                 posclass="Tumor_Stage_I",
                 return_conf_mx=T)
print(conf_mx)


test_genes <- c("RIN1", "RAF1")
it_model  <-  glm(as.formula(paste0("sample_type_factor ~ ",
                                    paste(test_genes,collapse="+"))),
                  data = ras_data, family = "binomial")
conf_mx <- runCV(ras_data, it_model,
                 depvar="category_factor",
                 depvar_levels=c("Normal", "Tumor_Stage_I"),
                 posclass="Tumor_Stage_I",
                 return_conf_mx=T)
print(conf_mx)

