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


phome <- "/media/anna/bigdisk/Anna/Enzimhome_local/Projects/StatisticsCore/Ras"
#phome <- "/bigdisk/users/alovrics/Projects/StatisticsCore/Ras"
datafolder <- file.path(phome, "Data")
figfolder <- file.path(phome, "Figures")
if (!file.exists(figfolder)) dir.create(figfolder)
calcfolder <- file.path(phome, "Data_calculated")
if (!file.exists(calcfolder)) dir.create(calcfolder)
resfolder <- file.path(phome, "Results")
if (!file.exists(resfolder)) dir.create(resfolder)

######################
for (rtype in c("BRCA", "COAD", "LUAD", "LUSC")){
  # original data
  # read data
  ras_list <- readRDS(file=file.path(calcfolder, rtype,
                                     paste("TCGA_", rtype, ".rds", sep="")))
  ras_data <- ras_list[["data"]]
  indep_vars <- ras_list[["indep_vars"]]
  
  ######################x
  # read obtained results
  resfile <- file.path(resfolder, 
                       paste(rtype,"best_three_gene_models.rds", sep=""))
  reslist <- readRDS(file=resfile)
    # list("mod_list_aic"=mod_list_aic,
  #            "aic_vals"=aic_vals,
  #            "mod_list_acc"=mod_list_acc,
  #            "acc_vals"=acc_vals)
  
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
  figfile <- file.path(figfolder, rtype, "best_model_aic_bic_acc_counts.pdf")
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
  
  print(rtype)
  print(summ_df)
  
#   # is the accuracy better for 3-gene than for 2-gene models?
#   test_genes <- c("KRAS", "ARHGAP20", "PIK3C2G")
#   it_model  <-  glm(as.formula(paste0("sample_type_factor ~ ",
#                                       paste(test_genes,collapse="+"))),
#                     data = ras_data, family = "binomial")
#   it_aic <- it_model$aic
#   it_acc <- runCV(ras_data, it_model, verbose=F)
#   print(c(it_aic, it_acc, BIC(it_model)))
#   
#   # try to add one more gene#   it_aic <- it_model$aic
  #   it_acc <- runCV(ras_data, it_model, verbose=F)
#   ras_data$sample_type_bool <- ifelse(ras_data$sample_type2=="Tumor", F, T)
#   log_min <-  glm(as.formula(paste0("sample_type_bool ~ ", 
#                                     paste(test_genes,collapse="+"))), 
#                   data = ras_data, family = "binomial")
#   log_full <- glm(as.formula(paste0("sample_type_bool ~ ", 
#                                     paste(indep_vars,collapse="+"))), 
#                   data = ras_data, family = "binomial")
#   # Variance Inflation Factor (VIF): High VIF values indicate multicollinearity among predictors, 
#   # which can cause convergence problems.
#   print("variance inflation factor is huge for the full model variables")
#   print(car::vif(log_full))
#   
#   log_best_list <- list()
#   for (omethod in c("AIC", "BIC")){
#     if (omethod=="AIC") kval <- 2
#     if (omethod=="BIC") kval <- log(nrow(ras_data))
#     suppressWarnings(log_best_list[[omethod]] <- log_min %>%
#                        MASS::stepAIC(
#                          scope = list(upper = log_full, lower = ~1),
#                          direction="forward",
#                          k=kval, 
#                          trace=1))
#     print(omethod)
#     print(car::vif(log_best_list[[omethod]])) # still huge vif values
#   }
#   
#   
#   # try with one more gene
#   test_genes <- c("KRAS", "ARHGAP20", "PIK3C2G", "MYO9B")
#   it_model  <-  glm(as.formula(paste0("sample_type_factor ~ ",
#                                       paste(test_genes,collapse="+"))),
#                     data = ras_data, family = "binomial")
#   it_aic <- it_model$aic
#   it_acc <- runCV(ras_data, it_model, verbose=F)
#   print(c(it_aic, it_acc, BIC(it_model)))
}

