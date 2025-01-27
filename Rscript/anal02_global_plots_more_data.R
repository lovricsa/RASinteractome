# 2024.10.01.
# Anna Lovrics

library(readxl) # another package to read xls files
library(dplyr)
library(tidyr) # to transform data for ggplot2 plots
library(data.table) # while I use tibble mostly, fread is best input option
library(tibble)
library(ggplot2) # for visualization
library(ggforce) # for multi-page pdfs
library(EnvStats) # for statistic tests
library(ggcorrplot) # for nice correlation plots
library(cowplot) # to plot q-q plots in a grid
library(gtools) # specifically to generate stars from significance values
library(ggsignif) # to add significance asterics to the plots
library(caret) # for Classification And REgression Training
library(pROC) # to easily plot ROC curves
library(limma) # for the strsplit2 function
library(gridExtra) # to plot a list of figures in a grid
library(ggpubr) # for the same purpose to arrange figures in a grid
library(plotly) # for 3D plot

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

file_list <- c("TCGA_BRCA_breast_new_20241014.tsv",
                  "TCGA_COAD_new_20241014.tsv",
                  "TCGA_LUAD_new_20241014.tsv",
                  "TCGA_LUSC_new_20241014.tsv",
                  "TCGA_PAAD pancreatic_new_20241014.tsv")

type_list <- c("Breast", "COAD", "LUAD", "LUSC", "PAAD")
ras_data_list <- list()
for (ii in 1:length(file_list)){
  infile <- file_list[ii]
  ctype <- type_list[ii]
  ras_data <- as_tibble(fread(file.path(datafolder, infile), sep="\t"))
  ras_data$samples <- NULL
  ras_data$sample_type2 <- "Tumor"
  ras_data$sample_type2[ras_data$sample_type=="Solid Tissue Normal"] <- "Normal"
  ras_data_list[[ctype]] <- ras_data
}
# check whether all columns are the same
for (ii in 1:(length(ras_data_list)-1)){
  for (jj in (ii+1):length(ras_data_list)){
    print(c(ii,jj))
    print(c(names(ras_data_list)[ii], names(ras_data_list)[jj]))
    print(all.equal(colnames(ras_data_list[[ii]]), 
                    colnames(ras_data_list[[jj]])))
    print(setdiff(colnames(ras_data_list[[ii]]), colnames(ras_data_list[[jj]])))
    print(setdiff(colnames(ras_data_list[[jj]]), colnames(ras_data_list[[ii]])))
  }
}

# bind the wide tibbles
corder <- colnames(ras_data_list[[1]])
ras_data <- ras_data_list[[1]]
ras_data$cancer_type <- names(ras_data_list)[1]
for (di in 2:length(ras_data_list)){
  add_data <- ras_data_list[[di]][,match(corder, colnames(ras_data_list[[di]]))]
  add_data$cancer_type <- names(ras_data_list)[di]
  ras_data <- bind_rows(ras_data, add_data)
}

# save ras data
ras_file <- file.path(calcfolder, "ras_all.rds")
saveRDS(ras_data, file=ras_file)

print(xtable::xtable(table(ras_data$cancer_type, ras_data$sample_type2)))
# remove outliers - first only take into account normal or tumor, not the cancer type


# add cancer type & create long data
ras_long_list <- list()
for (ctype in names(ras_data_list)){
  #ras_data_list[[ctype]]$cancer_type <- ctype
  ras_data_selected <- ras_data[ras_data$cancer_type==ctype,]
  # create long data for ggplot visualizations
  ras_long <- ras_data_selected  %>% 
    pivot_longer(!c(sample, sample_type, sample_type2, cancer_type), 
                 names_to = "gene", values_to = "expression")
  # # to keep variables order (use the order in the original file)
  # print(ctype)
  # print(length(interaction(
  #   colnames(ras_data_list[["COAD"]])[3:(ncol(ras_data_list[["COAD"]])-2)],
  #   unique(ras_long$gene)
  # )))
  # ras_long$gene <- factor(ras_long$gene, 
  #                         levels=colnames(ras_data_list[["COAD"]])[
  #                           3:(ncol(ras_data_list[["COAD"]])-2)])
  ras_long$gene <- factor(ras_long$gene, 
            levels=colnames(ras_data)[3:(ncol(ras_data)-2)])
  ras_long_list[[ctype]] <- ras_long
}

# bind all together
ras_long <- do.call(rbind, ras_long_list)



# plot boxplots for all independent variables
# also add result of t-test
pvalue_list <- list()
signif_list <- list()
for (ivar in unique(ras_long$gene)){
  ttest_res <- (
    t.test(ras_data[[ivar]][ras_data$sample_type2=="Normal"],
           ras_data[[ivar]][ras_data$sample_type2=="Tumor"], 
           na.action = "na.omit"))
  pvalue_list[[ivar]] <- ttest_res$p.value
  signif_list[[ivar]] <- stars.pval(ttest_res$p.value)
}
best_indep <- names(which(pvalue_list==min(unlist(pvalue_list))))
signif_tbl <- as_tibble(signif_list)
signif_long <- pivot_longer(signif_tbl, cols = everything(),
                            names_to="gene", values_to="significance")

# to keep variables order:
signif_long$gene <- factor(signif_long$gene, 
                        levels=colnames(ras_data)[3:(ncol(ras_data)-2)])
signif_long$sample_type2 <- "Tumor" # so that asterix is above tumor
signif_long$yval <- 0
for (ivar in unique(ras_long$gene)){
  maxval <- (median(ras_data[[ivar]], na.rm=T) + 
                1.5*IQR(ras_data[[ivar]], na.rm=T))
  signif_long[signif_long$gene==ivar, "yval"] <- 1.1*maxval
}

# plot finally
p <- ggplot(data=ras_long, aes(x=sample_type2, y=expression)) +
  geom_boxplot() + 
  # geom_text(data=signif_long, aes(label=significance, y=yval),
  #           fontface=2) +
  # geom_signif(comparisons = list(c("Normal", "Tumor")),
  #   map_signif_level = TRUE)+ # only works without facets
  facet_wrap(~gene, scale="free")+
  xlab("")
figfile <- file.path(figfolder, "boxplot_all_5samples.pdf")
ggsave(figfile, p, width=14, height=14)


# also plot seperately for each cancer type
ras_long_pages <- list()
gnames <- unique(ras_long$gene)
for (pi in 1:6){
  fi <- 10*(pi-1)+1
  li <- min(10*pi, length(gnames))
  genes <- gnames[fi:li]
  ras_long_pages[[pi]] <- ras_long[ras_long$gene %in% genes,]
  p <- ggplot(data=ras_long_pages[[pi]], aes(x=sample_type2, y=expression)) +
    geom_boxplot() + 
    facet_grid(rows=vars(gene), cols=vars(cancer_type),
               scale="free")
  figfile <- file.path(figfolder, paste("boxplot_all_5samples_grid",
                  as.character(pi), ".pdf", sep=""))
  ggsave(figfile, plot=p, width=7, height=14)
}

# plot seperately but more concisely
p <- ggplot(data=ras_long, 
      aes(x=sample_type2, y=expression, fill=cancer_type)) +
  geom_boxplot() + 
  facet_wrap(~gene, scale="free")+
  xlab("")
figfile <- file.path(figfolder, "boxplot_all_5samples_grouped.pdf")
ggsave(figfile, p, width=14, height=14)

# plot histograms for all independent variables
ras_long_tmp <- ras_long
ras_long_tmp$sample_type2 <- factor(ras_long_tmp$sample_type2,
                                       levels=c("Tumor", "Normal"))
p <- ggplot(data=ras_long_tmp, 
            mapping=aes(x=expression, fill=sample_type2))+
  geom_histogram()+
  facet_wrap(~gene, scales="free")
# plot(p)
figfile <- file.path(figfolder, "histogram_all_5samples.pdf")
ggsave(figfile, p, width=14, height=14)

# also plot seperately for each cancer type
for (pi in 1:6){
  ras_long_tmp <- ras_long_pages[[pi]]
  ras_long_tmp$sample_type2 <- factor(ras_long_tmp$sample_type2,
                                      levels=c("Tumor", "Normal"))
  p <- ggplot(data=ras_long_tmp, 
              mapping=aes(x=expression, fill=sample_type2))+
    geom_histogram()+
    # facet_grid(rows=vars(gene), cols=vars(cancer_type),
    #            scales="free")
    facet_wrap(~gene+cancer_type, scales = "free", ncol=5)
  figfile <- file.path(figfolder, paste("histogram_all_5samples_grid",
                                        as.character(pi), ".pdf", sep=""))
  ggsave(figfile, p, width=7, height=14)
}

# plot correlation among the variables
# do it for all samples, tumor samples and normal samples
# order always as would be for normal
indep_vars <- unique(ras_long$gene)
corr <- round(cor(ras_data[which(ras_data$sample_type2=="Normal"), indep_vars], 
                  use="pairwise.complete.obs",
                  method = "pearson"), digits=1)
# order by the "Normal samples clustering
hc <- hclust(dist(corr), method="ward.D2")
# now plot
for (stype in c("Normal","All", "Tumor")){
  sel_samples <- 1:nrow(ras_data)
  if (stype=="Normal"){
    sel_samples <- which(ras_data$sample_type2=="Normal")
  }
  if (stype=="Tumor"){
    sel_samples <- which(ras_data$sample_type2=="Tumor")
    # randomly select 10% of samples
    sel_samples <- sample(sel_samples, 
          size=round(length(sel_samples)/10))
  }
  if (stype=="All"){
    sel_samples <- caret::createDataPartition(factor(ras_data$sample_type2), 
                                              p=1/11)[[1]]
  }
  corr <- round(cor(ras_data[sel_samples, indep_vars], 
                    use="pairwise.complete.obs",
                    method = "pearson"), digits=1)
  figfile <- file.path(figfolder, 
                paste("correlation_", stype, "_orig_order_5samples.pdf", sep=""))
  p <- ggcorrplot(corr, lab=T,
                  hc.order=F, hc.method="ward.D2")
  ggsave(figfile, plot=p, width=20, height =20)
  corr <- corr[hc$order, hc$order]
  figfile <- file.path(figfolder, 
                paste("correlation_", stype, "_hc_order_5samples.pdf", sep=""))
  p <- ggcorrplot(corr, lab=T,
                  hc.order=F, hc.method="ward.D2")
  ggsave(figfile, plot=p, width=20, height =20)
}


############################################################################
# now remove outliers based on PCA analysis
# scaling does not do anything, as input data is already centered and scaled

# also for indep variables create a numeric value for the cancer type
ras_data$cancer_type2 <- as.numeric(factor(ras_data$cancer_type))
indata_mx <- as.matrix(ras_data[, c("cancer_type2", as.character(indep_vars))])
rownames(indata_mx) <- ras_data$sample
pca <- prcomp(indata_mx, center=T, scale=T,
              tol = sqrt(.Machine$double.eps), rank.=50)
# scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

sumvar <- 0
maxi <- 1
for (pcind in 1:length(pca.var)){
  sumvar <- sumvar + pca.var.per[pcind]
  #print(c(pcind, sumvar))
  if (sumvar > 95){
    # if (sumvar >= 100){
    maxi <- pcind
    #print(maxi)
    break
  }
}
vardf <- as.data.frame(cbind('PC'=1:maxi, 'variance'=pca.var.per[1:maxi]))

p <- ggplot(data=vardf) +
  geom_col(mapping=aes(x=PC, y=cumsum(variance)))+
  xlab("Principal Component") +  ylab("Cumulative percent variation") +
  geom_hline(yintercept=95, linetype="dashed",
             color = "red", linewidth=1)
plot(p)
# save scree plot as well
figfile <- file.path(figfolder, 
                     "ras_all_5samples_scree_plot.pdf")
ggsave(figfile, plot=p, height=5, width=5)
# now plot a PCA plot with annotation
pca_res <- as_tibble(pca$x, rownames = "sample")
pca_annot <- inner_join(pca_res, ras_data, by=c("sample"))

# OK, now do the outlier test
# check outliers for numeric parameters - only first 6
vardf <- vardf[1:6,]
plist <- list()
outlier_list <- list()
#for (ci in bfi:bli){
for (ctype in c("All", unique(ras_data$cancer_type))){
  outlier_list[[ctype]] <- list()
  for (stype in c("normal", "tumor")){
    plist[[stype]] <- list()
    outlier_list[[ctype]][[stype]] <- c()
    for (ci in 1:nrow(vardf)){
      # assume 5 outliers, test if it's true
      if (ctype=="All"){
        normal_samples <- pca_annot[pca_annot$sample_type2=="Normal", c(1,ci+1)]
        tumor_samples <- pca_annot[pca_annot$sample_type2=="Tumor", c(1,ci+1)]
      } else {
        normal_samples <- pca_annot[pca_annot$sample_type2=="Normal" &
                                      pca_annot$cancer_type==ctype, c(1,ci+1)]
        tumor_samples <- pca_annot[pca_annot$sample_type2=="Tumor" &
                                     pca_annot$cancer_type==ctype, c(1,ci+1)]
      }
      pca_samples <- list("normal"=normal_samples,
                          "tumor"=tumor_samples)
      if (ctype=="PAAD" & stype=="normal"){
        rosner_result <- rosnerTest(pca_samples[[stype]][,2,drop=T], k=1)
      } else {
        rosner_result <- rosnerTest(pca_samples[[stype]][,2,drop=T], k=5)
      }
      otl <- rosner_result$n.outliers
      # draw quantile plot, even if there are no outliers
      #if (otl > 0){
      varname <- colnames(pca_annot)[ci+1]
      plotvector <- pca_samples[[stype]][,2,drop=T]
      names(plotvector) <- pca_samples[[stype]][,1,drop=T]
      outlier_ids <- car::qqPlot(plotvector,
                                 ylab=varname, 
                                 id=list(method="y", n=otl, 
                                         cex=1, col=car::carPalette()[1], location="lr"))
      plist[[stype]][[varname]] <- recordPlot()
      outlier_list[[ctype]][[stype]] <- c(outlier_list[[stype]], outlier_ids)
      # print(varname)
      # print(outlier_ids)
      # print("=========")
      #}
      
    }
    figfile = file.path(figfolder, 
            paste("outliers_", stype, "_", ctype, ".pdf", sep=""))
    colfigs <- ceiling(sqrt(length(plist[[stype]])))
    p <- plot_grid(plotlist = plist[[stype]], 
                   ncol=colfigs,
                   align = 'vh', scale=0.8)
    ggsave(figfile, plot=p,  width=7*colfigs, height=5*colfigs)
  }
} # for each cancer type


outlier_all <- unique(c(names(outlier_list[["All"]][["normal"]]), 
                 names(outlier_list[["All"]][["tumor"]])))

# or alternatively
outlier_all <- names(unlist(outlier_list))
outlier_all <- limma::strsplit2(outlier_all, split="\\.")[,3]

pca_annot$outlier <- "NO"
pca_annot$outlier[which(pca_annot$sample %in% outlier_all)] <- "YES"
# show the outliers on the PCA plot

nicelabel <- list()
for (ii in 1:nrow(vardf)){
  nicelabel[[ii]] <- paste('PC',
                           as.character(ii), ' (', as.character(vardf$variance[ii]), '%)', sep="")
}

# plot the results of the PCA analysis
plist <- list()
for (ii in 1:(nrow(vardf)/2)){
  xname <- colnames(pca_annot)[2*ii]
  yname <- colnames(pca_annot)[2*ii+1]
  plist[[ii]] <- ggplot(data=pca_annot) +
    geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                           shape  = sample_type2, col=outlier),
               size = 1, stroke = 1, fill="white")+
    # color points as required
    scale_color_manual("outlier", values = c("NO" = "black", "YES" = "red"))+
    scale_shape_manual("sample type", values = c("Tumor" = 19, "Normal" = 22))+
    labs(x = nicelabel[[2*(ii-1)+1]], y = nicelabel[[2*ii]])+
    #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
    theme_gray(base_size = 20) # change the default font size in the grey theme  
}

pfinal <- plot_grid(plotlist=plist, ncol = 1)
outfig <- file.path(figfolder, "ras_all_5samples_pca_outliers.pdf")
ggsave(outfig, plot=pfinal, width=14, height=28)

# plot LDA with and without outliers
indata_mx <-  as.matrix(ras_data[, c("cancer_type2", as.character(indep_vars))])
ras_data$sample_type4 <- factor(ras_data$sample_type2)
lda_res <- MASS::lda(indata_mx, grouping=ras_data$sample_type4)
lda_values <- predict(lda_res)
all.equal(lda_values$class, ras_data$sample_type4) # "53 string mismatches"
lda_annot <- as_tibble(lda_values$x)
lda_annot$samples <- ras_data$sample
lda_annot$sample_type2 <- factor(ras_data$sample_type2,
                                 levels=c("Tumor", "Normal"))

p <- ggplot(data=lda_annot, mapping=aes(x=LD1, color=sample_type2))+
  geom_histogram(fill="white", position="dodge")
#facet_wrap(~sample_type2, nrow = 2)
outfig <- file.path(figfolder, "ras_all_5samples_lda_hist.pdf")
ggsave(outfig, plot=p, width=4, height=4)

# repeat lda without outliers
outlier_index <- which(ras_data$sample %in% outlier_all)
ras_omit <- ras_data[-outlier_index,]
indata_mx <- as.matrix(ras_omit[, c("cancer_type2", as.character(indep_vars))])
lda_res <- MASS::lda(indata_mx, grouping=ras_omit$sample_type4)
lda_values <- predict(lda_res)
all.equal(lda_values$class, ras_omit$sample_type4) # "20 string mismatches"
lda_annot <- as_tibble(lda_values$x)
lda_annot$samples <- ras_omit$sample
lda_annot$sample_type2 <- factor(ras_omit$sample_type2,
                                 levels=c("Tumor", "Normal"))

p <- ggplot(data=lda_annot, mapping=aes(x=LD1, color=sample_type2))+
  geom_histogram(fill="white", position="dodge")
#facet_wrap(~sample_type2, nrow = 2)
outfig <- file.path(figfolder, "ras_all_5samples_lda_hist_omit.pdf")
ggsave(outfig, plot=p, width=4, height=4)

# verdict: no need to omit outliers

# plot pca showing different cancer types
plist <- list()
for (ii in 1:(nrow(vardf)/2)){
  xname <- colnames(pca_annot)[2*ii]
  yname <- colnames(pca_annot)[2*ii+1]
  plist[[ii]] <- ggplot(data=pca_annot) +
    geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                           col= cancer_type, shape=sample_type2),
               size = 1, stroke = 1, fill="white")+
    # color points as required
    scale_shape_manual("sample type", values = c("Normal" = 1, "Tumor" = 4))+
    scale_color_manual("cancer type", values = c("Breast" = "black", 
                                                 "COAD" = "blue", "LUAD" = "cyan", 
                                                 "LUSC" = "orange", "PAAD" = "red"))+
    labs(x = nicelabel[[2*(ii-1)+1]], y = nicelabel[[2*ii]])+
    #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
    theme_gray(base_size = 20) # change the default font size in the grey theme  
}

pfinal <- plot_grid(plotlist=plist, ncol = 1)
outfig <- file.path(figfolder, "ras_all_5samples_pca_ctypes.pdf")
ggsave(outfig, plot=pfinal, width=14, height=28)  


# repeat for plotting different cancer types seperately
for (ctype in c("All", unique(ras_data$cancer_type))){
  if (ctype=="All"){
    pca_mod <- pca_annot
  } else {
    pca_mod <- pca_annot[pca_annot$cancer_type==ctype,]
  }
  
  plist <- list()
  for (ii in 1:(nrow(vardf)/2)){
    xname <- colnames(pca_annot)[2*ii]
    yname <- colnames(pca_annot)[2*ii+1]
    plist[[ii]] <- ggplot(data=pca_mod) +
      geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                             shape= cancer_type, color=sample_type2),
                 size = 1, stroke = 1, fill="white")+
      # color points as required
      scale_color_manual("sample type", values = c("Normal" = "black", "Tumor" = "red"))+
      scale_shape_manual("cancer type", values = c("Breast" = 0, 
                                      "COAD" = 1, "LUAD" = 2, 
                                      "LUSC" = 3, "PAAD" = 4))+
      labs(x = nicelabel[[2*(ii-1)+1]], y = nicelabel[[2*ii]])+
      #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
      theme_gray(base_size = 20) # change the default font size in the grey theme  
  }
  
  pfinal <- plot_grid(plotlist=plist, ncol = 1)
  outfig <- file.path(figfolder, 
      paste("ras_all_5samples_pca_ctypes_", ctype, ".pdf", sep=""))
  ggsave(outfig, plot=pfinal, width=14, height=28)  
}

#############################################################################

# for the glm create a boolean output
ras_data$sample_type3 <- ifelse(ras_data$sample_type2=="Tumor", F, T)
table(ras_data$cancer_type)
table(ras_data$cancer_type2)

# now the glm using forward selection
# log_min <- glm(as.formula(paste("sample_type3 ~ ", best_indep, sep="")),
#                 data=ras_data, family = "binomial")
log_min <- glm(sample_type3 ~ cancer_type2, 
               data=ras_data, family="binomial")
log_full <- glm(as.formula(paste0("sample_type3 ~ ", 
                paste(c(as.character(indep_vars), "cancer_type2"),collapse="+"))), 
                data = ras_data, family = "binomial")
# Variance Inflation Factor (VIF): High VIF values indicate multicollinearity among predictors, 
# which can cause convergence problems.
car::vif(log_full) # note, very high values!
# keep siginficant variables
keep_vars <- union("cancer_type2", rownames(summary(log_full)$coefficients)[
  which(summary(log_full)$coefficients[,4]<0.05)])
log_keep <- glm(as.formula(paste0("sample_type3 ~ ", 
                     paste(c(as.character(keep_vars)),collapse="+"))), 
                data = ras_data, family = "binomial")
car::vif(log_keep)

# keep 10 lowest
# keep_vars <- union("cancer_type", names(which(car::vif(log_full)<
#                 sort(car::vif(log_full), decreasing = F)[11])))
# log_keep <- glm(as.formula(paste0("sample_type3 ~ ", 
#               paste(c(as.character(keep_vars)),collapse="+"))), 
#                 data = ras_data, family = "binomial")
# car::vif(log_keep)

# # cluster analysis, select one gene from each cluster (indep_vars & cancer_type2)
# # rdata <- ras_data[, c(as.character(indep_vars), "cancer_type2")]
# rdata <- t(ras_data[, as.character(indep_vars)])
# gap_stat <- cluster::clusGap(rdata,
#                     FUN = kmeans,
#                     nstart = 2,
#                     K.max = nrow(rdata)/10,
#                     B = 50,
#                     verbose=T)
# numk <- cluster::maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"],
#                  method="firstmax")
# 
# #plot number of clusters vs. gap statistic
# factoextra::fviz_gap_stat(gap_stat)
# # now perform the clustering
# km <- kmeans(rdata, centers = numk, 
#              iter.max=20, nstart = 25)
# table(km$cluster)
# 
# # plot correlation for each cluster
# cluster_elements  <- rownames(rdata)[which(km$cluster==1)]
# corr <- round(cor(ras_data[, cluster_elements], 
#                   use="pairwise.complete.obs",
#                   method = "pearson"), digits=1)



log_best_list <- list()
for (omethod in c("AIC", "BIC")){
  if (omethod=="AIC") kval <- 2
  if (omethod=="BIC") kval <- log(nrow(ras_data))
  log_best_list[[omethod]] <- log_keep %>%
    MASS::stepAIC(
      scope = list(upper = log_full, lower = log_min),
      direction="both",
      k=kval, 
      trace=0)
  print(car::vif(log_best_list[[omethod]])) # not that bad for BIC
  }

# warnings: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# compar the best models obtained
vars1 <- limma::strsplit2(as.character(log_best_list[[1]]$formula[3]), split=" + ", fixed=T)[1,]
vars2 <- limma::strsplit2(as.character(log_best_list[[2]]$formula[3]), split=" + ", fixed=T)[1,]

# > intersect(vars1, vars2)
# [1] "RIN1"         "RADIL"        "cancer_type2" "PLCE1"        "RASSF3"       "RGS14"        "SNX27"       
# [8] "RIN3"        
# > setdiff(vars1, vars2)
# [1] "KRIT1"   "ARAP1"   "PIK3CB"  "MYO9B"   "APBB1IP" "GRB7"    "MLLT4"   "RAPH1"   "RASSF7"  "RGL2"    "NRAS"   
# [12] "GRB14"   "TIAM2"  
# > setdiff(vars2, vars1)
# [1] "BRAF"   "RAF1"   "PIK3CG"

# test model obtained
log_best_repr <- glm(log_best_list[[2]]$formula,
                     data=ras_data, family = "binomial")
# log_test <- glm(sample_type3 ~ PLCE1 + RGL1 ,
#                 data=ras_data, family = "binomial") # not it converges

################################################################################
# create training and test sets and check prediction power
set.seed(7004) # set seed for reproducibility
# remove missing values

# do a five-fold cross-validation
numCV_outer <- 5
hold_rows_list <- createFolds(
  as.factor(paste(ras_data$sample_type2, ras_data$cancer_type, sep="_")), 
                              k=numCV_outer, list = T)
# ctrlspecs <- trainControl(method="LOOCV", 
#                           savePredictions="all",
#                           classProbs=TRUE)



savefile <- file.path(calcfolder, "best_models.rds")
if (file.exists(savefile)){
  savelist <- readRDS(savefile)
  logreg_model0_list <- savelist[["models"]]
  logreg_model0_resp <- savelist[["responses"]]
} else {
  # repeat for AIC and BIC
  logreg_model0_list <- list()
  logreg_model0_resp <- list()
  for (omethod in c("AIC", "BIC", "full")){
    print(paste("omethod: ", omethod, sep=""))
    logreg_model0_list[[omethod]] <- list()
    logreg_model0_resp[[omethod]] <- factor(rep("Tumor", times=nrow(ras_data)), 
                                            levels=c("Normal", "Tumor"))
    # logreg_model0_prob <- rep(NA, times=nrow(ras_data))
    for (cvi in 1:numCV_outer){
      print(paste("cvi: ", as.character(cvi), sep=""))
      test_rows <- hold_rows_list[[cvi]]
      ras_train <- ras_data[-test_rows,]
      ras_test <- ras_data[test_rows,]
      
      # create the sampling data for the inner groups as well
      inner_train_list <- createFolds(
        as.factor(paste(ras_train$sample_type2, ras_train$cancer_type, sep="_")), 
        k=10, list = T, returnTrain=T)
      ctrlspecs <- trainControl(method="CV",
                                index=inner_train_list, 
                                savePredictions="all",
                                classProbs=TRUE)
      
      # multiple logistic regression the best model
      if (omethod=="full"){
        mformula <- as.formula(paste0("sample_type4 ~ ", 
              paste(c(as.character(indep_vars), "cancer_type2"),collapse="+")))
      } else {
          mformula <- as.formula(paste("sample_type4 ~ ",
                    as.character(log_best_list[[omethod]]$formula[3]), sep="")) 
      }

      logreg_model0 <- train(mformula, 
        data=ras_train,
        method="glm",
        preProcess=NULL,
        family=binomial,
        trControl=ctrlspecs)
      logreg_model0_list[[omethod]][[cvi]] <- logreg_model0
      summary(logreg_model0)
      varImp(logreg_model0)  
      # predict response
      logreg_model0_resp[[omethod]][test_rows] <- predict(logreg_model0, newdata=ras_test)
      # tmp <- predict(logreg_model0, newdata=ras_test)
      # logreg_model0_prob[test_rows] <- predict(logreg_model0, newdata=ras_test,
      #                                          type="prob")$Tumor
      # tmp2 <- predict(logreg_model0, newdata=ras_test, type="prob")$Tumor
    }
  }
  # save the obtained models & responses
  savelist <- list("models"=logreg_model0_list,
                   "responses"=logreg_model0_resp)
  saveRDS(savelist, file=savefile)  
}

for (omethod in names(logreg_model0_list)){
  coef_tbl <- list()
  varimp_tbl <- list()
  for (cvi in 1:numCV_outer){
    coef_tbl[[cvi]] <- as_tibble(
      summary(logreg_model0_list[[omethod]][[cvi]])$coefficients,
      rownames = "variable")
    colnames(coef_tbl[[cvi]]) <- paste("cvi", as.character(cvi), " ", 
                                colnames(coef_tbl[[cvi]]), sep="")
    colnames(coef_tbl[[cvi]])[1] <- "variable"
    
    varimp_tbl[[cvi]] <- as_tibble(
      as.data.frame(varImp(logreg_model0_list[[omethod]][[cvi]])$importance),
                                   rownames="variable")
    colnames(varimp_tbl[[cvi]]) <- paste("cvi", as.character(cvi), " ", 
                                       colnames(varimp_tbl[[cvi]]), sep="")
    colnames(varimp_tbl[[cvi]])[1] <- "variable"
  }
  coef_merged <- Reduce(inner_join, coef_tbl)
  # if (omethod=="AIC"){
  #   print(xtable::xtable(coef_merged))
  # }
  varimp_merged <- Reduce(inner_join, varimp_tbl)
  print(xtable::xtable(varimp_merged))
}

# write results as latex tables - add more thoughts
for (omethod in c("AIC", "BIC")){
  print(xtable::xtable(as.data.frame(varImp(logreg_model0_list[[omethod]], n=50)$importance)))
}


# plot the confusion matrices
plist <- list()
for (omethod in c("AIC", "BIC", "full")){
  conf_mx <- confusionMatrix(data=logreg_model0_resp[[omethod]], 
                             reference=ras_data$sample_type4,
                             positive="Tumor")
  TClass <- factor(c(0, 0, 1, 1))
  PClass <- factor(c(0, 1, 0, 1))
  # y <- c(FN, TP, TN, FP) 
  y <- rep(0, times=4)
  y[1] <- conf_mx$table[2]
  y[2] <- conf_mx$table[1]
  y[3] <- conf_mx$table[4]
  y[4] <- conf_mx$table[3]
  df <- data.frame(TClass, PClass, y)
  
  plist[[omethod]] <- ggplot(data =  df, mapping = aes(x = TClass, y = PClass)) +
    geom_tile(aes(fill = y), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1, size=10) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() + theme(legend.position = "none") 
    # geom_text(aes(x = 1.5, y = 2.5), parse = TRUE, size=2,
    #           label = as.expression(log_best_list[[omethod]]$formula))
  #plot(plist[[omethod]])
}
pfinal <- plot_grid(plotlist=plist, 
                    ncol = 1, labels = names(plist))
outfig <- file.path(figfolder, "conf_mx_5samples_best_models.pdf")
ggsave(outfig, plot=pfinal, width=5, height=14)

# detect incorrectly classified values
icindex <- which(!ras_data$sample_type4 == logreg_model0_resp[["AIC"]])
icsamples <- ras_data$sample[icindex]
isdata <- ras_data[ras_data$sample %in% icsamples,]
print(xtable::xtable(table(isdata$sample_type2, isdata$cancer_type)))

# add to PCA plot the wrongly classified samples
pca_annot$incorrect <- "NO"
pca_annot$incorrect[which(pca_annot$sample %in% icsamples)] <- "YES"
# show the outliers on the PCA plot

# plot the incorrectly classified samples
for (ctype in c("All", unique(ras_data$cancer_type))){
  if (ctype=="All"){
    pca_mod <- pca_annot
  } else {
    pca_mod <- pca_annot[pca_annot$cancer_type==ctype,]
  }
  
  plist <- list()
  for (ii in 1:(nrow(vardf)/2)){
    xname <- colnames(pca_annot)[2*ii]
    yname <- colnames(pca_annot)[2*ii+1]
    plist[[ii]] <- ggplot(data=pca_mod) +
      geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                             shape  = sample_type2, col=incorrect),
                 size = 1, stroke = 1, fill="white")+
      # color points as required
      scale_color_manual("incorrectly classified", values = c("NO" = "black", "YES" = "red"))+
      scale_shape_manual("sample type", values = c("Tumor" = 19, "Normal" = 22))+
      labs(x = nicelabel[[2*(ii-1)+1]], y = nicelabel[[2*ii]])+
      #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
      theme_gray(base_size = 20) # change the default font size in the grey theme  
  }
  
  pfinal <- plot_grid(plotlist=plist, ncol = 1)
  outfig <- file.path(figfolder, 
                      paste("ras_all_5samples_pca_incorrect_classif_ctypes_", ctype, ".pdf", sep=""))
  ggsave(outfig, plot=pfinal, width=14, height=28)  
}

#######################################################

# OK, repeat PCA with the identified variables
indep_vars <- names(log_best_list[["AIC"]]$coefficients)[2:
                       length(names(log_best_list[["AIC"]]$coefficients))]
# scaling does not do anything, as input data is already centered and scaled
indata_mx <- as.matrix(ras_data[, indep_vars])
rownames(indata_mx) <- ras_data$sample
pca <- prcomp(indata_mx, center=T, scale=T,
              tol = sqrt(.Machine$double.eps), rank.=30)
# scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

sumvar <- 0
maxi <- 1
for (pcind in 1:length(pca.var)){
  sumvar <- sumvar + pca.var.per[pcind]
  #print(c(pcind, sumvar))
  if (sumvar > 95){
    # if (sumvar >= 100){
    maxi <- pcind
    #print(maxi)
    break
  }
}
vardf <- as.data.frame(cbind('PC'=1:maxi, 'variance'=pca.var.per[1:maxi]))

p <- ggplot(data=vardf) +
  geom_col(mapping=aes(x=PC, y=cumsum(variance)))+
  xlab("Principal Component") +  ylab("Cumulative percent variation") +
  geom_hline(yintercept=95, linetype="dashed",
             color = "red", linewidth=1)
plot(p)
# save scree plot as well
figfile <- file.path(figfolder, 
                     "ras_all_5samples_scree_plot_AIC.pdf")
ggsave(figfile, plot=p, height=5, width=5)
# now plot a PCA plot with annotation
pca_res <- as_tibble(pca$x, rownames = "sample")
pca_annot <- inner_join(pca_res, ras_data, by=c("sample"))

vardf <- vardf[1:6,]
nicelabel <- list()
for (ii in 1:nrow(vardf)){
  nicelabel[[ii]] <- paste('PC',
                           as.character(ii), ' (', as.character(vardf$variance[ii]), '%)', sep="")
}
# repeat for plotting different cancer types seperately
for (ctype in c("All", unique(ras_data$cancer_type))){
  if (ctype=="All"){
    pca_mod <- pca_annot
  } else {
    pca_mod <- pca_annot[pca_annot$cancer_type==ctype,]
  }
  
  plist <- list()
  for (ii in 1:(nrow(vardf)/2)){
    xname <- colnames(pca_annot)[2*ii]
    yname <- colnames(pca_annot)[2*ii+1]
    plist[[ii]] <- ggplot(data=pca_mod) +
      geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                             shape= cancer_type, color=sample_type2),
                 size = 1, stroke = 1, fill="white")+
      # color points as required
      scale_color_manual("sample type", values = c("Normal" = "black", "Tumor" = "red"))+
      scale_shape_manual("cancer type", values = c("Breast" = 0, 
                                                   "COAD" = 1, "LUAD" = 2, 
                                                   "LUSC" = 3, "PAAD" = 4))+
      labs(x = nicelabel[[2*(ii-1)+1]], y = nicelabel[[2*ii]])+
      #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
      theme_gray(base_size = 20) # change the default font size in the grey theme  
  }
  
  pfinal <- plot_grid(plotlist=plist, ncol = 1)
  outfig <- file.path(figfolder, 
                      paste("ras_all_5samples_pca_ctypes_", ctype, "_AIC.pdf", sep=""))
  ggsave(outfig, plot=pfinal, width=14, height=28)  
}

# add to PCA plot the wrongly classified samples with AIC as well
pca_annot$incorrect <- "NO"
pca_annot$incorrect[which(pca_annot$sample %in% icsamples)] <- "YES"
# show the outliers on the PCA plot

# plot the incorrectly classified samples
for (ctype in c("All", unique(ras_data$cancer_type))){
  if (ctype=="All"){
    pca_mod <- pca_annot
  } else {
    pca_mod <- pca_annot[pca_annot$cancer_type==ctype,]
  }
  
  plist <- list()
  for (ii in 1:(nrow(vardf)/2)){
    xname <- colnames(pca_annot)[2*ii]
    yname <- colnames(pca_annot)[2*ii+1]
    plist[[ii]] <- ggplot(data=pca_mod) +
      geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                             shape  = sample_type2, col=incorrect),
                 size = 1, stroke = 1, fill="white")+
      # color points as required
      scale_color_manual("incorrectly classified", values = c("NO" = "black", "YES" = "red"))+
      scale_shape_manual("sample type", values = c("Tumor" = 19, "Normal" = 22))+
      labs(x = nicelabel[[2*(ii-1)+1]], y = nicelabel[[2*ii]])+
      #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
      theme_gray(base_size = 20) # change the default font size in the grey theme  
  }
  
  pfinal <- plot_grid(plotlist=plist, ncol = 1)
  outfig <- file.path(figfolder, 
                      paste("ras_all_5samples_pca_incorrect_classif_ctypes_", ctype, "_AIC.pdf", sep=""))
  ggsave(outfig, plot=pfinal, width=14, height=28)  
}


# # plot in 3D
# plot_ly(pca_annot,
#         x=~PC1, 
#         y=~PC2, 
#         z=~PC3, 
#         type="scatter3d", mode="markers", 
#         color=~sample_type2)

# perform MDS analysis as well
mds_res <- cmdscale(dist(indata_mx))
colnames(mds_res) <- c("MDS1", "MDS2")
mds_res <- as_tibble(mds_res, rownames = "sample")
mds_annot <- inner_join(mds_res, ras_data, by=c("sample"))

# plot the results of the MDS analysis
p <- ggplot(data=mds_annot) +
  geom_point(mapping=aes(x = MDS1, y = MDS2, col = sample_type2),
             size = 1, stroke = 1, fill="white")+
  # color points as required
  scale_color_manual("sample type", values = c("Tumor" = "red", "Normal" = "blue"))+
  theme_gray(base_size = 20) # change the default font size in the grey theme
plot(p)
resfile <- file.path(figfolder, 
                     "ras_all_5samples_mds_plot_AIC.pdf")
ggsave(resfile, plot=p, width=10, height=6) 

# lda analysis
lda_res <- MASS::lda(indata_mx, grouping=ras_data$sample_type4)
lda_values <- predict(lda_res)
all.equal(lda_values$class, ras_data$sample_type4) # "53 string mismatches"
lda_annot <- as_tibble(lda_values$x)
lda_annot$samples <- ras_data$sample
lda_annot$sample_type2 <- factor(ras_data$sample_type2,
                                 levels=c("Tumor", "Normal"))

p <- ggplot(data=lda_annot, mapping=aes(x=LD1, color=sample_type2))+
  geom_histogram(fill="white", position="dodge")
#facet_wrap(~sample_type2, nrow = 2)
outfig <- file.path(figfolder, "ras_all_5samples_lda_hist_AIC.pdf")
ggsave(outfig, plot=p, width=4, height=4)

# highlight the selected variables on the correlation matrix
omethod <- "BIC"
select_vars <- names(log_best_list[[omethod]]$coefficients)[2:
                  length(names(log_best_list[[omethod]]$coefficients))]
indep_vars <- unique(ras_long$gene)
ras_sel <- ras_data
for (ii in 1:length(unique(ras_long$gene))){
  ci <- match(unique(ras_long$gene)[ii], colnames(ras_sel))
  if (unique(ras_long$gene)[ii] %in% select_vars){
    colnames(ras_sel)[ci] <- paste("S", colnames(ras_sel)[ci], sep="_")
  }
}
indep_vars <- colnames(ras_sel)[3:(ncol(ras_sel)-5)]
corr <- round(cor(ras_sel[which(ras_sel$sample_type2=="Normal"), indep_vars], 
                  use="pairwise.complete.obs",
                  method = "pearson"), digits=1)
# order by the "Normal samples clustering
hc <- hclust(dist(corr), method="ward.D2")
# now plot
for (stype in c("Normal","All", "Tumor")){
  sel_samples <- 1:nrow(ras_sel)
  if (stype=="Normal"){
    sel_samples <- which(ras_sel$sample_type2=="Normal")
  }
  if (stype=="Tumor"){
    sel_samples <- which(ras_sel$sample_type2=="Tumor")
    # randomly select 10% of samples
    sel_samples <- sample(sel_samples, 
                          size=round(length(sel_samples)/10))
  }
  if (stype=="All"){
    sel_samples <- caret::createDataPartition(factor(ras_sel$sample_type2), 
                                              p=1/11)[[1]]
  }
  corr <- round(cor(ras_sel[sel_samples, indep_vars], 
                    use="pairwise.complete.obs",
                    method = "pearson"), digits=1)
  figfile <- file.path(figfolder, 
        paste("correlation_", stype, "_orig_order_5samples_", omethod, ".pdf", sep=""))
  p <- ggcorrplot(corr, lab=T,
                  hc.order=F, hc.method="ward.D2")
  ggsave(figfile, plot=p, width=20, height =20)
  corr <- corr[hc$order, hc$order]
  figfile <- file.path(figfolder, 
            paste("correlation_", stype, "_hc_order_5samples_", omethod, ".pdf", sep=""))
  p <- ggcorrplot(corr, lab=T,
                  hc.order=F, hc.method="ward.D2")
  ggsave(figfile, plot=p, width=20, height =20)
}

############################################################################
# # calculate the confusion matrix based on how well individual genes
# # would seperate tumor from normal
# plist <- list()
# for (vname in indep_vars){
#   sel_gene <- as.symbol(as.character(vname))
#   numCV_outer <- 5
#   hold_rows_list <- createFolds(
#     as.factor(paste(ras_data$sample_type2, ras_data$cancer_type, sep="_")), 
#     k=numCV_outer, list = T)
#   logreg_model0 <- train(
#     eval(bquote(sample_type2 ~ .(sel_gene))),
#     data=ras_train,
#     method="glm",
#     preProcess=NULL,
#     family=binomial,
#     trControl=ctrlspecs)
#   # predict response
#   logreg_model0_resp <- predict(logreg_model0, newdata=ras_test)
#   # confusion matrix from the folds
#   conf_mx <- confusionMatrix(data=logreg_model0_resp,
#                              reference=as.factor(ras_test$sample_type2),
#                              positive="Tumor")
#   TClass <- factor(c(0, 0, 1, 1))
#   PClass <- factor(c(0, 1, 0, 1))
#   # y <- c(FN, TP, TN, FP)
#   y <- rep(0, times=4)
#   y[1] <- conf_mx$table[2]
#   y[2] <- conf_mx$table[1]
#   y[3] <- conf_mx$table[4]
#   y[4] <- conf_mx$table[3]
#   df <- data.frame(TClass, PClass, y)
# 
#   plist[[vname]] <- ggplot(data =  df, mapping = aes(x = TClass, y = PClass)) +
#     geom_tile(aes(fill = y), colour = "white") +
#     geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1, size=10) +
#     scale_fill_gradient(low = "blue", high = "red") +
#     theme_bw() + theme(legend.position = "none")
# }
# 
# n <- length(plist)
# nCol <- floor(sqrt(n))
# pfinal <- plot_grid(plotlist=plist,
#                     ncol = nCol, labels = names(plist))
# outfig <- file.path(figfolder, "conf_mx_single_genes.pdf")
# ggsave(outfig, plot=pfinal, width=14, height=14)
