# 2024.10.01.
# Anna Lovrics

library(readxl) # another package to read xls files
library(dplyr)
library(tidyr) # to transform data for ggplot2 plots
library(data.table) # while I use tibble mostly, fread is best input option
library(tibble)
library(ggplot2) # for visualization
library(EnvStats) # for statistic tests
library(ggcorrplot) # for nice correlation plots
library(cowplot) # to plot q-q plots in a grid
library(gtools) # specifically to generate stars from significance values
library(caret) # for Classification And REgression Training
library(pROC) # to easily plot ROC curves
library(limma) # for the strsplit2 function
library(gridExtra) # to plot a list of figures in a grid
library(ggpubr) # for the same purpose to arrange figures in a grid

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

ras_data <- fread(file.path(datafolder, 
                               "TCGA_COAD_new_20241014.tsv"),
                     sep="\t")
ras_data <- as_tibble(ras_data)
# View(ras_data)
# first two columns are the same -> remove
# all.equal(ras_data$sample, ras_data$samples)
ras_data$samples <- NULL
# add a column for tumor or normal
ras_data$sample_type2 <- ras_data$sample_type
ras_data$sample_type2[ras_data$sample_type2 %in% 
    c("Metastatic", "Primary Tumor", "Recurrent Tumor")] <- "Tumor"
ras_data$sample_type2[ras_data$sample_type2=="Solid Tissue Normal"] <- "Normal"

# check for constant columns and remove if any
cvar <- apply(ras_data[, 3:ncol(ras_data)], 2, var, na.rm=T)
ci <- which(cvar==0)

# create long data for ggplot visualizations
ras_long <- ras_data %>% 
       pivot_longer(!c(sample, sample_type, sample_type2), 
                    names_to = "gene", values_to = "expression")
# to keep variables order:
ras_long$gene <- factor(ras_long$gene, 
                            levels=colnames(ras_data)[3:(ncol(ras_data)-1)])

# plot boxplots for all independent variables
# also add result of t-test
pvalue_list <- list()
signif_list <- list()
for (ivar in unique(ras_long$gene)){
  ttest_res <- (
    t.test(ras_data[[ivar]][ras_data$sample_type2=="Normal"],
           ras_data[[ivar]][ras_data$sample_type2=="Tumor"]))
  pvalue_list[[ivar]] <- ttest_res$p.value
  signif_list[[ivar]] <- stars.pval(ttest_res$p.value)
}
best_indep <- names(which(pvalue_list==min(unlist(pvalue_list))))
signif_tbl <- as_tibble(signif_list)
signif_long <- pivot_longer(signif_tbl, cols = everything(),
                            names_to="gene", values_to="significance")
# to keep variables order:
signif_long$gene <- factor(signif_long$gene, 
                        levels=colnames(ras_data)[3:(ncol(ras_data)-1)])
signif_long$sample_type2 <- "Tumor" # so that asterix is above tumor
signif_long$yval <- 0
for (ivar in unique(ras_long$gene)){
  maxval <- median(ras_data[[ivar]]) + 1.5*IQR(ras_data[[ivar]])
  signif_long[signif_long$gene==ivar, "yval"] <- 1.1*maxval
}

# plot finally
p <- ggplot(data=ras_long, aes(x=sample_type2, y=expression)) +
  geom_boxplot() + 
  geom_text(data=signif_long, aes(label=significance, y=yval),
            fontface=2) +
  facet_wrap(~gene, scale="free")
figfile <- file.path(figfolder, "boxplot_all.pdf")
ggsave(figfile, p, width=14, height=14)

# plot histograms for all independent variables
ras_long_tmp <- ras_long
ras_long_tmp$sample_type2 <- factor(ras_long_tmp$sample_type2,
                                       levels=c("Tumor", "Normal"))
p <- ggplot(data=ras_long_tmp, 
            mapping=aes(x=expression, fill=sample_type2))+
  geom_histogram()+
  facet_wrap(~gene, scales="free")
plot(p)
figfile <- file.path(figfolder, "histogram_all.pdf")
ggsave(figfile, p, width=14, height=14)

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
  }
  corr <- round(cor(ras_data[sel_samples, indep_vars], 
                    use="pairwise.complete.obs",
                    method = "pearson"), digits=1)
  #corr <- corr[hc$order, hc$order]
  figfile <- file.path(figfolder, 
                paste("correlation_", stype, "_orig_order.pdf", sep=""))
  p <- ggcorrplot(corr, lab=T,
                  hc.order=F, hc.method="ward.D2")
  ggsave(figfile, plot=p, width=20, height =20)
}

ras_data$sample_type3 <- ifelse(ras_data$sample_type2=="Tumor", F, T)



# now the glm using forward selection
# log_min <- glm(as.formula(paste("sample_type3 ~ ", best_indep, sep="")),
#                 data=ras_data, family = "binomial")
log_min <- glm(sample_type3 ~ 1, 
               data=ras_data, family="binomial")
log_full <- glm(as.formula(paste0("sample_type3 ~ ", 
                                  paste(indep_vars,collapse="+"))), 
                data = ras_data, family = "binomial")
# Variance Inflation Factor (VIF): High VIF values indicate multicollinearity among predictors, 
# which can cause convergence problems.
car::vif(log_full)

log_best <- log_min %>%
              MASS::stepAIC(
                scope = list(upper = log_full, lower = ~1),
                direction="both",
                k=2, 
                trace=1)
# Warning messages:
#  1: glm.fit: algorithm did not converge 
# 2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# k - the multiple of the number of degrees of freedom used for the penalty
# Only k = 2 gives the genuine AIC
# k = log(n) is sometimes referred to as BIC or SBC.

# test model obtained
log_best_repr <- glm(log_best$formula,
                     data=ras_data, family = "binomial")
log_test <- glm(sample_type3 ~ PLCE1 + RGL1 ,
                data=ras_data, family = "binomial") # not it converges

################################################################################
# create training and test sets and check prediction power
set.seed(7004) # set seed for reproducibility
train_rows <- createDataPartition(as.factor(ras_data$sample_type2), 
                                  p = 0.75, list = FALSE)
ras_train <- ras_data[train_rows,]
ras_test <- ras_data[-train_rows,]
ctrlspecs <- trainControl(method="LOOCV", 
                          savePredictions="all",
                          classProbs=TRUE)
# multiple logistic regression the best model
logreg_model0 <- train(
  as.formula(paste("sample_type2 ~ ",
                   as.character(log_best$formula[3]), sep="")), 
  data=ras_train,
  method="glm",
  preProcess=NULL,
  family=binomial,
  trControl=ctrlspecs)
summary(logreg_model0)
varImp(logreg_model0)
# predict response
logreg_model0_resp <- predict(logreg_model0, newdata=ras_test)
logreg_model0_prob <- predict(logreg_model0, newdata=ras_test, type="prob")
# calculate roc curve
roc_res <- roc(response = ras_test$sample_type,
               predictor = logreg_model0_prob$Tumor)
# plot in ggplot instead
roc_res_plot <- as_tibble(cbind(roc_res$sensitivities,
                                1-roc_res$specificities))
colnames(roc_res_plot) <- c("TPR", "FPR")

# sort by TPR
roc_res_plot <- roc_res_plot[order(roc_res_plot$TPR),]
mylabel <- paste('AUC = ', 
                 as.character(format(round(roc_res$auc, digits=2), nsmall = 2)), 
                 sep="")
# p <- ggplot(data=roc_res_plot) +
#   geom_line(mapping=aes(x = FPR, y = TPR))+
#   geom_text(mapping = aes(x=0.7, y=0.15), label=mylabel,
#             size=8, family="Courier")+
#   ggtitle("sample type ~ Tks4 + CD2AP + CAPZA1 + \n SRC + WASL + GRB2 + CTTN")+
#   # change the default font size in the grey theme
#   theme_gray(base_size = 16)+
#   theme(plot.title = element_text(size = 12, face = "bold", hjust=0.5))
# plot(p)

# confusion matrix from the folds
conf_mx <- confusionMatrix(data=logreg_model0_resp, 
                           reference=as.factor(ras_test$sample_type2),
                           positive="Tumor")
print("########################")
print("multiple logistic regression with best model")
print(conf_mx, digits = 2) 
xtable::xtable(conf_mx$table, digits = 2)
print("########################")

#######################################################
# OK, plot a PCA plot as well
# scaling does not do anything, as input data is already centered and scaled
indata_mx <- ras_data[, indep_vars]
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
  if (sumvar > 99.99){
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
                     "ras_all_scree_plot.pdf")
ggsave(figfile, plot=p, height=5, width=5)
# now plot a PCA plot with annotation
pca_res <- as_tibble(pca$x, rownames = "sample")
pca_annot <- inner_join(pca_res, ras_data, by=c("sample"))

nicelabel <- list()
for (ii in 1:4){
 nicelabel[[ii]] <- paste('PC',
        as.character(ii), ' (', as.character(vardf$variance[ii]), '%)', sep="")
}
nicexlabel <- paste('PC1 (', as.character(vardf$variance[1]), '%)', sep="")
niceylabel <- paste('PC2 (', as.character(vardf$variance[2]), '%)', sep="")

# plot the results of the PCA analysis
p <- ggplot(data=pca_annot) +
  geom_point(mapping=aes(x = PC1, y = PC2, col = sample_type2),
             size = 1, stroke = 1, fill="white")+
  # color points as required
  scale_color_manual("sample type", values = c("Tumor" = "red", "Normal" = "blue"))+
  labs(x = nicelabel[[1]], y = nicelabel[[2]])+
  #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
  theme_gray(base_size = 20) # change the default font size in the grey theme
plot(p)
resfile <- file.path(figfolder, 
                     "ras_all_pca_plot.pdf")
#gsave(resfile, plot=p, width=20, height=12)
ggsave(resfile, plot=p, width=10, height=6)  

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
                     "ras_all_mds_plot.pdf")
ggsave(resfile, plot=p, width=10, height=6) 

# repeat for reduced model
vars_best <- strsplit2(as.character(log_best$formula[3]), split=" \\+ ")[1,,drop=T]
indata_mx <- ras_data[, vars_best]
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
  if (sumvar > 99.99){
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
                     "ras_best_scree_plot.pdf")
ggsave(figfile, plot=p, height=5, width=5)
# now plot a PCA plot with annotation
pca_res <- as_tibble(pca$x, rownames = "sample")
pca_annot <- inner_join(pca_res, ras_data, by=c("sample"))

nicelabel <- list()
for (ii in 1:2){
  nicelabel[[ii]] <- paste('PC',
                           as.character(ii), ' (', as.character(vardf$variance[ii]), '%)', sep="")
}

# plot the results of the PCA analysis
p <- ggplot(data=pca_annot) +
  geom_point(mapping=aes(x = PC1, y = PC2, col = sample_type2),
             size = 1, stroke = 1, fill="white")+
  # color points as required
  scale_color_manual("sample type", values = c("Tumor" = "red", "Normal" = "blue"))+
  labs(x = nicelabel[[1]], y = nicelabel[[2]])+
  #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
  theme_gray(base_size = 20) # change the default font size in the grey theme
plot(p)
resfile <- file.path(figfolder, 
                     "ras_best_pca_plot.pdf")
#gsave(resfile, plot=p, width=20, height=12)
ggsave(resfile, plot=p, width=10, height=6)  

############################################################################
# calculate the confusion matrix based on how well individual genes
# would seperate tumor from normal
plist <- list()
for (vname in indep_vars){
  sel_gene <- as.symbol(as.character(vname))
  logreg_model0 <- train(
    eval(bquote(sample_type2 ~ .(sel_gene))),
    data=ras_train,
    method="glm",
    preProcess=NULL,
    family=binomial,
    trControl=ctrlspecs)
  # predict response
  logreg_model0_resp <- predict(logreg_model0, newdata=ras_test)
  logreg_model0_prob <- predict(logreg_model0, newdata=ras_test, type="prob")
  # confusion matrix from the folds
  conf_mx <- confusionMatrix(data=logreg_model0_resp, 
                             reference=as.factor(ras_test$sample_type2),
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
  
  plist[[vname]] <- ggplot(data =  df, mapping = aes(x = TClass, y = PClass)) +
    geom_tile(aes(fill = y), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1, size=10) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() + theme(legend.position = "none")
}

n <- length(plist)
nCol <- floor(sqrt(n))
pfinal <- plot_grid(plotlist=plist, 
                    ncol = nCol, labels = names(plist))
outfig <- file.path(figfolder, "conf_mx_single_genes.pdf")
ggsave(outfig, plot=pfinal, width=14, height=14)
