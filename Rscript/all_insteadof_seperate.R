## ----setup, include=F----------------------------------------------------------------------------------------
knitr::opts_chunk$set(message=FALSE, warnings=FALSE, fig.width=12, fig.height=5)


## ----libraries-----------------------------------------------------------------------------------------------
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
library(ggpubr) # for the same purpose to arrange figures in a grid and also for significance stars
library(arsenal) # Functions for Large-Scale Statistical Summaries
library(dr4pl) # for ROUT outlier finding
library(gridGraphics) # Functions to convert a page of plots drawn with the 'graphics' package into identical output drawn with the 'grid' package
library(ggpubr) # Add manually p-values to a ggplot
library(rstatix) # Autocompute P-value Positions For Plotting Significance
library(factoextra) # to plot elbow plot for PCA eigenvalues

## ----folders-------------------------------------------------------------------------------------------------
# phome <- "/media/anna/bigdisk/Anna/Enzimhome_local/Projects/StatisticsCore/Ras"
phome <- "/home/anna/Anna/Enzimhome_local/Projects/StatisticsCore/Ras"
datafolder <- file.path(phome, "Data")
figfolder <- file.path(phome, "Figures")
resfolder <- file.path(phome, "Results")
if (!file.exists(figfolder)) dir.create(figfolder)
calcfolder <- file.path(phome, "Data_calculated")
if (!file.exists(calcfolder)) dir.create(calcfolder)
# source content of utility file
rscriptfolder <- file.path(phome, "Rscript")
util_file <- file.path(rscriptfolder, "anal_utils.R")
source(util_file)

# define colors
stage1col <- "#14f2e0"
stage2col <- "#b5d6d6"
stage3col <- "#ceb5b7"
stage4col <- "#F320FA"


# TO DO: define ras_theme

## ----cancer type---------------------------------------------------------------------------------------------
# rtype <- "COAD"
ras_data_list <- list()
ras_scaled_list <- list()
for (rtype in c("COAD", "LUAD", "LUSC")){
  # create folders for each seperate tumor type
  if (!file.exists(file.path(figfolder, rtype))) dir.create(figfolder)
  if (!file.exists(file.path(calcfolder, rtype))) dir.create(calcfolder)
  if (!file.exists(file.path(resfolder, rtype))) dir.create(resfolder)
  
  
## ----read data, results="asis"-------------------------------------------------------------------------------
  ras_data <- fread(file.path(datafolder, 
                              paste("TCGA_", rtype, "_plus_stage_20241028.tsv", sep="")),
                    sep="\t")
  # print dimensions of ras data
  print(rtype)
  print(dim(ras_data))
  ras_data <- as_tibble(ras_data)
  # remove FFPE Scrolls sample
  tmpi <- which(ras_data$sample_type=="FFPE Scrolls")
  if (length(tmpi) > 0) ras_data <- ras_data[-tmpi,]
  # exclude some of the genes
  exclude <- c("RAPGEF3", "RAPGEF5", "RASIP1", "MYO9B", "PIK3C2A", "RGL3", "RGL4",
               "TIAM1", "ARHGAP20","RASSF2", "RASSF3", "RASSF4", "RASSF6", "HRAS", 
               "KRAS", "NRAS")
  keep_cols <- setdiff(colnames(ras_data), exclude)
  ras_data <- ras_data[, keep_cols]
  # View(ras_data)
  # first two columns are the same -> remove
  # all.equal(ras_data$sample, ras_data$samples)
  ras_data$samples <- NULL
  indep_vars <- setdiff(colnames(ras_data), 
                        c( "sample","pathologic_stage","sample_type"))
  # keep original data, replace columns with what we need in the computations
  ras_data_orig <- ras_data
  # add a column for tumor or normal
  ras_data$sample_type2 <- ras_data$sample_type
  ras_data$sample_type2[ras_data$sample_type2 %in% 
                          c("Metastatic", "Primary Tumor", "Recurrent Tumor")] <- "Tumor"
  ras_data$sample_type2[ras_data$sample_type2=="Solid Tissue Normal"] <- "Normal"
  # add a column for stages I-IV or X
  ras_data$pathologic_stage2 <- ras_data$pathologic_stage
  ras_data$pathologic_stage2[ras_data$pathologic_stage2 %in% 
                               c("Stage I","Stage IA", "Stage IB", "Stage IC")] <- "Stage_I"
  
  ras_data$pathologic_stage2[ras_data$pathologic_stage2 %in% 
                               c("Stage II","Stage IIA", "Stage IIB", "Stage IIC")] <- "Stage_II"
  ras_data$pathologic_stage2[ras_data$pathologic_stage2 %in% 
                               c("Stage III","Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage_III"
  ras_data$pathologic_stage2[ras_data$pathologic_stage2 %in% 
                               c("Stage IV", "Stage IVA", "Stage IVB", "Tumor_Stage IVB")] <- "Stage_IV"
  ras_data$pathologic_stage2[ras_data$pathologic_stage2 %in% 
                               c("[Discrepancy]","Stage X")] <- "Stage_X"
  
  # sample type plus stage together
  ras_data$category <- paste(ras_data$sample_type2, ras_data$pathologic_stage2, sep="_")
  ras_data$category2 <- ras_data$category
  ras_data$category2[ras_data$category2 %in% c("Normal_Stage_X", "Normal_Stage_I",
                                               "Normal_Stage_II", "Normal_Stage_III", "Normal_Stage_IV")] <- "Normal"
  
  # check for constant columns and remove if any
  cvar <- apply(ras_data[, indep_vars], 2, var, na.rm=T)
  ci <- which(cvar==0)
  if (length(ci)>0) ras_data <- ras_data[-ci,]
  
  # add some more response variables
  ras_data$sample_type_factor <- factor(ras_data$sample_type2,
                                        levels=c("Normal", "Tumor"))
  ras_data$sample_type_bool <- ifelse(ras_data$sample_type2=="Tumor", F, T)
  ras_data$category_factor <- factor(ras_data$category2)
  
  # remove "Tumor_Stage_X" samples
  print(paste(
      as.character(length(which(ras_data$category_factor=="Tumor_Stage_X"))),
       " number of Tumor_Stage_X samples removed in ", rtype, " data.", sep=""))
  ras_data <- ras_data[which(!ras_data$category_factor=="Tumor_Stage_X"),]
  ras_data$category_factor <- droplevels(ras_data$category_factor)
  
  # save data
  saveRDS(list("data"=ras_data, "indep_vars"=indep_vars),
        file=file.path(calcfolder, rtype, paste("TCGA_", rtype, ".rds", sep="")))
  
  # # print out tables for check
  # print(table(ras_data$sample_type_factor))
  # print(table(ras_data$category_factor))
  # xtable::xtable(table(ras_data$category_factor))
  
  # draw a pie chart
  # TO DO - forgatás, vagy legend: Fraction of total (%)
  data <- as.data.frame(table(ras_data$category_factor))
  sumsamples <- sum(data$Freq)
  data$percent <- 100*(data$Freq/sumsamples)
  data$percent <- paste("(", as.character(round(data$percent,digits=1)), "%)", sep="")
  data$freq_percent <- paste(as.character(data$Freq), data$percent, sep=" ")
  p <- ggplot(data, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(stat="identity", width=1, color="white") +
    geom_text(aes(x=1.2, label = freq_percent),
              position = position_stack(vjust = 0.5),
              color="white", size=8) +
    coord_polar("y", start=0) +
    scale_fill_manual("", 
                       values = c("Normal" = "black",
                                  "Tumor_Stage_I"=stage1col,
                                  "Tumor_Stage_II"=stage2col,
                                  "Tumor_Stage_III"=stage3col,
                                  "Tumor_Stage_IV"=stage4col),
                      labels= c("Normal" = "Normal",
                                "Tumor_Stage_I"="Stage I",
                                "Tumor_Stage_II"="Stage II",
                                "Tumor_Stage_III"="Stage III",
                                "Tumor_Stage_IV"="Stage IV"))+
    theme_void() # remove background, grid, numeric labels
  plot(p)
  
  outfile <- file.path(resfolder, rtype,  "stages_summary.docx")
  if (!file.exists(outfile)){
    tab1 <- tableby( ~ category_factor, data=ras_data)
    # summary(tab1)
    # Create HTML output
    html_output <- summary(tab1, text = FALSE)
    outfile <- file.path(resfolder, rtype, "stages_summary.docx")
    write2word(summary(tab1), outfile)
    outfile <- gsub(".docx", ".html", outfile)
    write2html(summary(tab1), outfile, title="summary")
    #browseURL(outfile)
  }
  ras_data_list[[rtype]] <- ras_data
  print(dim(ras_data))
  
  # also calculate scaled data
  ras_scaled <- ras_data
  ras_scaled[, indep_vars] <- scale(
    ras_scaled[, indep_vars], center=T, scale=T)
  ras_scaled_list[[rtype]] <- ras_scaled
}
  

## ----PCA & Rosner test for outlier detection  ---------------------------------
for (rtype in c("COAD", "LUAD", "LUSC")){
  ras_data <- ras_data_list[[rtype]]
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
               color = "red", linewidth=1)+
    theme(axis.title=element_text(size=20),
  panel.grid.minor=element_blank(),
  panel.grid.major=element_blank(),
  axis.text=element_text(size=16),
  legend.text=element_text(size=16),
  strip.text.x = element_text(size=12), # facet titles
  panel.background=element_rect(fill="white", colour = "black"))
  plot(p)
  # save scree plot as well
  figfile <- file.path(figfolder, rtype,
                       "ras_all_scree_plot.pdf")
  
  # also plot percent variation
  p <- ggplot(data=vardf) +
    geom_col(mapping=aes(x=PC, y=variance))+
    xlab("Principal Component") +  ylab("Percent variation") +
    geom_point(mapping=aes(x=PC, y=variance))+
    geom_line(mapping=aes(x=PC, y=variance))+
    theme(axis.title=element_text(size=20),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.text=element_text(size=16),
          legend.text=element_text(size=16),
          strip.text.x = element_text(size=12), # facet titles
          panel.background=element_rect(fill="white", colour = "black"))
    # geom_hline(yintercept=95, linetype="dashed",
    #            color = "red", linewidth=1)
  plot(p)
  # save scree plot as well
  figfile <- file.path(figfolder, rtype,
                       "ras_all_scree_plot_percent_variation.pdf")
  ggsave(figfile, plot=p, height=5, width=5)
  
  # now plot a PCA plot with annotation
  pca_res <- as_tibble(pca$x, rownames = "sample")
  pca_annot <- inner_join(pca_res, ras_data, by=c("sample"))
  # keep the first 4 PCs only (based on knee plot)
  vardf <- vardf[1:4,]
  nicelabel <- list()
  for (ii in 1:nrow(vardf)){
    nicelabel[[ii]] <- paste('PC',
                       as.character(ii), ' (', as.character(vardf$variance[ii]), '%)', sep="")
  }
  
  plist <- list()
  for (ii in 1:(nrow(vardf)/2)){
    xname <- colnames(pca_annot)[2*ii]
    yname <- colnames(pca_annot)[2*ii+1]
    plist[[ii]] <- ggplot(data=pca_annot) +
      geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                             col  = sample_type2),
                 size = 1, stroke = 1, fill="white")+
      # color points as required
      scale_color_manual("sample type", values = c("Tumor" = "red", "Normal" = "black"))+
      labs(x = nicelabel[[2*(ii-1)+1]], y = nicelabel[[2*ii]])+
      #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
      theme_gray(base_size = 20) + # change the default font size in the grey theme 
    theme(axis.title=element_text(size=20),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.text=element_text(size=16),
          legend.text=element_text(size=16),
          strip.text.x = element_text(size=12), # facet titles
          panel.background=element_rect(fill="white", colour = "black"))
  }
  
  pfinal <- plot_grid(plotlist=plist, ncol = 2)
  plot(pfinal)
  resfile <- file.path(figfolder, rtype,
                       "ras_all_pca_plot.pdf")
  ggsave(resfile, plot=pfinal, width=20, height=12)
  
  # write out labels
  xname <- "PC1"
  yname <- "PC2"
  p_labels <- ggplot(data=pca_annot) +
    geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                           col  = sample_type2),
               size = 1, stroke = 1, fill="white")+
    # color points as required
    scale_color_manual("sample type", values = c("Tumor" = "red", "Normal" = "black"))+
    labs(x = nicelabel[[1]], y = nicelabel[[2]])+
    geom_text(mapping=aes(x = PC1, y = PC2, label=sample), size=1)+
    theme_gray(base_size = 20) + # change the default font size in the grey theme 
  theme(axis.title=element_text(size=20),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        strip.text.x = element_text(size=12), # facet titles
        panel.background=element_rect(fill="white", colour = "black"))
  plot(p_labels)
  resfile <- file.path(figfolder, rtype,
                       "ras_labels_pca_plot.pdf")
  ggsave(resfile, plot=p_labels, width=20, height=12)
  
  # plot the different stages as well
  plist <- list()
  for (ii in 1:(nrow(vardf)/2)){
    xname <- colnames(pca_annot)[2*ii]
    yname <- colnames(pca_annot)[2*ii+1]
    plist[[ii]] <- ggplot(data=pca_annot) +
      geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                             col  = category2),
                 size = 1, stroke = 1, fill="white")+
      # color points as required
      scale_color_manual("category", 
                         values = c("Normal" = "black",
                                    "Tumor_Stage_I"=stage1col,
                                    "Tumor_Stage_II"=stage2col,
                                    "Tumor_Stage_III"=stage3col,
                                    "Tumor_Stage_IV"=stage4col),
                         labels= c("Normal" = "Normal",
                                   "Tumor_Stage_I"="Stage I",
                                   "Tumor_Stage_II"="Stage II",
                                   "Tumor_Stage_III"="Stage III",
                                   "Tumor_Stage_IV"="Stage IV"))+
      labs(x = nicelabel[[2*(ii-1)+1]], y = nicelabel[[2*ii]])+
      #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
      theme_gray(base_size = 20) + # change the default font size in the grey theme 
    theme(axis.title=element_text(size=20),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.text=element_text(size=16),
          legend.text=element_text(size=16),
          strip.text.x = element_text(size=12), # facet titles
          panel.background=element_rect(fill="white", colour = "black"))
  }
  
  pfinal <- plot_grid(plotlist=plist, ncol = 2)
  plot(pfinal)
  resfile <- file.path(figfolder, 
                       "ras_all_pca_categories.pdf")
  ggsave(resfile, plot=pfinal, width=20, height=12)
  
  ## ----find outliers with Rosner test-------------------------------------------------------------------------------------------
  
  plist <- list()
  outlier_list <- list()
  #for (ci in bfi:bli){
  
  outlier_list <- list()
  plist <- list()
  for (ctype in levels(ras_data$category_factor)){
    plist[[ctype]] <- list()
    outlier_list[[ctype]] <- c()
    for (ci in 1:nrow(vardf)){
      # assume 3 outliers, test if it's true
      selected_samples <- pca_annot[pca_annot$category_factor==ctype, c(1,ci+1)]
      rosner_result <- rosnerTest(selected_samples[,2,drop=T], k=3)
      otl <- rosner_result$n.outliers
      # draw quantile plot, even if there are no outliers
      #if (otl > 0){
      varname <- colnames(pca_annot)[ci+1]
      # plotvector <- pca_samples[[stype]][,2,drop=T]
      plotvector <- selected_samples[,2,drop=T]
      names(plotvector) <- selected_samples[,1,drop=T]
      outlier_ids <- car::qqPlot(plotvector,
                                 ylab=varname, 
                                 id=list(method="y", n=otl, 
                                         cex=1, col=car::carPalette()[1], location="lr"))
      plist[[ctype]][[varname]] <- recordPlot()
      outlier_list[[ctype]] <- c(outlier_list[[ctype]], outlier_ids)
      # print(varname)
      # print(outlier_ids)
      # print("=========")
      #}
      
    }
    figfile <- file.path(figfolder, rtype,
                         paste("outliers_", ctype, ".pdf", sep=""))
    colfigs <- ceiling(sqrt(length(plist[[ctype]])))
    p <- plot_grid(plotlist = plist[[ctype]], 
                   ncol=colfigs,
                   align = 'vh', scale=0.8)
    ggsave(figfile, plot=p,  width=7*colfigs, height=5*colfigs)
  }
  
  # list outliers
  outlier_all <- names(unlist(outlier_list))
  outlier_all <- limma::strsplit2(outlier_all, split="\\.")[,2]
  
  pca_annot$outlier <- "NO"
  pca_annot$outlier[which(pca_annot$sample %in% outlier_all)] <- "YES"
  
  # plot PCA with outliers
  plist <- list()
  for (ii in 1:(nrow(vardf)/2)){
    xname <- colnames(pca_annot)[2*ii]
    yname <- colnames(pca_annot)[2*ii+1]
    plist[[ii]] <- ggplot(data=pca_annot) +
      geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                             shape  = sample_type2, col=outlier),
                 size = 1, stroke = 1, fill="white")+
      # color points as required
      scale_color_manual("potential outlier", values = c("NO" = "black", "YES" = "red"))+
      scale_shape_manual("sample type", values = c("Tumor" = 19, "Normal" = 22))+
      labs(x = nicelabel[[2*(ii-1)+1]], y = nicelabel[[2*ii]])+
      #geom_text(mapping=aes(x = PC1, y = PC2, label=samples), size=1)+
      theme(axis.title=element_text(size=20),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            axis.text=element_text(size=16),
            legend.text=element_text(size=16),
            strip.text.x = element_text(size=12), # facet titles
            panel.background=element_rect(fill="white", colour = "black"))
  }
  
  pfinal <- plot_grid(plotlist=plist, ncol = 2)
  plot(pfinal)
  resfile <- file.path(figfolder, rtype,
                       "ras_all_pca_outliers.pdf")
  ggsave(resfile, plot=pfinal, width=20, height=9)
  
  # if we need the outliers
  # pca_annot[pca_annot$outlier=="YES", c("sample", "outlier")]
}

## ----correlation    (Figure 2) -----------------------------------------------
# we need legend (without the work Corr) for each - maradjon 0.2
# Figure 2 - külön ábra legyen a C, D, a tetején ne legyen külön COAD, LUAD, etc.
# uganilyen lost connection ábra: fekete - normál, (zöld különbség), stage I - colstageI,
# stage I -re is (F, G)
#div_colors <- corrplot::COL2('RdBu', 200)

# loop for comparisons
comptype_vector <- c("normal_vs_tumor", "normal_vs_stageI", 
                     "stageI_vs_stageIV")
depvar_vector <- c("sample_type_factor", "category_factor", "category_factor")
depvar_levels_list <- list(c("Normal", "Tumor"),
                           c("Normal", "Tumor_Stage_I"),
                           c("Tumor_Stage_I", "Tumor_Stage_IV"))
depvar_colors <- list("Normal"="black", 
                      "Tumor"="red",
                      "Tumor_Stage_I"=stage1col)

for (ci in c(1:2)){
  comptype <- comptype_vector[ci]
  depvar <- depvar_vector[ci]
  depvar_levels <- depvar_levels_list[[ci]]
  corr_list <- list()
  plot_list <- list()
  for (rtype in c("COAD", "LUAD", "LUSC")){
    ras_data <- ras_data_list[[rtype]]
    set.seed(1004)
    cfigfolder <- file.path(figfolder, rtype, "correlation")
    if (!file.exists(cfigfolder)) dir.create(cfigfolder)
    # plot correlation among the variables
    # do it for all samples, tumor samples and normal samples
    # order always as would be for normal
    corr <- round(cor(ras_data[which(ras_data[[depvar]]=="Normal"), indep_vars], 
                      use="pairwise.complete.obs",
                      method = "spearman"), digits=2)
    # order by the "Normal samples clustering
    hc <- hclust(dist(corr), method="ward.D2")
    # use the number of samples that are minimal in the stages
    slist <- setdiff(unique(ras_data$category2), "Tumor_Stage_X")
    #num_samples_min <- min(table(ras_data[ras_data$category2 %in% slist,]$category2))
    num_samples_min <- 40
    # 
    corr_list[[rtype]] <- list()
    plot_list[[rtype]] <- list()
    for (stype in depvar_levels){
      sel_samples <- which(ras_data[[depvar]]==stype)
      corr <- round(cor(ras_data[sel_samples, indep_vars], 
                        use="pairwise.complete.obs",
                        method = "pearson"), digits=2)
      # 600 dpi TIFF
      p <- ggcorrplot(corr[hc$order, hc$order], method="circle", type="upper",
                      hc.order = F, sig.level=0.1, lab=T, lab_size = 1.2,
                      colors = c("#234387", "white", "#b81007"),
                      show.legend = T,
                      legend.title = "")
      # mean absolute correlation
      corr_list[[rtype]][[stype]] <- unlist(corr[lower.tri(corr, diag = FALSE)])
      plot_list[[rtype]][[stype]] <- p
      # also save each figure separately
      outfig <- file.path(cfigfolder, 
                          paste(rtype, "_correlation_",  stype, ".pdf", sep=""))
      ggsave(outfig, plot=p, height=10, width=9)
      
      # and save once with legend as well
      p <- ggcorrplot(corr[hc$order, hc$order], method="circle", type="upper",
                      hc.order = F, sig.level=0.1, lab=T, lab_size = 1.2,
                      colors = c("#234387", "white", "#b81007"),
                      show.legend = T,
                      legend.title = "")
      outfig <- file.path(cfigfolder, 
                          paste(rtype, "_correlation_",  stype, "_with_legend.pdf", sep=""))
      ggsave(outfig, plot=p, height=10, width=9)
    } # loop for stype
  } # loop for rtype
  
  # now create connection plots
  num_connections <- list()
  num_connections_pos_norm <- list()
  num_connections_pos_tum <- list()
  num_connections_neg_norm <- list()
  num_connections_neg_tum <- list()
  mean_connections <- list()
  mean_connections_diff <- list()
  for (rtype in names(corr_list)){
    for (stype in depvar_levels){
      mean_connections[[rtype]][[stype]] <- mean(
        abs(corr_list[[rtype]][[stype]]))
    }
    mean_connections_diff[[rtype]] <- round(
               mean_connections[[rtype]][[depvar_levels[1]]] -
                    mean_connections[[rtype]][[depvar_levels[2]]], digits=2)
    num_connections[[rtype]] <- list()
    num_connections_pos_norm[[rtype]] <- list()
    num_connections_pos_tum[[rtype]] <- list()
    num_connections_neg_norm[[rtype]] <- list()
    num_connections_neg_tum[[rtype]] <- list()
    for (thr in seq(from=0.05, to=0.9, by=0.05)){
      num_norm <- length(which(
        abs(corr_list[[rtype]][[depvar_levels[1]]])>=thr))
      num_connections_pos_norm[[rtype]][[as.character(thr)]] <- length(which(
        corr_list[[rtype]][[depvar_levels[1]]]>=thr))
      num_connections_neg_norm[[rtype]][[as.character(-thr)]] <- length(which(
        corr_list[[rtype]][[depvar_levels[1]]]<=-thr))
      num_tum <- length(which(
        abs(corr_list[[rtype]][[depvar_levels[2]]])>=thr))
      num_connections_pos_tum[[rtype]][[as.character(thr)]] <- length(which(
        corr_list[[rtype]][[depvar_levels[2]]]>=thr))
      num_connections_neg_tum[[rtype]][[as.character(-thr)]] <- length(which(
        corr_list[[rtype]][[depvar_levels[2]]]<=-thr))
      num_diff <- num_norm - num_tum
      num_connections[[rtype]][[as.character(thr)]] <- num_diff
    }
    num_connections[[rtype]] <- unlist(num_connections[[rtype]])
    num_connections_pos_norm[[rtype]] <- unlist(num_connections_pos_norm[[rtype]])
    num_connections_pos_tum[[rtype]] <- unlist(num_connections_pos_tum[[rtype]])
    num_connections_neg_norm[[rtype]] <- unlist(num_connections_neg_norm[[rtype]])
    num_connections_neg_tum[[rtype]] <- unlist(num_connections_neg_tum[[rtype]])
  } # loop for rtype
  
  # plot number of positive and negative correlations above threshold - 
  # TO DO
  num_list <- list()
  for (rtype in c("COAD", "LUAD", "LUSC")){
    num_list[[rtype]] <- list()
    num_list[[rtype]][[depvar_levels[1]]] <- c(rev(num_connections_neg_norm[[rtype]]),
                                       num_connections_pos_norm[[rtype]])
    num_list[[rtype]][[depvar_levels[2]]] <- c(rev(num_connections_neg_tum[[rtype]]),
                                      num_connections_pos_tum[[rtype]])
    num_list[[rtype]] <- as.data.frame(num_list[[rtype]])
    num_list[[rtype]]$rtype <- rtype
    num_list[[rtype]] <- as_tibble(num_list[[rtype]], rownames = "corr")
    num_list[[rtype]]$corr <- as.numeric(num_list[[rtype]]$corr)
    num_list[[rtype]]$diff <-  (num_list[[rtype]][[depvar_levels[1]]] - 
                                   num_list[[rtype]][[depvar_levels[2]]])
  } # loop for rtype
  
  #rlim <- 0; ylim1 <- 600; ylim2 <- 200
  #rlim <- 0.4; ylim1 <- 200; ylim2 <- 150
  rlim <- 0.2; ylim1 <- 400; ylim2 <- 200
  plist <- list()
  pindex <- 1
  psize <- 6
  linewidth_vector <- rep(0.8, times=6)
  linetype_vector <- c(rep("dashed", times=4), rep("dashed", times=2))
  for (rtype in names(corr_list)){
    plist[[pindex]] <- ggplot(data=num_list[[rtype]])+
      geom_line(mapping=aes(x=corr, y=.data[[depvar_levels[1]]]), 
                color=depvar_colors[[depvar_levels[1]]],
                linewidth=linewidth_vector[pindex], linetype=linetype_vector[pindex])+
      geom_point(mapping=aes(x=corr, y=.data[[depvar_levels[1]]]), , 
                 color=depvar_colors[[depvar_levels[1]]], size=psize)+
      geom_line(mapping=aes(x=corr, y=.data[[depvar_levels[2]]]), , 
                color=depvar_colors[[depvar_levels[2]]],
                linewidth=linewidth_vector[pindex], linetype=linetype_vector[pindex])+
      geom_point(mapping=aes(x=corr, y=.data[[depvar_levels[2]]]), 
                 color=depvar_colors[[depvar_levels[2]]], size=psize)+
      #xlim(-0.8,-rlim)+ 
      ylim(0,ylim1) +
      scale_x_continuous(limits=c(-0.8,-rlim),
                         breaks=seq(-0.8, -rlim, 0.2))+
      ylab("Number of \n connections")+
      xlab("correlation threshold")+
      # theme(axis.text=element_text(size=16),
      #       axis.title=element_text(size=20),
      #       axis.title.x=element_blank(),
      #       axis.text.x=element_blank(),
      #       axis.ticks.x=element_blank())+
      #       # plot.title = element_text(hjust = 0.8, vjust = -10, size = 16),
      #       # plot.title.position = "plot")+
      theme() %+replace%
      la_replace
    #labs(title = rtype)
    #plot(plist[[pindex]])
    
    pindex <- pindex + 1
    plist[[pindex]] <- ggplot(data=num_list[[rtype]])+
      geom_line(mapping=aes(x=corr, y=.data[[depvar_levels[1]]]), 
                color=depvar_colors[[depvar_levels[1]]],
                linewidth=linewidth_vector[pindex], linetype=linetype_vector[pindex])+
      geom_point(mapping=aes(x=corr, y=.data[[depvar_levels[1]]]), 
                 color=depvar_colors[[depvar_levels[1]]], size=psize)+
      geom_line(mapping=aes(x=corr, y=.data[[depvar_levels[2]]]), 
                color=depvar_colors[[depvar_levels[2]]],
                linewidth=linewidth_vector[pindex], linetype=linetype_vector[pindex])+
      geom_point(mapping=aes(x=corr, y=.data[[depvar_levels[2]]]), 
                 color=depvar_colors[[depvar_levels[2]]], size=psize)+
      xlim(rlim,0.8)+ ylim(0,ylim1)+
      xlab("correlation threshold")+
      # theme(axis.text=element_text(size=16),
      #       axis.title=element_text(size=20),
      #   axis.title.x=element_blank(),
      #   axis.text.x=element_blank(),
      #   axis.ticks.x=element_blank(),
      #   axis.title.y=element_blank(),
      #   axis.text.y=element_blank(),
      #   axis.ticks.y=element_blank())+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) %+replace%
      la_replace
    #labs(title = rtype)
    
    pindex <- pindex+1
  }
  p_top <- cowplot::plot_grid(plotlist=plist, nrow=1, 
                              rel_widths = rep(c(1.25,1), times=3))
  plot(p_top)
  # save
  outfig <- file.path(figfolder, 
                      paste("lost_connections_",
                      comptype, "_minus_plus_seperately_top.pdf", sep=""))
  ggsave(outfig, plot=p_top, width=21, height=7)
  
  plist_bottom <- list()
  pindex <- 1
  for (rtype in names(corr_list)){
    plist_bottom[[pindex]] <- ggplot(data=num_list[[rtype]])+
      geom_line(mapping=aes(x=corr, y=diff), color="darkgreen",
                linewidth=linewidth_vector[pindex], linetype=linetype_vector[pindex])+
      geom_point(mapping=aes(x=corr, y=diff), color="darkgreen", size=psize)+
      xlim(-0.8,-rlim)+ ylim(0,ylim2)+
      ylab("Number of \n lost connections")+ xlab("correlation threshold")+
      # theme(axis.text=element_text(size=16),
      #       axis.title=element_text(size=20))
      theme() %+replace%
      la_replace
    pindex <- pindex + 1
    plist_bottom[[pindex]] <- ggplot(data=num_list[[rtype]])+
      geom_line(mapping=aes(x=corr, y=diff), color="darkgreen",
                linewidth=linewidth_vector[pindex], linetype=linetype_vector[pindex])+
      geom_point(mapping=aes(x=corr, y=diff), color="darkgreen", size=psize)+
      xlim(rlim,0.8)+ ylim(0,ylim2)+
      xlab("correlation threshold")+
      # theme(axis.text=element_text(size=16),
      #       axis.title=element_text(size=20),
      #       axis.title.y=element_blank(),
      #       axis.text.y=element_blank(),
      #       axis.ticks.y=element_blank())
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) %+replace%
      la_replace
    pindex <- pindex+1
  }
  p_bottom <- cowplot::plot_grid(plotlist=plist_bottom, nrow=1, 
                                 rel_widths = rep(c(1.25,1), times=3))
  plot(p_bottom)
  # save
  outfig <- file.path(figfolder, paste("lost_connections_",
                              comptype, "_minus_plus_seperately_bottom.pdf", sep=""))
  ggsave(outfig, plot=p_bottom, width=21, height=7)
  
  plist_all <- do.call(c, list(plist, plist_bottom))
  p_all <- cowplot::plot_grid(plotlist=plist_all, nrow=2, 
                              rel_widths = rep(c(1.25,1), times=6))
  plot(p_all)
  outfig <- file.path(figfolder, paste("lost_connections_",
                                       comptype, "_minus_plus_seperately.pdf", sep=""))
  ggsave(outfig, plot=p_all, width=21, height=14)
  
  # write out to file
  outfile <- file.path(resfolder, paste("lost_connections_",
                                    comptype, "_minus_plus_seperately.xlsx", sep=""))
  #openxlsx::write.xlsx((Reduce(bind_rows, num_list)), file=outfile)
  openxlsx::write.xlsx(num_list, file=outfile)
  
  print(comptype)
  print(round(unlist(mean_connections), digits=2))
  print(unlist(mean_connections_diff))
  
}



# repeat for Normal vs type I
corr_list <- list()
plot_list <- list()
for (rtype in c("COAD", "LUAD", "LUSC")){
  ras_data <- ras_data_list[[rtype]]
  set.seed(1004)
  cfigfolder <- file.path(figfolder, rtype, "correlation")
  if (!file.exists(cfigfolder)) dir.create(cfigfolder)
  # plot correlation among the variables
  # do it for all samples, tumor samples and normal samples
  # order always as would be for normal
  corr <- round(cor(ras_data[which(ras_data$sample_type2=="Normal"), indep_vars], 
                    use="pairwise.complete.obs",
                    method = "spearman"), digits=2)
  # order by the "Normal samples clustering
  hc <- hclust(dist(corr), method="ward.D2")
  # use the number of samples that are minimal in the stages
  slist <- setdiff(unique(ras_data$category2), "Tumor_Stage_X")
  # num_samples_min <- min(table(ras_data[ras_data$category2 %in% slist,]$category2))
  num_samples_min <- 40
  # 
  corr_list[[rtype]] <- list()
  plot_list[[rtype]] <- list()
  for (stype in c("Normal", "Tumor_Stage_I")){
    sel_samples <- 1:nrow(ras_data)
    if (stype=="Normal"){
      sel_samples <- which(ras_data$category2=="Normal")
      #sel_samples <- sample(sel_samples, size=num_samples_min)
    }
    if (stype=="Tumor_Stage_I"){
      sel_samples <- which(ras_data$category2=="Tumor_Stage_I")
      #sel_samples <- sample(sel_samples, size=num_samples_min)
    }
    corr <- round(cor(ras_data[sel_samples, indep_vars], 
                      use="pairwise.complete.obs",
                      method = "pearson"), digits=2)
    # p <- ggcorrplot(corr[hc$order, hc$order], method="circle", type="upper",
    #                 hc.order = F, sig.level=0.1, lab=T, lab_size = 1.2)
    p <- ggcorrplot(corr[hc$order, hc$order], method="circle", type="upper",
                    hc.order = F, sig.level=0.1, lab=T, lab_size = 1.2,
                    colors = c("#234387", "white", "#b81007"),
                    show.legend = F,
                    legend.title = "")
    # mean absolute correlation
    corr_list[[rtype]][[stype]] <- unlist(corr[lower.tri(corr, diag = FALSE)])
    plot_list[[rtype]][[stype]] <- p
    # also save each figure separately
    outfig <- file.path(cfigfolder, 
                        paste(rtype, "_correlation_",  stype, ".pdf", sep=""))
    ggsave(outfig, plot=p, height=10, width=9)
    
    # and save once with legend as well
    p <- ggcorrplot(corr[hc$order, hc$order], method="circle", type="upper",
                    hc.order = F, sig.level=0.1, lab=T, lab_size = 1.2,
                    colors = c("#234387", "white", "#b81007"),
                    show.legend = T,
                    legend.title = "")
    outfig <- file.path(cfigfolder, 
                        paste(rtype, "_correlation_",  stype, "with_legend.pdf", sep=""))
    ggsave(outfig, plot=p, height=10, width=9)
  } # loop for stype
} # loop for rtype


# now create connections plots
num_connections <- list()
mean_connections <- list()
for (rtype in names(corr_list)){
  mean_connections[[rtype]][["Normal"]] <- mean(
    abs(corr_list[[rtype]][["Normal"]]))
  mean_connections[[rtype]][["Tumor_Stage_I"]] <- mean(
    abs(corr_list[[rtype]][["Tumor_Stage_I"]]))
  num_connections[[rtype]] <- list()
  for (thr in seq(from=0.2, to=0.9, by=0.05)){
    num_norm <- length(which(
      abs(corr_list[[rtype]][["Normal"]])>=thr))
    num_tum <- length(which(
      abs(corr_list[[rtype]][["Tumor_Stage_I"]])>=thr))
    num_diff <- num_norm - num_tum
    num_connections[[rtype]][[as.character(thr)]] <- num_diff
  }
  num_connections[[rtype]] <- unlist(num_connections[[rtype]])
}
num_connections_df <- as_tibble(as.data.frame(num_connections), rownames="threshold")
num_longer <- pivot_longer(num_connections_df, names_to="cancer_type", 
                           values_to="lost_connections", cols=!threshold)
num_longer$cancer_type <- factor(num_longer$cancer_type)
num_longer$threshold <- as.numeric(num_longer$threshold)
p <- ggplot(data=num_longer, mapping=aes(x=threshold, y=lost_connections,
                                         color = cancer_type))+
  geom_line() + geom_point()+
  guides(color  = guide_legend(position = "inside")) +
  theme(legend.title=element_blank(), # legend title not needed
        legend.justification.top = "left",
        legend.justification.left = "top",
        legend.justification.bottom = "right",
        legend.justification.inside = c(0.95, 0.95),
        legend.location = "plot",
        plot.background=element_rect(color="grey"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=20))+
  ylab("lost connections") + xlab("Pearson correlation threshold")
plot(p)
figfile <- file.path(figfolder, "correlation_lost_connections_normal_vs_stageI.pdf")
ggsave(figfile, plot=p, width=5, height=4)



##-----create violinplots of gene expressions of different stages--------##
# minden plotnál legyen x tengely felirat,
# ne legyen ábra cím (fehér háttér plusz megfelelő színek)
# y tengelyfelirat meg csak egyszer (gén nevek)

# TO DO: violinplots with p-stars
# normal vs tumor (fekete és piros) 
# stages, significance: normal vs stage1 & stage1 vs stage4
# ugyanolyan színek
plotlist <- list()
for (rtype in c("COAD", "LUAD", "LUSC")){
  ras_data <- ras_scaled_list[[rtype]]
  ras_tmp <- ras_scaled_list[[rtype]][
    ras_scaled_list[[rtype]]$category_factor %in% 
      c("Normal", "Tumor_Stage_I", "Tumor_Stage_II", 
        "Tumor_Stage_III", "Tumor_Stage_IV"),]
  ras_tmp$category_factor <- droplevels(ras_tmp$category_factor)
  ras_tmp$category_factor <- factor(ras_tmp$category_factor, levels=
                                      c("Normal", "Tumor_Stage_I", "Tumor_Stage_II", 
                                        "Tumor_Stage_III", "Tumor_Stage_IV"))
  plotlist[[rtype]] <- list()
  for (gene in indep_vars){
    stagelm <- lm(as.formula(paste0(gene, " ~ category_factor")), data=ras_tmp)
    stageav <- aov(stagelm)
    summary(stageav)
    tukeydf <- as_tibble(TukeyHSD(stageav, conf.level=.95)$category_factor, 
                         rownames="comparison")
    tukeydf$p.signif <- stars.pval(tukeydf$`p adj`)
    tukeydf$p.signif[tukeydf$p.signif==" "] <- "ns"
    tukeydf$group1 <- limma::strsplit2(tukeydf$comparison, split="-")[,2]
    tukeydf$group2 <- limma::strsplit2(tukeydf$comparison, split="-")[,1]
    
    p1 <- ggplot(ras_tmp,
                 aes(x=category_factor, y=.data[[gene]], 
                     colour = category_factor)) +
      geom_violin(trim=F) +
      geom_jitter(size=0.5) +
      geom_boxplot(width=0.2)+
      scale_color_manual("category",
                         values = c("Normal" = "black",
                                    "Tumor_Stage_I"="cyan",
                                    "Tumor_Stage_II"="purple",
                                    "Tumor_Stage_III"="pink",
                                    "Tumor_Stage_IV"="red"),
                         labels= c("Normal" = "Normal",
                                   "Tumor_Stage_I"="Stage I",
                                   "Tumor_Stage_II"="Stage II",
                                   "Tumor_Stage_III"="Stage III",
                                   "Tumor_Stage_IV"="Stage IV"))+
      xlab("") + ylab(gene)+ guides(colour="none") + 
      ylim(-6, 6)+
      #ylim(-11, 6)+
      # scale_x_discrete(labels=c("Normal", "Stage I", "Stage II", 
      #                           "Stage III", "Stage IV"))+
      scale_x_discrete(labels=c("", "", "", "", ""))+
      theme(axis.title.x=element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=1),
            plot.margin = unit(c(0.05, 0.05, 0.3, 0.05), "inches"))
    y.position <- 4.5
    # testtbl <- as_tibble(list("group1"="Normal", "group2"="Tumor",
    #                           "p"=stars.pval(tres$p.value), "y.position"=y.position))
    testtbl <- tukeydf[1,]
    colnames(testtbl)[colnames(testtbl)=="p.signif"] <- "p"
    testtbl$y.position <- y.position
    
    p2 <- p1+ stat_pvalue_manual(testtbl)
    plotlist[[rtype]][[gene]] <- p2
    # plot to see warnings
    tt <- tryCatch(plot(p2),error=function(e) e, warning=function(w) w)
    if(is(tt,"warning")){
      print(paste("warning in plot: ", rtype, ", ", gene, sep=""))
    }
  }
  p_all <- cowplot::plot_grid(plotlist=plotlist[[rtype]], ncol=11)
  outfig <- file.path(figfolder, rtype, 
                      paste("tcga_", rtype, "violinplot_stages_with_stars.pdf", sep=""))
  ggsave(outfig, plot=p_all, width=28, height=14)
}


## ----LDA (Figure 3) stages but use Normal and Stages I and IV only for prediction---------------------------------------

for (rtype in c("COAD", "LUAD", "LUSC")){
  ras_data <- ras_scaled_list[[rtype]]
  ri <- which(ras_data$category2 %in% c("Normal", 
                                        "Tumor_Stage_I","Tumor_Stage_IV"))
  # also delete unused factors
  ras_tmp <- ras_data[ri, ]
  ras_tmp$category_factor <- droplevels(ras_tmp$category_factor)
  indata_mx <-  as.matrix(ras_tmp[, indep_vars])
  lda_model <- MASS::lda(indata_mx, grouping=ras_tmp$category_factor)
  lda_values <- predict(lda_model, newdata = ras_tmp[, indep_vars])
  all.equal(lda_values$class, ras_tmp$category_factor) # "3 string mismatches"
  lda_annot <- as_tibble(lda_values$x)
  lda_annot$sample <- ras_tmp$sample
  lda_annot$category_factor <- ras_tmp$category_factor
  
  # obtain the “proportion of trace” which is the proportion of between-class 
  # variance that is explained by successive discriminant functions. 
  prop_trace <- prop.table(lda_model$svd^2)
  nicelabel <- list()
  for (ii in 1:length(prop_trace)){
    nicelabel[[ii]] <- paste('LD',
                             as.character(ii), ' (', 
                             as.character(round(100*prop_trace[ii], digits = 2)), '%)', sep="")
  }
  
  # predict position for stageII & stageIII
  lda_values <- predict(lda_model, newdata = ras_data[, indep_vars])
  lda_annot <- as_tibble(lda_values$x)
  lda_annot$sample <- ras_data$sample
  lda_annot$category_factor <- ras_data$category_factor
  
  xname <- "LD1"
  yname <- "LD2"
  pfinal <- ggplot(data=lda_annot) +
    geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                           color=category_factor),
               size = 3, stroke = 1, fill="white")+
    # color points as required
    scale_color_manual("", 
                       values = c("Normal" = "black",
                                  "Tumor_Stage_I"=stage1col,
                                  "Tumor_Stage_II"=stage2col,
                                  "Tumor_Stage_III"=stage3col,
                                  "Tumor_Stage_IV"=stage4col),
                       labels= c("Normal" = "Normal",
                                 "Tumor_Stage_I"="Stage I",
                                 "Tumor_Stage_II"="Stage II",
                                 "Tumor_Stage_III"="Stage III",
                                 "Tumor_Stage_IV"="Stage IV"))+
    # scale_shape_manual("", guide = "none",
    #                    values = c("Normal" = 19,
    #                               "Tumor_Stage_I"=19,
    #                               "Tumor_Stage_II"=1,
    #                               "Tumor_Stage_III"=1,
    #                               "Tumor_Stage_IV"=19))+
    xlab(nicelabel[[1]]) + ylab(nicelabel[[2]])+
    # theme(plot.background=element_rect(color="white"),
    #            axis.text=element_text(size=16),
    #            axis.title=element_text(size=20),
    #       legend.text=element_text(size=16),
    #       panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(),
    #       panel.background=element_rect(fill="white", colour = "black"))
    theme_la()

  plot(pfinal)
  outfig <- file.path(figfolder, rtype, "ras_lda_normal_stageI_stageIV.pdf")
  ggsave(outfig, plot=pfinal, width=7, height=5)
  outfig <- file.path(figfolder, rtype, "ras_lda_normal_stageI_stageIV.tiff")
  ggsave(outfig, plot=pfinal, width=7, height=7, dpi=600)
  
  # density plots for all stages
  lda_values <- predict(lda_model, newdata = ras_data[, indep_vars])
  lda_annot <- as_tibble(lda_values$x)
  lda_annot$sample <- ras_data$sample
  lda_annot$category_factor <- ras_data$category_factor
  
  # plot LD1 density plot
  lda_tmp <- lda_annot[!lda_annot$category_factor %in% 
                         c("Tumor_Stage_X"),]
  lda_tmp$category_factor <- droplevels(lda_tmp$category_factor)
  
  p <- ggplot(data=lda_tmp, mapping=aes(x=LD1, color=category_factor))+
    geom_density(fill="white",  alpha=0.25, linewidth=1.5)+
    scale_color_manual("",
                       values = c("Normal" = "black",
                                  "Tumor_Stage_I"=stage1col,
                                  "Tumor_Stage_II"=stage2col,
                                  "Tumor_Stage_III"=stage3col,
                                  "Tumor_Stage_IV"=stage4col),
                       labels= c("Normal" = "Normal",
                                 "Tumor_Stage_I"="Stage I",
                                 "Tumor_Stage_II"="Stage II",
                                 "Tumor_Stage_III"="Stage III",
                                 "Tumor_Stage_IV"="Stage IV"))+
    theme(plot.background=element_rect(color="white"),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill="white", colour = "black"))
  
  plot(p)
  outfig <- file.path(figfolder, rtype, "ras_LD1_density_categories.pdf")
  ggsave(outfig, plot=p, width=7, height=5)
  
  # LD2 density plot
  lda_tmp <- lda_annot[!lda_annot$category_factor %in% 
                         c("Tumor_Stage_X"),]
  lda_tmp$category_factor <- droplevels(lda_tmp$category_factor)
  
  p <- ggplot(data=lda_tmp, mapping=aes(x=LD2, color=category_factor))+
    geom_density(fill="white",  alpha=0.25, linewidth=2)+
    scale_color_manual("",
                       values = c("Normal" = "black",
                                  "Tumor_Stage_I"=stage1col,
                                  "Tumor_Stage_II"=stage2col,
                                  "Tumor_Stage_III"=stage3col,
                                  "Tumor_Stage_IV"=stage4col),
                       labels= c("Normal" = "Normal",
                                 "Tumor_Stage_I"="Stage I",
                                 "Tumor_Stage_II"="Stage II",
                                 "Tumor_Stage_III"="Stage III",
                                 "Tumor_Stage_IV"="Stage IV"))+
    theme(plot.background=element_rect(color="white"),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill="white", colour = "black"))
  plot(p)
  outfig <- file.path(figfolder, rtype, "ras_LD2_density_categories.pdf")
  ggsave(outfig, plot=p, width=7, height=5)
  
  # LD2 density plot - repeat without normal
  lda_tmp <- lda_annot[!lda_annot$category_factor %in% 
                         c("Normal", "Tumor_Stage_X"),]
  lda_tmp$category_factor <- droplevels(lda_tmp$category_factor)
  
  p <- ggplot(data=lda_tmp, mapping=aes(x=LD2, color=category_factor))+
    geom_density(fill="white",  alpha=0.25, linewidth=2)+
                # key_glyph = draw_key_path)+
    scale_color_manual("",
                       values = c("Tumor_Stage_I"=stage1col,
                                  "Tumor_Stage_II"=stage2col,
                                  "Tumor_Stage_III"=stage3col,
                                  "Tumor_Stage_IV"=stage4col),
                       labels= c("Tumor_Stage_I"="Stage I",
                                 "Tumor_Stage_II"="Stage II",
                                 "Tumor_Stage_III"="Stage III",
                                 "Tumor_Stage_IV"="Stage IV"))+
    theme(plot.background=element_rect(color="white"),
          axis.text=element_text(size=16),
          axis.title=element_text(size=20),
          legend.text=element_text(size=16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill="white", colour = "black")) 
  plot(p)
  outfig <- file.path(figfolder, rtype, "ras_LD2_density_categories_noNormal.pdf")
  ggsave(outfig, plot=p, width=7, height=5)
}


## ----k gene balanced accuracy accross all comparisons (Figure 4)-----------------------------------------------
# figure 4 barplot - színes betű - ITT TARTOK!

comptype_vector <- c("normal_vs_tumor", "normal_vs_stageI", 
                     "stageI_vs_stageIV")
depvar_vector <- c("sample_type_factor", "category_factor", "category_factor")
depvar_levels_list <- list(c("Normal", "Tumor"),
                           c("Normal", "Tumor_Stage_I"),
                           c("Tumor_Stage_I", "Tumor_Stage_IV"))

# for (rtype in c("COAD", "LUAD", "LUSC")){
#   ras_data <- ras_scaled_list[[rtype]]
# }
rtype <- "COAD"
ras_data <- ras_scaled_list[[rtype]]
bathreshold <- c(0.9, 0.95, 0.65)

# save results when comparing stageI with stageIV
stageIvsIVplots <- list()
for (comparison_index in c(1:3)){
  plotlist <- list()
  measures_list <- list()
  depvar <- depvar_vector[comparison_index]
  depvar_levels <- depvar_levels_list[[comparison_index]]
  ras_tmp <- ras_data[ras_data[[depvar]] %in%  depvar_levels,
                                      c(indep_vars, depvar)]
  ras_tmp[[depvar]] <- droplevels(ras_tmp[[depvar]])
  comptype <- comptype_vector[comparison_index]
  # calculate all measures
  for (k in c(1:2)){
    print(paste("start of calculation ", as.character(k), " gene measures: ",
                comptype, sep=""))
    acc_file <- file.path(resfolder, rtype,
                          paste("k", as.character(k), "_gene_model_measures_",
                                comptype, ".rds", sep=""))
    all_measures <- calculate_k_gene_measures(acc_file,
                                              k,
                                              indep_vars,
                                              ras_tmp,
                                              depvar,
                                              depvar_levels_list[[comparison_index]],
                                              depvar_levels_list[[comparison_index]][2])
    # change '.' to '+'
    all_measures$genes <- gsub("\\.", "+", all_measures$genes)
    all_measures <- all_measures[order(all_measures$balanced_accuracy, decreasing = T),]
    measures_list[[k]] <- all_measures
    
    # plot 
    
    # in the case of k==1: color is based on gene expression change
    if (k==1){
      ras_long <- pivot_longer(ras_tmp, !all_of(depvar),
                               names_to="genes", values_to="scaled_expr")
      ras_long[[depvar]] <- factor(ras_long[[depvar]],
                                   levels=depvar_levels)
      stat.test <- ras_long %>%
        group_by(genes) %>%
        t_test(as.formula(paste0("scaled_expr ~ ", depvar))) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance("p.adj")
      
      # test_ARAF <- aggregate(ras_tmp[, c("ARAF", "sample_type_factor")], 
      #           list(ras_tmp$sample_type_factor), mean)
      # t.test(ras_tmp$ARAF[ras_tmp$sample_type_factor=="Normal"], 
      #        ras_tmp$ARAF[ras_tmp$sample_type_factor=="Tumor"])
      
      # add the result of stat.test of indata
      indata <- all_measures
      indata$balanced_accuracy <- round(indata$balanced_accuracy, digits=2)
      indata <- inner_join(indata, stat.test[, c("genes",  "statistic", "p.adj")],
                           by = join_by(genes))
      # for some reason groups are joined the other way round
      indata$statistic <- -indata$statistic
      indata$barcolor <- "#d4d4d4"
      indata$barcolor[which(indata$p.adj < 0.05 & indata$statistic > 0)] <- "#CD3232"
      indata$barcolor[which(indata$p.adj < 0.05 & indata$statistic < 0)] <- "#33ce33"
      
      # order for plot
      indata <- indata[order(indata$balanced_accuracy, decreasing = F),]
      indata$genes <- factor(indata$genes, levels=indata$genes)
      
      plotlist[[k]] <- ggplot(indata, aes(x = balanced_accuracy, y = genes)) +
        #geom_col(fill = "gray70") +
        geom_col(aes(fill = barcolor), color="black", width = 0.8) +
        scale_fill_manual("category",
                          values = c("#d4d4d4" = "#d4d4d4",
                                     "#CD3232"="#CD3232",
                                     "#33ce33"="#33ce33"))+
        ## add percentage labels
        geom_text(aes(label = genes, x=0.2), nudge_x = -0.05) +
        # geom_vline(xintercept=0.95, color="red")+
        #scale_x_continuous(breaks= seq(0.2, 0.8, by = 0.2)) +
        # add accuracy values
        geom_text(aes(label = '-', x=balanced_accuracy, hjust=0)) +
        geom_text(aes(label = '-', x=balanced_accuracy, hjust=1)) +
        geom_text(aes(label = balanced_accuracy, x=0.01+balanced_accuracy, hjust=0),) +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text=element_text(size=14),
              axis.title=element_text(size=18),
              legend.text=element_text(size=14),
              legend.position = "none",
              panel.grid.major =element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_rect(fill="white", colour = "black"),
              plot.background=element_rect(color="white")
              )+
        scale_x_continuous(limits=c(0,1), breaks= seq(0, 1, by = 0.2))+
        xlab("Balanced accuracy")
    } else { # TO DO for k=2: 0.95 felett; bar-on kívülre: balanced accuracy számmal
      # multicolor text: https://andrewwhitby.com/2017/09/18/multi-color-text-ggplot2/
      # stage1 vs stage4: for k=3-5: first 10 bars
      k1data <- indata
      
      #indata <- indata[1:length(indep_vars), ]
      indata <- all_measures
      if (comptype=="stageI_vs_stageIV"){
        indata <- indata[order(indata$balanced_accuracy, decreasing = T),]
        indata <- indata[1:10,]
      } else {
        indata <- indata[indata$balanced_accuracy>=bathreshold[comparison_index],]
      }
      indata$balanced_accuracy <- round(indata$balanced_accuracy, digits=2)
      indata$gene1 <- limma::strsplit2(indata$genes, split="\\+")[,1]
      indata$gene2 <- limma::strsplit2(indata$genes, split="\\+")[,2]
      # increasing order for annotation
      indata <- indata[order(indata$balanced_accuracy, decreasing = T),]
      indata$genes <- factor(indata$genes, levels=indata$genes)
      
      # add annotation - I could only manage it using a loop
      annot_add <- list()
      for (ii in 1:nrow(indata)){
        var1 <- indata$gene1[ii]
        var2 <- indata$gene2[ii]
        col1 <- k1data$barcolor[match(var1, k1data$genes)]
        if (col1=="#d4d4d4") col1 <- 'black'
        col2 <- k1data$barcolor[match(var2, k1data$genes)]
        if (col2=="#d4d4d4") col2 <- "black"
        jj <- 3*(ii-1)+1
        annot_add[[jj]] <-annotate("text", x=0.02, y=nrow(indata)+1-ii, hjust=0, parse=F,
                                  label=bquote(.(var1)*phantom("+")*phantom(.(indata$gene2[ii]))), color=col1)
        annot_add[[jj+1]] <-annotate("text", x=0.02, y=nrow(indata)+1-ii, hjust=0, parse=F,
                                  label=bquote(phantom(.(var1))*"+"*phantom(.(indata$gene2[ii]))))
        annot_add[[jj+2]] <-annotate("text", x=0.02, y=nrow(indata)+1-ii, hjust=0, parse=F,
                                  label=bquote(phantom(.(var1))*phantom("+")*.(indata$gene2[ii])), color=col2)
      }
      
      # opposite order for plot
      indata <- indata[order(indata$balanced_accuracy, decreasing = F),]
      indata$genes <- factor(indata$genes, levels=indata$genes)
      
      plotlist[[k]] <- ggplot(indata, aes(x = balanced_accuracy, y = genes)) +
        geom_col(fill = "white", color="black", width = 0.8) +
        ## add gene name labels - could not make it work like this
        # geom_text(aes(label=bquote(phantom(.(gene1))*.(gene2)), x=0.02, hjust=0), 
        #           parse=T)+
        annot_add+
        # add accuracy values
        geom_text(aes(label = '-', x=balanced_accuracy, hjust=0)) +
        geom_text(aes(label = '-', x=balanced_accuracy, hjust=1)) +
        geom_text(aes(label = balanced_accuracy, x=0.01+balanced_accuracy, hjust=0),) +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text=element_text(size=14),
              axis.title=element_text(size=18),
              legend.text=element_text(size=14),
              panel.grid.major =element_blank(),
              panel.grid.minor=element_blank(),
              panel.background=element_rect(fill="white", colour = "black"),
              plot.background=element_rect(color="white")
              ) +
        scale_x_continuous(limits=c(0,1.05), breaks= seq(0, 1.2, by = 0.2)) +
        xlab("Balanced accuracy") 
      #plot(plotlist[[k]])
    } # if k==1
    # save stage I vs stage IV plots
    if (comparison_index==3) stageIvsIVplots[[k]] <- plotlist[[k]]
  } # end of k loop
  # save to text file as well
  xls_file <-  file.path(resfolder, 
                         paste("k", "_gene_model_measures_",
                               comptype, ".xls", sep=""))
  openxlsx::write.xlsx(measures_list, file=xls_file)
  # save the figures together
  p <- cowplot::plot_grid(plotlist=plotlist)
  #plot(p)
  figfile <- file.path(figfolder, rtype, paste("k", "_gene_model_measures_",
                                        comptype, ".pdf", sep=""))
  ggsave(figfile, plot=p, width=10, height=7)
}



## ----three and four genes balanced accuracy at different comparisons: results from server----------------------------------------------------------
num_genes_vector <- c(3:5)
comptype_vector <- c("normal_vs_tumor", "normal_vs_stageI", 
                     "stageI_vs_stageIV")

plotlist <- list()
for (num_genes in num_genes_vector){
  for (comptype in comptype_vector){
    acc_file <- file.path(gsub(paste("/", rtype, sep=""), "", resfolder),
                          paste(rtype, "best_", as.character(num_genes), 
                                "gene_models_", comptype,".rds", sep=""))
    if (!file.exists(acc_file)){
      print(paste("file ", acc_file, " does not exist", sep=""))
      next
    } else {
      print(paste(acc_file, " plot is being created", sep=""))
    }
    all_measures_list <- readRDS(acc_file)
    
    # create a dataframe
    all_measures <- as.data.frame(all_measures_list[c("aic_vals",
                                                      "bic_vals", "acc_vals")])
    all_measures$gene_trio <- "NA"
    for (ii in 1:length(all_measures_list[["mod_list"]])){
      coefs <- all_measures_list[["mod_list"]][[ii]]$coefficients
      coefs <- paste0(names(coefs)[2:length(coefs)], collapse="+")
      all_measures$gene_trio[ii] <- coefs
    }
    
    # add individual genes & keep indep vars_only
    keep_rows_list <- list()
    for (gi in 1:num_genes){
      all_measures$newgene <- limma::strsplit2(all_measures$gene_trio, 
                                               split="\\+")[,gi]
      keep_rows_list[[gi]] <- which(all_measures$newgene %in% indep_vars)
      colnames(all_measures)[colnames(all_measures)=="newgene"] <- (
            paste("gene", as.character(gi), sep=""))
    }
    
    # keep indep_vars only
    keep_rows <- Reduce(intersect, keep_rows_list)
    all_measures <- all_measures[keep_rows,]
    
    # save to text file as well
    names(all_measures) <- c("AIC", "BIC", "balanced_accuracy", "gene_trio")
    all_measures <- all_measures[c("gene_trio", "AIC", "BIC", "balanced_accuracy")]
    xls_file <- gsub(".rds", ".xlsx", acc_file)
    all_measures <- all_measures[
      order(all_measures$balanced_accuracy, decreasing = T),]
    openxlsx::write.xlsx(all_measures, file=xls_file)
    
    # plot first 10 hits only
    all_measures <- all_measures[order(all_measures$balanced_accuracy, decreasing = T),]
    all_measures <- all_measures[1:10,]
    # for the plot, reverse order is needed
    all_measures <- all_measures[
      order(all_measures$balanced_accuracy, decreasing = F),]
    all_measures$gene_trio <- factor(all_measures$gene_trio, 
                                     levels=all_measures$gene_trio)
    
    
    # copied
    indata <- all_measures
    indata$balanced_accuracy <- round(indata$balanced_accuracy, digits=2)
    colnames(indata)[colnames(indata)=="gene_trio"] <- "genes"

    plotlist[[num_genes]] <- ggplot(indata, aes(x = balanced_accuracy, y = genes)) +
      geom_col(fill = "white", color="black", width = 0.8) +
      geom_text(aes(label=genes, x=0.02, hjust=0))+
      # add accuracy values
      geom_text(aes(label = '-', x=balanced_accuracy, hjust=0)) +
      geom_text(aes(label = '-', x=balanced_accuracy, hjust=1)) +
      geom_text(aes(label = balanced_accuracy, x=0.01+balanced_accuracy, hjust=0),) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text=element_text(size=14),
            axis.title=element_text(size=18),
            legend.text=element_text(size=14),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background=element_rect(fill="white", colour = "black")
      )+
      scale_x_continuous(limits=c(0,1.05), breaks= seq(0, 1.2, by = 0.2)) +
      xlab("Balanced accuracy") 
    #plot(plotlist[[k]])
    # outfig <- file.path(figfolder, 
    #                     paste(rtype, "best_", as.character(num_genes), 
    #                           "gene_models_", comptype,".pdf", sep=""))
    # ggsave(outfig, plot=plotlist[[k]], width=4, height=14)
  } # if k==1
  # save stage I vs stage IV plots
  if (comptype=="stageI_vs_stageIV"){
    stageIvsIVplots[[num_genes]] <- plotlist[[num_genes]]
  }
  # end of copied
}

# change xlim
for (kk in 1:5){
  stageIvsIVplots[[kk]] <- stageIvsIVplots[[kk]]+
    scale_x_continuous(limits=c(0,1), breaks= seq(0, 1, by = 0.2))
}
# put together
p_top <- cowplot::plot_grid(plotlist=stageIvsIVplots[1:2], nrow=1)
p_bottom <- cowplot::plot_grid(plotlist=stageIvsIVplots[3:5], nrow=1)
p_all <- cowplot::plot_grid(p_top, p_bottom, nrow=2)
plot(p_all)

p_left <- cowplot::plot_grid(plotlist=stageIvsIVplots[1], nrow=1)
p_right <- cowplot::plot_grid(plotlist=stageIvsIVplots[2:5], nrow=2)
p_all <- cowplot::plot_grid(p_left, p_right, nrow=1, 
                            rel_widths = c(1,2))
plot(p_all)
figfile <- file.path(figfolder, rtype, 
                     "k_gene_model_measures_stageI_vs_stageIV.pdf")
ggsave(figfile, plot=p_all, width=21, height=14)

## ----calculate anova for stages------------------------------------------------------------------------------
ras_tmp <- ras_data[ras_data$category2 %in% 
                      c("Normal", "Tumor_Stage_I", "Tumor_Stage_IV"),]
plotlist <- list()
for (gene in indep_vars){
  stagelm <- lm(as.formula(paste0(gene, " ~ category_factor")), data=ras_tmp)
  stageav <- aov(stagelm)
  summary(stageav)
  tukeydf <- as_tibble(TukeyHSD(stageav, conf.level=.95)$category_factor, 
                       rownames="comparison")
  tukeydf$p.signif <- stars.pval(tukeydf$`p adj`)
  tukeydf$p.signif[tukeydf$p.signif==" "] <- "ns"
  tukeydf$group1 <- limma::strsplit2(tukeydf$comparison, split="-")[,2]
  tukeydf$group2 <- limma::strsplit2(tukeydf$comparison, split="-")[,1]
  
  # use Tukey later, but test with the example given
  # tmpdf <- compare_means(as.formula(paste0(gene, " ~ category_factor")), 
  #                        data=ras_tmp)
  my_comparisons <- list( c("Normal", "Tumor_Stage_I"), 
                          c("Normal", "Tumor_Stage_IV"), 
                          c("Tumor_Stage_I", "Tumor_Stage_IV") )
  # set y positions on top of highest notch
  p1 <- ggplot(data=ras_tmp)+
    geom_boxplot(mapping=aes(x=category_factor, y=.data[[gene]]))+
    scale_x_discrete(labels=c("Normal","Stage I","Stage IV"))+
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=1),
          plot.margin = unit(c(0.05, 0.05, 0.3, 0.05), 
                                "inches"))
  ymax <- median(ras_tmp[[gene]]) + 2 * IQR(ras_tmp[[gene]])
  ymin <- median(ras_tmp[[gene]]) - 2 * IQR(ras_tmp[[gene]])
  ystep <- (ymax-ymin)/10
  p2 <- p1+ stat_pvalue_manual(tukeydf, label = "p.signif", 
    y.position = c(ymax+0.1*ystep, ymax+1.1*ystep, ymax+2.1*ystep))
  plotlist[[gene]] <- p2
}
# create a plot with normal, tumor and stage I
plot_all <- cowplot::plot_grid(plotlist=plotlist,
                                nrow=ceiling(sqrt(length(plotlist))))
figfile <- file.path(figfolder, "tukey_grid.pdf")
ggsave(figfile, plot=plot_all, width=20, height=20, limitsize = F)



## ----lower sample size---------------------------------------------------------------------------------------
gene_list <- c("PLCE1", "RGL1", "RAF1")
ras_tmp <- ras_data[ras_data$category2 %in% 
                      c("Normal", "Tumor_Stage_I", "Tumor_Stage_IV"),
                    c(gene_list, "category_factor")]
# now select the given amount of samples
# option  1) 8 normal, 5 stage I, 10 stage IV
# option 2) 5 normal, 10 stage I, 6 stage IV
# randomly select these number of samples
normal_index <- which(ras_tmp$category_factor=="Normal")
normal_selected <- sample(normal_index, size=5)
stageI_index <- which(ras_tmp$category_factor=="Tumor_Stage_I")
stageI_selected <- sample(stageI_index, size=10)
stageIV_index <- which(ras_tmp$category_factor=="Tumor_Stage_IV")
stageIV_selected <- sample(stageIV_index, size=6)

ras_tmp2 <- ras_tmp[c(normal_selected, stageI_selected, stageIV_selected),]
# scale the new dataframe
ras_tmp2[, gene_list] <- scale(ras_tmp2[, gene_list])
# repeat the anova plots
plotlist <- list()
for (gene in gene_list){
  stagelm <- lm(as.formula(paste0(gene, " ~ category_factor")), data=ras_tmp2)
  stageav <- aov(stagelm)
  summary(stageav)
  tukeydf <- as_tibble(TukeyHSD(stageav, conf.level=.95)$category_factor, 
                       rownames="comparison")
  tukeydf$p.signif <- stars.pval(tukeydf$`p adj`)
  tukeydf$p.signif[tukeydf$p.signif==" "] <- "ns"
  tukeydf$group1 <- limma::strsplit2(tukeydf$comparison, split="-")[,2]
  tukeydf$group2 <- limma::strsplit2(tukeydf$comparison, split="-")[,1]
  
  # use Tukey later, but test with the example given
  tmpdf <- compare_means(as.formula(paste0(gene, " ~ category_factor")), data=ras_tmp)
  my_comparisons <- list( c("Normal", "Tumor_Stage_I"), 
                          c("Normal", "Tumor_Stage_IV"), 
                          c("Tumor_Stage_I", "Tumor_Stage_IV") )
  # set y positions on top of highest notch
  p1 <- ggplot(data=ras_tmp2)+
    geom_boxplot(mapping=aes(x=category_factor, y=.data[[gene]]))+
    scale_x_discrete(labels=c("Normal","Stage I","Stage IV"))+
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=1),
          plot.margin = unit(c(0.05, 0.05, 0.3, 0.05), 
                                "inches"))
  ymax <- median(ras_tmp2[[gene]]) + 2 * IQR(ras_tmp2[[gene]])
  ymin <- median(ras_tmp2[[gene]]) - 2 * IQR(ras_tmp2[[gene]])
  ystep <- (ymax-ymin)/10
  p2 <- p1+ stat_pvalue_manual(tukeydf, label = "p.signif", 
    y.position = c(ymax+0.1*ystep, ymax+1.1*ystep, ymax+2.1*ystep))
  plotlist[[gene]] <- p2
}
# create a plot with normal, tumor and stage I
plot_all <- cowplot::plot_grid(plotlist=plotlist, nrow=1)
figfile <- file.path(figfolder, "best_genes_reduced_tukey_normal_vs_stageI.pdf")
ggsave(figfile, plot=plot_all, width=14, height=5)

# can we seperate normal from tumor stage I using these genes?
ras_tmp3 <- ras_tmp2[ras_tmp2$category_factor %in% c("Normal", "Tumor_Stage_I"),]
ras_tmp3$category_factor <- droplevels(ras_tmp3$category_factor)
it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                    paste(gene_list,collapse="+"))),
                  data = ras_data, family = "binomial")
conf_mx <- runCV(ras_tmp3, it_model,
                 depvar="category_factor",
                 depvar_levels=c("Normal", "Tumor_Stage_I"),
                 posclass="Tumor_Stage_I",
                 return_conf_mx=T)
#print(conf_mx)
# plot conf_mx
TClass <- factor(c(0, 0, 1, 1))
PClass <- factor(c(0, 1, 0, 1))
# y <- c(FN, TP, TN, FP)
y <- rep(0, times=4)
y[1] <- conf_mx$table[2]
y[2] <- conf_mx$table[1]
y[3] <- conf_mx$table[4]
y[4] <- conf_mx$table[3]
df <- data.frame(TClass, PClass, y)

p1 <- ggplot(data =  df, mapping = aes(x = TClass, y = PClass)) +
  geom_tile(aes(fill = y), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1, size=10) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none")

outfig <- file.path(figfolder, "reduced_data_best_conf_mx_normal_vs_stageI.pdf")
ggsave(outfig, plot=p1, width=4, height=4)


## ----lower sample size stage I vs stage IV-------------------------------------------------------------------
#gene_list <- c("RASIP1", "RASSF1", "ARHGAP20")
gene_list <- names(log_best_list[[2]]$coefficients)[
                       2:length(names(log_best_list[[2]]$coefficients))]
ras_tmp <- ras_data[ras_data$category2 %in% 
                      c("Normal", "Tumor_Stage_I", "Tumor_Stage_IV"),
                    c(gene_list, "category_factor")]
# now select the given amount of samples
# option  1) 8 normal, 5 stage I, 10 stage IV
# option 2) 5 normal, 10 stage I, 6 stage IV
# randomly select these number of samples
normal_index <- which(ras_tmp$category_factor=="Normal")
normal_selected <- sample(normal_index, size=5)
stageI_index <- which(ras_tmp$category_factor=="Tumor_Stage_I")
stageI_selected <- sample(stageI_index, size=10)
stageIV_index <- which(ras_tmp$category_factor=="Tumor_Stage_IV")
stageIV_selected <- sample(stageIV_index, size=6)

ras_tmp2 <- ras_tmp[c(normal_selected, stageI_selected, stageIV_selected),]
# scale the new dataframe
ras_tmp2[, gene_list] <- scale(ras_tmp2[, gene_list])
# repeat the anova plots
plotlist <- list()
for (gene in gene_list){
  stagelm <- lm(as.formula(paste0(gene, " ~ category_factor")), data=ras_tmp2)
  stageav <- aov(stagelm)
  summary(stageav)
  tukeydf <- as_tibble(TukeyHSD(stageav, conf.level=.95)$category_factor, 
                       rownames="comparison")
  tukeydf$p.signif <- stars.pval(tukeydf$`p adj`)
  tukeydf$p.signif[tukeydf$p.signif==" "] <- "ns"
  tukeydf$group1 <- limma::strsplit2(tukeydf$comparison, split="-")[,2]
  tukeydf$group2 <- limma::strsplit2(tukeydf$comparison, split="-")[,1]
  
  # use Tukey later, but test with the example given
  tmpdf <- compare_means(as.formula(paste0(gene, " ~ category_factor")), data=ras_tmp)
  my_comparisons <- list( c("Normal", "Tumor_Stage_I"), 
                          c("Normal", "Tumor_Stage_IV"), 
                          c("Tumor_Stage_I", "Tumor_Stage_IV") )
  # set y positions on top of highest notch
  p1 <- ggplot(data=ras_tmp2)+
    geom_boxplot(mapping=aes(x=category_factor, y=.data[[gene]]))+
    scale_x_discrete(labels=c("Normal","Stage I","Stage IV"))+
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=1),
          plot.margin = unit(c(0.05, 0.05, 0.3, 0.05), 
                                "inches"))
  ymax <- median(ras_tmp2[[gene]]) + 2 * IQR(ras_tmp2[[gene]])
  ymin <- median(ras_tmp2[[gene]]) - 2 * IQR(ras_tmp2[[gene]])
  ystep <- (ymax-ymin)/10
  p2 <- p1+ stat_pvalue_manual(tukeydf, label = "p.signif", 
    y.position = c(ymax+0.1*ystep, ymax+1.1*ystep, ymax+2.1*ystep))
  plotlist[[gene]] <- p2
}
# create a plot with normal, tumor and stage I
plot_all <- cowplot::plot_grid(plotlist=plotlist, nrow=1)
figfile <- file.path(figfolder, "best_genes_reduced_tukey_stageI_vs_stageIV.pdf")
ggsave(figfile, plot=plot_all, width=14, height=5)

# can we seperate tumor stage I from stage IV using these genes?
ras_tmp3 <- ras_tmp2[ras_tmp2$category_factor %in% c("Tumor_Stage_I", "Tumor_Stage_IV"),]
ras_tmp3$category_factor <- droplevels(ras_tmp3$category_factor)
it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                    paste(gene_list,collapse="+"))),
                  data = ras_data, family = "binomial")
conf_mx <- runCV(ras_tmp3, it_model,
                 depvar="category_factor",
                 depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                 posclass="Tumor_Stage_IV",
                 return_conf_mx=T)
#print(conf_mx)
# plot conf_mx
TClass <- factor(c(0, 0, 1, 1))
PClass <- factor(c(0, 1, 0, 1))
# y <- c(FN, TP, TN, FP)
y <- rep(0, times=4)
y[1] <- conf_mx$table[2]
y[2] <- conf_mx$table[1]
y[3] <- conf_mx$table[4]
y[4] <- conf_mx$table[3]
df <- data.frame(TClass, PClass, y)

p2 <- ggplot(data =  df, mapping = aes(x = TClass, y = PClass)) +
  geom_tile(aes(fill = y), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1, size=10) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none")

outfig <- file.path(figfolder, "reduced_data_best_conf_mx_stageI_vs_stageIV.pdf")
ggsave(outfig, plot=p2, width=4, height=4)


## ----lower sample size repeated with more genes normal vs stageI---------------------------------------------
#gene_list <- c("PLCE1", "RGL1", "RAF1")
#gene_list <- c("PLCE1", "RGL2", "RAF1")
gene_list <- c("PLCE1", "RAF1")
test_char <- paste0(gene_list, collapse="_")
ras_tmp <- ras_data[ras_data$category2 %in% 
                      c("Normal", "Tumor_Stage_I", "Tumor_Stage_IV"),
                    c(gene_list, "category_factor")]
# now select the given amount of samples
# option  1) 8 normal, 5 stage I, 10 stage IV
# option 2) 5 normal, 10 stage I, 6 stage IV
# randomly select these number of samples
num_boot <- 1000
acc_list <- list()
for (ii in 1:num_boot){
  if (ii %% 10==0) print(ii)
  normal_index <- which(ras_tmp$category_factor=="Normal")
  normal_selected <- sample(normal_index, size=5, replace=T)
  stageI_index <- which(ras_tmp$category_factor=="Tumor_Stage_I")
  stageI_selected <- sample(stageI_index, size=10, replace=T)
  stageIV_index <- which(ras_tmp$category_factor=="Tumor_Stage_IV")
  stageIV_selected <- sample(stageIV_index, size=6, replace=T)
  
  ras_tmp2 <- ras_tmp[c(normal_selected, stageI_selected, stageIV_selected),]
  # scale the new dataframe
  ras_tmp2[, gene_list] <- scale(ras_tmp2[, gene_list])
  # can we seperate normal from tumor stage I using these genes?
  ras_tmp3 <- ras_tmp2[ras_tmp2$category_factor %in% c("Normal", "Tumor_Stage_I"),]
  ras_tmp3$category_factor <- droplevels(ras_tmp3$category_factor)
  it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                      paste(gene_list,collapse="+"))),
                    data = ras_tmp3, family = "binomial")
  acc <- runCV_LOOCV(ras_tmp3, it_model,
                   depvar="category_factor",
                   depvar_levels=c("Normal", "Tumor_Stage_I"),
                   posclass="Tumor_Stage_I",
                   return_conf_mx=F)
  acc_list[[ii]] <- acc
}
ba_vector_normal_vs_stageI <- unlist(acc_list)
ptitle <- paste0(it_model$formula)[3]
p <- ggplot(data=as.data.frame(ba_vector_normal_vs_stageI), 
            mapping=aes(x=ba_vector_normal_vs_stageI))+
  geom_histogram(binwidth=0.05)+
  xlab("balanced accuracy")+
  ggtitle(ptitle)
suppressWarnings(plot(p))
figfile <- file.path(figfolder, 
        paste("histogram_bootsrap_normal_vs_stageI_", test_char,".pdf"))
suppressWarnings(ggsave(figfile, p, width=3, height=3))


## ----lower sample size repeated with more genes stageI vs stageIV--------------------------------------------
# make sure the gene list is obtained from running BIC for stage I vs stage IV
# still TO DO properly!!!
gene_list <- names(
log_best_list_stageI_vs_stageIV <- log_best_list[[2]]$coefficients)[
                       2:length(names(log_best_list[[2]]$coefficients))]
#gene_list <- c("RGL2", "RASSF1", "RADIL", "RASSF7", "RAPGEF4", "KRAS", "RIN3")
# gene_list <- c("RGL2", "RASSF1", "RASSF7", "KRAS")
# gene_list <- c("NRAS", "PIK3CA", "RALGDS", "RGL2", "RASSF1", "RASSF7", 
#                   "MLLT4")
test_char <- paste0(gene_list, collapse="_")
ras_tmp <- ras_data[ras_data$category2 %in% 
                      c("Normal", "Tumor_Stage_I", "Tumor_Stage_IV"),
                    c(gene_list, "category_factor")]
# now select the given amount of samples
# option  1) 8 normal, 5 stage I, 10 stage IV (104)
# option 3) 5 normal, 10 stage I, 6 stage IV (102)
# randomly select these number of samples
num_boot <- 10000
acc_list <- list()
for (ii in 1:num_boot){
  if (ii %% 100==0) print(ii)
  normal_index <- which(ras_tmp$category_factor=="Normal")
  normal_selected <- sample(normal_index, size=5, replace=T)
  stageI_index <- which(ras_tmp$category_factor=="Tumor_Stage_I")
  stageI_selected <- sample(stageI_index, size=10, replace=T)
  stageIV_index <- which(ras_tmp$category_factor=="Tumor_Stage_IV")
  stageIV_selected <- sample(stageIV_index, size=6, replace=T)
  
  ras_tmp2 <- ras_tmp[c(normal_selected, stageI_selected, stageIV_selected),]
  # scale the new dataframe
  ras_tmp2[, gene_list] <- scale(ras_tmp2[, gene_list])
  # can we seperate tumor stage I from stage IV using these genes?
  ras_tmp3 <- ras_tmp2[ras_tmp2$category_factor %in% c("Tumor_Stage_I", "Tumor_Stage_IV"),]
  ras_tmp3$category_factor <- droplevels(ras_tmp3$category_factor)
  it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                      paste(gene_list,collapse="+"))),
                    data = ras_tmp3, family = "binomial")
  acc <- runCV_LOOCV(ras_tmp3, it_model,
                   depvar="category_factor",
                   depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
                   posclass="Tumor_Stage_IV",
                   return_conf_mx=F)
  acc_list[[ii]] <- acc
}
ba_vector_stageI_vs_stageIV <- unlist(acc_list)
ptitle <- paste0(it_model$formula)[3]
pmedian <- median(ba_vector_stageI_vs_stageIV)
p <- ggplot(data=as.data.frame(ba_vector_stageI_vs_stageIV), 
            mapping=aes(x=ba_vector_stageI_vs_stageIV))+
  geom_histogram(binwidth=0.05)+
  geom_vline(xintercept=pmedian, color="red")+
  xlab("balanced accuracy")+
  ggtitle(ptitle)
suppressWarnings(plot(p))
figfile <- file.path(figfolder, 
        paste("histogram_bootsrap_stageI_vs_stageIV_", test_char,".pdf"))
suppressWarnings(ggsave(figfile, p, width=6, height=3))

# plot results together?


## ----reduce existing tables----------------------------------------------------------------------------------
filename1 <- "best_model_2genes_normal_vs_stageI.xlsx"
filename2 <- "best_model_2genes_stageI_vs_stageIV.xlsx"
for (filename in c(filename1, filename2)){
  newfile <- file.path(phome, "Results" ,filename)
  if (!file.exists(newfile)){
    oldfile <-  file.path(phome, "Results", "archive_new", filename)
    oldtable <- read_excel(oldfile, sheet=rtype)
    omitrows <- c()
    for (ri in 1:nrow(oldtable)){
      if (oldtable$gene1[ri] %in% exclude | oldtable$gene2[ri] %in% exclude){
        omitrows <- c(omitrows, ri)
      }
    }
    newtable <- oldtable[-omitrows,]
    outlist <- list()
    outlist[[rtype]] <- newtable
    openxlsx::write.xlsx(outlist, file=newfile)
  }
}


## ----best optimized model------------------------------------------------------------------------------------
ras_tmp <- ras_data[ras_data$category2 %in% 
                      c("Tumor_Stage_I", "Tumor_Stage_IV"),
                    c(indep_vars, "category_factor")]
# also delete unused factors
ras_tmp$category_factor <- droplevels(ras_tmp$category_factor)
# these genes were obtained by the optimization
# test_genes <- c("RGL2", "RASSF1", "RADIL", "RASSF7", "RAPGEF4", "KRAS", "RIN3",
#                 "PLCE1", "RAF1")
test_genes <- c("RGL2", "RASSF1", "RADIL", "RASSF7", "RAPGEF4", "KRAS", "RIN3")
test_char <- paste0(test_genes, collapse="_")
it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                    paste(test_genes,collapse="+"))),
                  data = ras_tmp, family = "binomial")
red_conf_mx <- runCV_LOOCV(ras_tmp, it_model,
         depvar="category_factor",
         depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
         posclass="Tumor_Stage_IV",
         return_conf_mx=T)
print(red_conf_mx)
# also plot
  conf_mx <- red_conf_mx
  TClass <- factor(c(0, 0, 1, 1))
  PClass <- factor(c(0, 1, 0, 1))
  # y <- c(FN, TP, TN, FP)
  y <- rep(0, times=4)
  y[1] <- conf_mx$table[2]
  y[2] <- conf_mx$table[1]
  y[3] <- conf_mx$table[4]
  y[4] <- conf_mx$table[3]
  df <- data.frame(TClass, PClass, y)
  
  p <- ggplot(data =  df, mapping = aes(x = TClass, y = PClass)) +
    geom_tile(aes(fill = y), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", y)), vjust = 1, size=10) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() + theme(legend.position = "none")
  
  plot(p)
  print(conf_mx)
  outfig <- file.path(figfolder, 
                      paste("best_conf_mx_", test_char, ".pdf", sep=""))
  ggsave(outfig, plot=p, width=7, height=7)


# also seperate tumor stage I from stage IV
indata_mx <-  as.matrix(ras_tmp[, test_genes])
lda_model <- MASS::lda(indata_mx, grouping=ras_tmp$category_factor)
lda_values <- predict(lda_model, newdata = ras_tmp[, test_genes])
all.equal(lda_values$class, ras_tmp$category_factor) # "14 string mismatches"
lda_annot <- as_tibble(lda_values$x)
lda_annot$category_factor <- ras_tmp$category_factor

# obtain the “proportion of trace” which is the proportion of between-class variance that is explained by successive discriminant functions. 
prop_trace <- prop.table(lda_model$svd^2)
nicelabel <- list()
for (ii in 1:length(prop_trace)){
 nicelabel[[ii]] <- paste('LD',
        as.character(ii), ' (', 
        as.character(round(100*prop_trace[ii], digits = 2)), '%)', sep="")
}
p <- ggplot(data=lda_annot,  mapping=aes(x=LD1, fill=category_factor)) +
  geom_histogram()

# plot density curve for all stages
lda_values <- predict(lda_model, newdata = ras_data[, test_genes])
lda_annot <- as_tibble(lda_values$x)
lda_annot$sample <- ras_data$sample
lda_annot$category_factor <- ras_data$category_factor

lda_tmp <- lda_annot[!lda_annot$category_factor %in% 
                       c("Normal", "Tumor_Stage_X"),]
lda_tmp$category_factor <- droplevels(lda_tmp$category_factor)

p <- ggplot(data=lda_tmp, mapping=aes(x=LD1, color=category_factor))+
  geom_density(fill="white",  alpha=0.25)+
      scale_color_manual("category",
            values = c("Normal" = "black",
                       "Tumor_Stage_I"="cyan",
                       "Tumor_Stage_II"="purple",
                       "Tumor_Stage_III"="pink",
                       "Tumor_Stage_IV"="red"),
            labels= c("Normal" = "Normal",
                       "Tumor_Stage_I"="Tumor Stage I",
                       "Tumor_Stage_II"="Tumor Stage II",
                       "Tumor_Stage_III"="Tumor Stage III",
                       "Tumor_Stage_IV"="Tumor Stage IV"))
plot(p)
outfig <- file.path(figfolder, 
                    paste("ras_all_lda_density",test_char,".pdf", sep=""))
ggsave(outfig, plot=p, width=7, height=7)

