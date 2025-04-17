## ----setup, include=F----------------------------------------------------------------------------------------
knitr::opts_chunk$set(message=FALSE, warnings=FALSE, fig.width=12, fig.height=5)


## ----libraries-----------------------------------------------------------------------------------------------
#library(readxl) # another package to read xls files
library(dplyr)
library(tidyr) # to transform data for ggplot2 plots
library(data.table) # while I use tibble mostly, fread is best input option
library(tibble)
library(ggplot2) # for visualization
library(EnvStats) # for statistic tests (i.e. Rosner test)
library(ggcorrplot) # for nice correlation plots
library(cowplot) # to plot q-q plots in a grid
library(gtools) # specifically to generate stars from significance values
library(caret) # for Classification And REgression Training
library(limma) # for the strsplit2 function
library(gridExtra) # to plot a list of figures in a grid, also for multipage pdf
#library(ggpubr) # for the same purpose to arrange figures in a grid and also for significance stars
library(arsenal) # Functions for Large-Scale Statistical Summaries
#library(dr4pl) # for ROUT outlier finding
#library(gridGraphics) # Functions to convert a page of plots drawn with the 'graphics' package into identical output drawn with the 'grid' package
library(ggpubr) # Add manually p-values to a ggplot
library(rstatix) # Autocompute P-value Positions For Plotting Significance & convert p-values to stars
#library(factoextra) # to plot elbow plot for PCA eigenvalues


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


## ----read data for each cancer type---------------------------------------------------------------------------------------------
# rtype <- "COAD"
ras_data_list <- list()
ras_scaled_list <- list()
pie_chart_plots <- list()
for (rtype in c("COAD", "LUAD", "LUSC")){
  # create folders for each seperate tumor type
  if (!file.exists(file.path(figfolder, rtype))) dir.create(figfolder)
  if (!file.exists(file.path(calcfolder, rtype))) dir.create(calcfolder)
  if (!file.exists(file.path(resfolder, rtype))) dir.create(resfolder)
  
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
  data$samples <- as.character(data$Var1)
  for (ii in 1:nrow(data)){
    if (length(grep("Tumor", data$samples[ii]))>0){
      data$samples[ii] <- gsub("Tumor_", "", data$samples[ii])
      data$samples[ii] <- gsub("_", " ", data$samples[ii])
    }
  }
  data$allinfo <- paste(data$samples, data$percent, sep=" ")
  p <- ggplot(data, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(stat="identity", width=1, color="white") +
    geom_text(aes(x=1.2, label = Freq),
              position = position_stack(vjust = 0.5),
              color="white", size=5) +
    coord_polar("y", start=0) +
    scale_fill_manual("", 
                       values = c("Normal" = "black",
                                  "Tumor_Stage_I"=stage1col,
                                  "Tumor_Stage_II"=stage2col,
                                  "Tumor_Stage_III"=stage3col,
                                  "Tumor_Stage_IV"=stage4col),
                      # labels= c("Normal" = "Normal",
                      #           "Tumor_Stage_I"="Stage I",
                      #           "Tumor_Stage_II"="Stage II",
                      #           "Tumor_Stage_III"="Stage III",
                      #           "Tumor_Stage_IV"="Stage IV")
                      labels=data$allinfo
                      )+
    theme_void() # remove background, grid, numeric labels
  #plot(p)
  # outfig <- file.path(figfolder, rtype, "pie_chart.pdf")
  # ggsave(outfig, p, width=8, height = 4)
  pie_chart_plots[[rtype]] <- p
  
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

# plot pie-charts together
pie_charts <- cowplot::plot_grid(plotlist = pie_chart_plots, 
                  nrow=1, labels=names(pie_chart_plots))
outfig <- file.path(figfolder, "pie_charts.pdf")
ggsave(outfig, plot=pie_charts, width=11.7, height=3.3)
  

## ----PCA & Rosner test for outlier detection  ---------------------------------
outlier_plots <- list()
for (rtype in c("COAD", "LUAD", "LUSC")){
  outlier_plots[[rtype]] <- list()
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
  #plot(p)
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
  #plot(p)
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
  #plot(pfinal)
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
  #plot(p_labels)
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
  #plot(pfinal)
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
    # ggsave(figfile, plot=p,  width=7*colfigs, height=5*colfigs)
    outlier_plots[[rtype]][[ctype]] <- p
  } # loop for category factors
  # save figures together
  p_rtype <- cowplot::plot_grid(plotlist=outlier_plots[[rtype]], 
                  ncol=1, labels=names(outlier_plots[[rtype]]))
  figfile <- file.path(figfolder, rtype, 
                       paste(rtype, "_outliers_all.pdf", sep=""))
  ggsave(figfile, plot=p_rtype, height=5*8.3, width=11.7)
                       
  
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
  #plot(pfinal)
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
      outfig2 <- file.path(cfigfolder, 
              paste(rtype, "_correlation_",  stype, "_with_legend.tiff", sep=""))
      ggsave(outfig2, plot=p, height=10, width=9, dpi=600)
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
      #ylab("Number of \n connections")+ xlab("correlation threshold")+
      xlab("")+ ylab("")+
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
    
    # save each small figure seperately
    outfig  <- file.path(figfolder, paste("lost_connections", 
                                          as.character(pindex), "_", rtype, "_",
                                          comptype, "_minus_plus_seperately_top.tiff", sep=""))
    ggsave(outfig, plot=plist[[pindex]], dpi=600, width=3, height = 7)
    
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
      # xlab("correlation threshold")+
      xlab("")+ ylab("")+
      # theme(axis.text=element_text(size=16),
      #       axis.title=element_text(size=20),
      #   axis.title.x=element_blank(),
      #   axis.text.x=element_blank(),
      #   axis.ticks.x=element_blank(),
      #   axis.title.y=element_blank(),
      #   axis.text.y=element_blank(),
      #   axis.ticks.y=element_blank())+
      # theme(axis.title.y=element_blank(),
      #       axis.text.y=element_blank(),
      #       axis.ticks.y=element_blank()) %+replace%
      theme() %+replace%
      la_replace
      la_replace
    #labs(title = rtype)
    
    # save each small figure seperately
    outfig  <- file.path(figfolder, paste("lost_connections", 
                                          as.character(pindex), "_", rtype, "_",
                                          comptype, "_minus_plus_seperately_top.tiff", sep=""))
    ggsave(outfig, plot=plist[[pindex]], dpi=600, width=3, height = 7)
    
    pindex <- pindex+1
  }
  # p_top <- cowplot::plot_grid(plotlist=plist, nrow=1, 
  #                             rel_widths = rep(c(1.25,1), times=3))
  # plot(p_top)
  # # save
  # outfig <- file.path(figfolder, 
  #                     paste("lost_connections_",
  #                     comptype, "_minus_plus_seperately_top.pdf", sep=""))
  # ggsave(outfig, plot=p_top, width=21, height=7)
  
  plist_bottom <- list()
  pindex <- 1
  for (rtype in names(corr_list)){
    plist_bottom[[pindex]] <- ggplot(data=num_list[[rtype]])+
      geom_line(mapping=aes(x=corr, y=diff), color="darkgreen",
                linewidth=linewidth_vector[pindex], linetype=linetype_vector[pindex])+
      geom_point(mapping=aes(x=corr, y=diff), color="darkgreen", size=psize)+
      xlim(-0.8,-rlim)+ ylim(0,ylim2)+
      #ylab("Number of \n lost connections")+xlab("correlation threshold")+
      ylab("")+xlab("")+
      # theme(axis.text=element_text(size=16),
      #       axis.title=element_text(size=20))
      theme() %+replace%
      la_replace
    
    # save each small figure seperately
    outfig  <- file.path(figfolder, paste("lost_connections", 
                                          as.character(pindex), "_", rtype, "_",
                                          comptype, "_minus_plus_seperately_bottom.tiff", sep=""))
    ggsave(outfig, plot=plist_bottom[[pindex]], dpi=600, width=3, height = 7)
    
    pindex <- pindex + 1
    plist_bottom[[pindex]] <- ggplot(data=num_list[[rtype]])+
      geom_line(mapping=aes(x=corr, y=diff), color="darkgreen",
                linewidth=linewidth_vector[pindex], linetype=linetype_vector[pindex])+
      geom_point(mapping=aes(x=corr, y=diff), color="darkgreen", size=psize)+
      xlim(rlim,0.8)+ ylim(0,ylim2)+
      #xlab("correlation threshold")+
      xlab("")+ylab("")+
      # theme(axis.text=element_text(size=16),
      #       axis.title=element_text(size=20),
      #       axis.title.y=element_blank(),
      #       axis.text.y=element_blank(),
      #       axis.ticks.y=element_blank())
      # theme(axis.title.y=element_blank(),
      #       axis.text.y=element_blank(),
      #       axis.ticks.y=element_blank()) %+replace%
      theme() %+replace%
      la_replace
      la_replace
    
    # save each small figure seperately
    outfig  <- file.path(figfolder, paste("lost_connections", 
                                          as.character(pindex), "_", rtype, "_",
                                          comptype, "_minus_plus_seperately_bottom.tiff", sep=""))
    ggsave(outfig, plot=plist_bottom[[pindex]], dpi=600, width=3, height = 7)
    
    pindex <- pindex+1
  }
  # p_bottom <- cowplot::plot_grid(plotlist=plist_bottom, nrow=1, 
  #                                rel_widths = rep(c(1.25,1), times=3))
  # plot(p_bottom)
  # # save
  # outfig <- file.path(figfolder, paste("lost_connections_",
  #                             comptype, "_minus_plus_seperately_bottom.pdf", sep=""))
  # ggsave(outfig, plot=p_bottom, width=21, height=7)
  
  # plist_all <- do.call(c, list(plist, plist_bottom))
  # p_all <- cowplot::plot_grid(plotlist=plist_all, nrow=2, 
  #                             rel_widths = rep(c(1.25,1), times=6))
  # plot(p_all)
  # outfig <- file.path(figfolder, paste("lost_connections_",
  #                                      comptype, "_minus_plus_seperately.pdf", sep=""))
  # ggsave(outfig, plot=p_all, width=21, height=14)
  
  # write out to file
  outfile <- file.path(resfolder, paste("lost_connections_",
                                    comptype, "_minus_plus_seperately.xlsx", sep=""))
  #openxlsx::write.xlsx((Reduce(bind_rows, num_list)), file=outfile)
  openxlsx::write.xlsx(num_list, file=outfile)
  
  print(comptype)
  print(round(unlist(mean_connections), digits=2))
  print(unlist(mean_connections_diff))
  
}




# ##-----create violinplots of gene expressions of different stages--------##
# 
# plotlist <- list()
# for (rtype in c("COAD", "LUAD", "LUSC")){
#   ras_data <- ras_scaled_list[[rtype]]
#   ras_tmp <- ras_scaled_list[[rtype]][
#     ras_scaled_list[[rtype]]$category_factor %in% 
#       c("Normal", "Tumor_Stage_I", "Tumor_Stage_II", 
#         "Tumor_Stage_III", "Tumor_Stage_IV"),]
#   ras_tmp$category_factor <- droplevels(ras_tmp$category_factor)
#   ras_tmp$category_factor <- factor(ras_tmp$category_factor, levels=
#                                       c("Normal", "Tumor_Stage_I", "Tumor_Stage_II", 
#                                         "Tumor_Stage_III", "Tumor_Stage_IV"))
#   plotlist[[rtype]] <- list()
#   for (gene in indep_vars){
#     stagelm <- lm(as.formula(paste0(gene, " ~ category_factor")), data=ras_tmp)
#     stageav <- aov(stagelm)
#     summary(stageav)
#     tukeydf <- as_tibble(TukeyHSD(stageav, conf.level=.95)$category_factor, 
#                          rownames="comparison")
#     tukeydf$p.signif <- stars.pval(tukeydf$`p adj`)
#     tukeydf$p.signif[tukeydf$p.signif==" "] <- "ns"
#     tukeydf$group1 <- limma::strsplit2(tukeydf$comparison, split="-")[,2]
#     tukeydf$group2 <- limma::strsplit2(tukeydf$comparison, split="-")[,1]
#     
#     p1 <- ggplot(ras_tmp,
#                  aes(x=category_factor, y=.data[[gene]], 
#                      colour = category_factor)) +
#       geom_violin(trim=F) +
#       geom_jitter(size=0.5) +
#       geom_boxplot(width=0.2)+
#       scale_color_manual("category",
#                          values = c("Normal" = "black",
#                                     "Tumor_Stage_I"="cyan",
#                                     "Tumor_Stage_II"="purple",
#                                     "Tumor_Stage_III"="pink",
#                                     "Tumor_Stage_IV"="red"),
#                          labels= c("Normal" = "Normal",
#                                    "Tumor_Stage_I"="Stage I",
#                                    "Tumor_Stage_II"="Stage II",
#                                    "Tumor_Stage_III"="Stage III",
#                                    "Tumor_Stage_IV"="Stage IV"))+
#       xlab("") + ylab(gene)+ guides(colour="none") + 
#       ylim(-6, 6)+
#       #ylim(-11, 6)+
#       # scale_x_discrete(labels=c("Normal", "Stage I", "Stage II", 
#       #                           "Stage III", "Stage IV"))+
#       scale_x_discrete(labels=c("", "", "", "", ""))+
#       theme(axis.title.x=element_blank(),
#             axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=1),
#             plot.margin = unit(c(0.05, 0.05, 0.3, 0.05), "inches"))
#     y.position <- 4.5
#     # testtbl <- as_tibble(list("group1"="Normal", "group2"="Tumor",
#     #                           "p"=stars.pval(tres$p.value), "y.position"=y.position))
#     testtbl <- tukeydf[1,]
#     colnames(testtbl)[colnames(testtbl)=="p.signif"] <- "p"
#     testtbl$y.position <- y.position
#     
#     p2 <- p1+ stat_pvalue_manual(testtbl)
#     plotlist[[rtype]][[gene]] <- p2
#     # plot to see warnings
#     tt <- tryCatch(plot(p2),error=function(e) e, warning=function(w) w)
#     if(is(tt,"warning")){
#       print(paste("warning in plot: ", rtype, ", ", gene, sep=""))
#     }
#   }
#   p_all <- cowplot::plot_grid(plotlist=plotlist[[rtype]], ncol=11)
#   outfig <- file.path(figfolder, rtype, 
#                       paste("tcga_", rtype, "violinplot_stages_with_stars.pdf", sep=""))
#   ggsave(outfig, plot=p_all, width=28, height=14)
# }


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
  outfig <- file.path(figfolder, rtype, "ras_LD1_density_categories.tiff")
  ggsave(outfig, plot=p, width=7, height=5, dpi=600)
  
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
  outfig <- file.path(figfolder, rtype, "ras_LD2_density_categories.tiff")
  ggsave(outfig, plot=p, width=7, height=5, dpi=600)
  
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
  outfig <- file.path(figfolder, rtype, "ras_LD2_density_categories_noNormal.tiff")
  ggsave(outfig, plot=p, width=7, height=5, dpi=600)
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
        # geom_col(aes(fill = barcolor), color="black", width = 0.8) +
        # scale_fill_manual("category",
        #                   values = c("#d4d4d4" = "#d4d4d4",
        #                              "#CD3232"="#CD3232",
        #                              "#33ce33"="#33ce33"))+
        geom_col(fill = "white", color="black", width = 0.8) +
        ## add percentage labels
        # geom_text(aes(label = genes, x=0.2), nudge_x = -0.05) +
        geom_text(aes(label = genes, x=0.02, color=barcolor), hjust=0) +
        scale_color_manual("category",
                          values = c("#d4d4d4" = "black",
                                     "#CD3232"="#CD3232",
                                     "#33ce33"="#33ce33"))+
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
      plot(plotlist[[k]])
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
  figfile2 <- file.path(figfolder, rtype, paste("k", "_gene_model_measures_",
                                               comptype, ".tiff", sep=""))
  ggsave(figfile2, plot=p, width=12, height=9, dpi=600)
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
figfile2 <- file.path(figfolder, rtype, 
                            "k_gene_model_measures_stageI_vs_stageIV.tiff")
ggsave(figfile2, plot=p_all, width=21, height=14, dpi=600)

## ----create violinplots for all genes considering all stages ---------------------------------
# read in order file
gorder_file <- file.path(datafolder, "gene_list.csv")
gorder <- read.csv(gorder_file, header = F)[, 1, drop=T]
# # check
# print(length(gorder)-length(union(gorder, indep_vars)))
# print(length(gorder)-length(intersect(gorder, indep_vars)))

# replace with proper use of facet_grid

#rtype <- "COAD"
for (rtype in names(ras_data_list)){
  ras_tmp <- ras_data_list[[rtype]]
  ras_tmp <- pivot_longer(ras_tmp[, c("sample", "category_factor", indep_vars)], 
                          all_of(indep_vars),
                          names_to="genes", values_to="scaled_expr")
  
  # order of comparisons is crucial!
  ras_tmp$genes <- factor(ras_tmp$genes, levels = gorder)
  
  
  p <- ggplot(data=ras_tmp, aes(x=category_factor, y=scaled_expr, color=category_factor))+
    # geom_violin(bw=0.3, trim=F)+
    # geom_point(size=0.5, aes(x=category_factor, y=expression, color=category_factor),
    #            position=position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2))+
    geom_violin(trim=F) +
    geom_jitter(size=0.5, width = 0.2) +
    #geom_boxplot(width=0.2)+
    ylab("scaled gene expression") + xlab("")+
    guides(color=guide_legend(title=""))+
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
    scale_x_discrete(labels= c("Normal" = "Normal",
                               "Tumor_Stage_I"="Stage I",
                               "Tumor_Stage_II"="Stage II",
                               "Tumor_Stage_III"="Stage III",
                               "Tumor_Stage_IV"="Stage IV")) +
    theme(
      #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
      axis.text.x = element_text(angle = 90),
      legend.position="none",
      axis.title=element_text(size=20),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      axis.text=element_text(size=16),
      legend.text=element_text(size=16),
      strip.text.x = element_text(size=12), # facet titles
      panel.background=element_rect(fill="white", colour = "black"),
      plot.background=element_rect(color="black"))+
    #facet_grid(gene~experiment)
    facet_wrap(~genes, ncol=5, scale="free_y")
  #plot(p)
  # outfig <- file.path(figfolder, rtype, "violinplot_with_ttest.pdf")
  # ggsave(outfig, plot=p, width=21, height=14)
  
  
  # also calculate p-values - 
  # format needed for stat_pvalue_manual: group1 | group2 | p | y.position | etc.
  
  stat.test <- ras_tmp %>%
    group_by(genes) %>%
    t_test(scaled_expr ~ category_factor) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  #stat.test
  # keep Normal vs StageI only
  stat.test <- stat.test[which(stat.test$group1 %in% c("Normal") &
                                 stat.test$group2 %in% c("Tumor_Stage_I")),]
  
  # Add p-values onto the box plots
  stat.test <- stat.test %>%
    add_xy_position(x = "category_factor",
                    dodge = 0.8, step.increase = 0.06,
    )
  
  p_final <- p + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.005,
    bracket.nudge.y=-1
  )
  
  #plot(p_final)
  outfig <- file.path(figfolder, rtype, 
                paste("violinplot_", rtype, "_with_ttest.pdf", sep=""))
  ggsave(outfig, plot=p_final, width=14, height=14)
}





## ----create violinplots for all genes for normal vs tumor ----------------------------------

# replace with proper use of facet_grid

#rtype <- "COAD"
violinplot_list <- list()
violinplot_list_BH <- list()
for (rtype in names(ras_data_list)){
  ras_tmp <- ras_data_list[[rtype]]
  ras_tmp <- pivot_longer(ras_tmp[, c("sample", "sample_type_factor", indep_vars)], 
                          all_of(indep_vars),
                          names_to="genes", values_to="scaled_expr")
  
  # order of comparisons is crucial (change to predefined order - TO DO)
  ras_tmp$genes <- factor(ras_tmp$genes, levels = gorder)
  
  
  p <- ggplot(data=ras_tmp, aes(x=sample_type_factor, y=scaled_expr, 
                                color=sample_type_factor))+
    # geom_violin(bw=0.3, trim=F)+
    # geom_point(size=0.5, aes(x=category_factor, y=expression, color=category_factor),
    #            position=position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2))+
    geom_violin(trim=F) +
    geom_jitter(size=0.5, width = 0.2) +
    #geom_boxplot(width=0.2)+
    ylab("scaled gene expression") + xlab("")+
    guides(color=guide_legend(title=""))+
    scale_color_manual("category",
                       values = c("Normal" = "black",
                                  "Tumor"="red")) +
    theme(
      #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
      axis.text.x = element_text(angle = 90),
      #legend.position="none",
      axis.title=element_text(size=20),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      axis.text=element_text(size=16),
      legend.text=element_text(size=16),
      strip.text.x = element_text(size=12), # facet titles
      panel.background=element_rect(fill="white", colour = "black"),
      plot.background=element_rect(color="black"))+
    #facet_grid(gene~experiment)
    facet_wrap(~genes, ncol=5, scale="free_y")
  #plot(p)
  #grid.draw(shift_legend(p))
  
  # outfig <- file.path(figfolder, rtype, "violinplot_with_ttest.pdf")
  # ggsave(outfig, plot=p, width=21, height=14)
  
  
  # also calculate p-values - 
  # format needed for stat_pvalue_manual: group1 | group2 | p | y.position | etc.
  
  stat.test <- ras_tmp %>%
    group_by(genes) %>%
    t_test(scaled_expr ~ sample_type_factor) %>%
    adjust_pvalue(method = "none") %>% # try BH as well
    add_significance("p.adj")
  #stat.test
  
  # Add p-values onto the box plots
  stat.test <- stat.test %>%
    add_xy_position(x = "sample_type_factor",
                    dodge = 0.8, step.increase = 0.06,
    )
  stat.test$y.position <- 1.05*stat.test$y.position
  
  p_final <- p + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.005,
    bracket.nudge.y=-1
  )
  
  #plot(p_final)
  # outfig <- file.path(figfolder, rtype, 
  #                     paste("violinplot_", rtype, "_normal_vs_tumor_with_ttest.pdf", sep=""))
  # ggsave(outfig, plot=p_final, width=14, height=14)
  violinplot_list[[rtype]] <- p_final
  
  # add BH correction
  stat.test <- ras_tmp %>%
    group_by(genes) %>%
    t_test(scaled_expr ~ sample_type_factor) %>%
    adjust_pvalue(method = "BH") %>% # try BH as well
    add_significance("p.adj")
  #stat.test
  
  # Add p-values onto the box plots
  stat.test <- stat.test %>%
    add_xy_position(x = "sample_type_factor",
                    dodge = 0.8, step.increase = 0.06,
    )
  
  p_final <- p + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.005,
    bracket.nudge.y=-1
  )
  
  #plot(p_final)
  # outfig <- file.path(figfolder, rtype, 
  #                     paste("violinplot_", rtype, "_normal_vs_tumor_with_ttest_BH.pdf", sep=""))
  # ggsave(outfig, plot=p_final, width=14, height=14)
  violinplot_list_BH[[rtype]] <- p_final
}

# save in a multipage pdf
plot_none <- plot_grid(plotlist=violinplot_list, ncol=1, 
                     labels=names(violinplot_list))
outfig <- file.path(figfolder, "violinplot_normal_vs_tumor_with_ttest.pdf")
ggsave(outfig, plot=plot_none, height = 2*3*8.3, width=2*11.7)


plot_BH <- plot_grid(plotlist=violinplot_list_BH, ncol=1, 
                     labels=names(violinplot_list_BH))
outfig <- file.path(figfolder, "violinplot_normal_vs_tumor_with_ttest_BH.pdf")
ggsave(outfig, plot=plot_BH, height = 2*3*8.3, width=2*11.7)

