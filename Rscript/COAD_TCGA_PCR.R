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
library(arsenal) # Functions for Large-Scale Statistical Summaries
library(ggpubr) # Add manually p-values to a ggplot
library(rstatix) # Autocompute P-value Positions For Plotting Significance


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

# define color
stage1col <- "#14f2e0"
stage2col <- "#b5d6d6"
stage3col <- "#ceb5b7"
stage4col <- "#F320FA"


## ----cancer type---------------------------------------------------------------------------------------------
rtype <- "COAD"
# use a loop here

figfolder <- file.path(figfolder, rtype)
if (!file.exists(figfolder)) dir.create(figfolder)
calcfolder <- file.path(calcfolder, rtype)
if (!file.exists(calcfolder)) dir.create(calcfolder)
resfolder <- file.path(resfolder, rtype)
if (!file.exists(resfolder)) dir.create(resfolder)

##----read in TCGA data already saved in rds file-------
ras_list <- readRDS(
        file=file.path(calcfolder, paste("TCGA_", rtype, ".rds", sep="")))
ras_data <- ras_list[["data"]]
indep_vars <- ras_list[["indep_vars"]]

##----read PCR data to use negative delta CT---------------------------
pcr_genes <- c("GRB7", "PLCE1", "RAF1", "RGL1", "RIN1")
pcr_list <- list()
for (pgene in pcr_genes){
  fname <- file.path(datafolder, paste(pgene, "deltaCT.csv", sep="_"))
  indata <- as_tibble(read.csv2(fname))
  indata$X <- NULL
  indata[[pgene]] <- -indata$delta.CT
  indata$delta.CT <- NULL
  pcr_list[[pgene]] <- indata
}

for (pgene in names(pcr_list)){
  print(head(pcr_list[[pgene]]))
  print("##########")
}
ras_pcr <- Reduce(inner_join, pcr_list)
colnames(ras_pcr)[1:3] <- c("sample", "sample_type_factor", "category_factor")
ras_pcr$category_factor[ras_pcr$category_factor %in% c("0")] <- "Normal"
ras_pcr$category_factor[ras_pcr$category_factor %in% c("I")] <- "Tumor_Stage_I"
ras_pcr$category_factor[ras_pcr$category_factor %in% 
                        c("II", "IIA", "IIB")] <- "Tumor_Stage_II"
ras_pcr$category_factor[ras_pcr$category_factor %in% 
                        c("III", "IIIA", "IIIB", "IIIC")] <- "Tumor_Stage_III"
ras_pcr$category_factor[ras_pcr$category_factor %in% c("IV")] <- "Tumor_Stage_IV"
ras_pcr$category_factor <- factor(ras_pcr$category_factor)
ras_pcr$sample_type_factor <- factor(ras_pcr$sample_type_factor)

##--create differenty formats of the data----------
# create long data for ggplot visualizations
ras_tcga <- ras_data[, c("sample", "sample_type_factor", "category_factor", 
                         pcr_genes)]
# from now repeat steps for both types of data
ras_both <- list("tcga"=ras_tcga, "pcr"=ras_pcr)
ras_both_scaled <- list()
ras_both_long <- list()
for (rname in names(ras_both)){
  ras_both_scaled[[rname]] <- ras_both[[rname]]
  ras_both_scaled[[rname]][, pcr_genes] <- scale(
    ras_both[[rname]][, pcr_genes], center=T, scale=T)
  ras_both_long[[rname]] <- ras_both_scaled[[rname]] %>% 
    pivot_longer(all_of(pcr_genes), 
                 names_to = "gene", values_to = "expression")
  ras_both_long[[rname]]$experiment <- rname
}
ras_both_long_tbl <- bind_rows(ras_both_long[[1]], ras_both_long[[2]])

# to keep variables order:
ras_both_long_tbl$gene <- factor(ras_both_long_tbl$gene, levels=pcr_genes)


##----find outliers with Rosner test---------------
outlier_list <- list()
plist <- list()
p_all_list <- list()
ras_tmp <- ras_both_scaled[["pcr"]]
for (ctype in levels(ras_tmp$category_factor)){
  plist[[ctype]] <- list()
  outlier_list[[ctype]] <- c()
  for (ci in 1:length(pcr_genes)){
    # assume 1 outlier, test if it's true
    pcr_gene <- pcr_genes[ci]
    selected_samples <- ras_tmp[ras_tmp$category_factor==ctype, 
                c(1, match(pcr_gene, colnames(ras_tmp)))]
    rosner_result <- rosnerTest(selected_samples[,2,drop=T], k=1)
    otl <- rosner_result$n.outliers
    # draw quantile plot, even if there are no outliers
    #if (otl > 0){
    plotvector <- selected_samples[[pcr_gene]]
    names(plotvector) <- selected_samples$sample
    outlier_ids <- car::qqPlot(plotvector,
                               ylab=pcr_gene, 
                               id=list(method="y", n=otl, 
                                       cex=1, col=car::carPalette()[1], location="lr"))
    plist[[ctype]][[pcr_gene]] <- recordPlot()
    outlier_list[[ctype]] <- c(outlier_list[[ctype]], outlier_ids)
    # print(varname)
    # print(outlier_ids)
    # print("=========")
    #}
    
  }
  figfile <- file.path(figfolder, 
                       paste("outliers_pcr_", ctype, ".pdf", sep=""))
  colfigs <- ceiling(sqrt(length(plist[[ctype]])))
  p <- plot_grid(plotlist = plist[[ctype]], 
                 ncol=colfigs,
                 align = 'vh', scale=0.8)
  # ggsave(figfile, plot=p,  width=7*colfigs, height=5*colfigs)
  p_all_list[[ctype]] <- p
}
# plot together
p_all <- cowplot::plot_grid(plotlist=p_all_list, 
                              ncol=1, labels=names(p_all_list))
figfile <- file.path(figfolder, 
                     paste(rtype, "_pcr_outliers_all.pdf", sep=""))
ggsave(figfile, plot=p_all, height=5*8.3, width=11.7)

# list outliers
outlier_all <- names(unlist(outlier_list))
outlier_all <- limma::strsplit2(outlier_all, split="\\.")[,2]

## ----create violinplots for PCR genes only in both dataset----------------------------------------

# replace with proper use of facet_grid

ras_tmp <- ras_both_long_tbl
ras_tmp$experiment <- factor(ras_tmp$experiment, levels=c("tcga", "pcr"))

# order of comparisons is crucial!
# I could only figure out to use one group
ras_tmp$gene_experiment <- paste(ras_tmp$gene, ras_tmp$experiment, sep="_")
ras_tmp$gene_experiment <- factor(ras_tmp$gene_experiment,
            levels=paste(rep(levels(ras_tmp$gene), each=2), 
                     rep(c("tcga", "pcr"), times=5), sep="_"))

# # define facet names
# facet_names <- levels(ras_tmp$gene_experiment)
# names(facet_names)  <- levels(ras_tmp$gene_experiment)
# for (ii in 1:length(facet_names)){
#   facet_names[ii] <- gsub("_", " ", toupper(facet_names[ii]), perl="fixed")
# }
# 
# p <- ggplot(data=ras_tmp, aes(x=category_factor, y=expression, color=category_factor))+
#   # geom_violin(bw=0.3, trim=F)+
#   # geom_point(size=0.5, aes(x=category_factor, y=expression, color=category_factor),
#   #            position=position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2))+
#   geom_violin(trim=F) +
#   geom_jitter(size=0.5, width = 0.2) +
#   #geom_boxplot(width=0.2)+
#   ylab("") + xlab("")+
#   guides(color=guide_legend(title=""))+
#   scale_color_manual("category", 
#                     values = c("Normal" = "black",
#                                "Tumor_Stage_I"="cyan",
#                                "Tumor_Stage_II"="purple",
#                                "Tumor_Stage_III"="pink",
#                                "Tumor_Stage_IV"="red"),
#                     labels= c("Normal" = "Normal",
#                               "Tumor_Stage_I"="Stage I",
#                               "Tumor_Stage_II"="Stage II",
#                               "Tumor_Stage_III"="Stage III",
#                               "Tumor_Stage_IV"="Stage IV"))+
#   scale_x_discrete(labels= c("Normal" = "Normal",
#                              "Tumor_Stage_I"="Stage I",
#                              "Tumor_Stage_II"="Stage II",
#                              "Tumor_Stage_III"="Stage III",
#                              "Tumor_Stage_IV"="Stage IV")) + 
#   theme(
#     #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
#     axis.text.x = element_text(angle = 90),
#         legend.position="none",
#         axis.title=element_text(size=20),
#         panel.grid.minor=element_blank(),
#         panel.grid.major=element_blank(),
#         axis.text=element_text(size=16),
#         legend.text=element_text(size=16),
#         strip.text.x = element_text(size=12), # facet titles
#         panel.background=element_rect(fill="white", colour = "black"),
#         plot.background=element_rect(color="black"))+
#   #facet_grid(gene~experiment)
#   facet_wrap(~gene_experiment, ncol=2,
#              labeller = as_labeller(facet_names))
# plot(p)


# also calculate p-values - 
# format needed for stat_pvalue_manual: group1 | group2 | p | y.position | etc.
# we need a proxy for Cond_Biorep, because fill is based on this


stat.test <- ras_tmp %>%
  group_by(gene_experiment) %>%
  t_test(expression ~ category_factor) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
#stat.test
# keep Normal vs StageI only
stat.test <- stat.test[which(stat.test$group1 %in% c("Normal") &
                               stat.test$group2 %in% c("Tumor_Stage_I")),]

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "category_factor",
                  dodge = 0.8, step.increase = 0.06,)

# p_final <- p + stat_pvalue_manual(
#   stat.test,  label = "p.adj.signif", tip.length = 0.005
# )
# 
# plot(p_final)
# outfig <- file.path(figfolder, "violinplot_with_stars.pdf")
# ggsave(outfig, plot=p_final, width=7, height=14)


##------------------------- old solution for violinplots with stars-------------
# add gene and experiment to stat.test
stat.test$gene <- limma::strsplit2(stat.test$gene_experiment, split="_")[,1]
stat.test$experiment <- limma::strsplit2(stat.test$gene_experiment, split="_")[,2]
# keep Normal vs StageI only
stat.test <- stat.test[which(stat.test$group1 %in% c("Normal") &
                               stat.test$group2 %in% c("Tumor_Stage_I")),]

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "category_factor",
                  dodge = 0.8, step.increase = 0.06,)
# change y.position completely manually
stat.test$y.position <- 4.5

theme_tcga <- theme(axis.title.x=element_blank(),
                    axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=0.5),
                    plot.margin = unit(c(0.05, 0.05, 0.3, 0.05), "inches"),
                    legend.position="none",
                    axis.title=element_text(size=20),
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_blank(),
                    axis.text=element_text(size=16),
                    legend.text=element_text(size=16),
                    strip.text.x = element_text(size=12), # facet titles
                    panel.background=element_rect(fill="white", colour = "black"),
                    plot.background=element_rect(color="white"))
  
  
theme_pcr <- theme(axis.title.x=element_blank(),
                   axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=0.5),
                   plot.margin = unit(c(0.05, 0.05, 0.3, 0.05), "inches"),
                   legend.position="none",
                   axis.title=element_text(size=20),
                   panel.grid.minor=element_blank(),
                   panel.grid.major=element_blank(),
                   axis.text=element_text(size=16),
                   legend.text=element_text(size=16),
                   strip.text.x = element_text(size=12), # facet titles
                   panel.background=element_rect(fill="white", colour = "black"),
                   plot.background=element_rect(color="white"),
                   axis.title.y=element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank())
  



plotlist <- list()
for (rname in names(ras_both_scaled)){
  ras_tmp <- ras_both_scaled[[rname]][
    ras_both_scaled[[rname]]$category_factor %in% 
      c("Normal", "Tumor_Stage_I", "Tumor_Stage_II", 
        "Tumor_Stage_III", "Tumor_Stage_IV"),]
  ras_tmp$category_factor <- droplevels(ras_tmp$category_factor)
  plotlist[[rname]] <- list()
  for (gene in pcr_genes){
    # stagelm <- lm(as.formula(paste0(gene, " ~ category_factor")), data=ras_tmp)
    # stageav <- aov(stagelm)
    # summary(stageav)
    # tukeydf <- as_tibble(TukeyHSD(stageav, conf.level=.95)$category_factor, 
    #                      rownames="comparison")
    # tukeydf$p.signif <- stars.pval(tukeydf$`p adj`)
    # tukeydf$p.signif[tukeydf$p.signif==" "] <- "ns"
    # tukeydf$group1 <- limma::strsplit2(tukeydf$comparison, split="-")[,2]
    # tukeydf$group2 <- limma::strsplit2(tukeydf$comparison, split="-")[,1]
    
    # use adjusted t-test value
    
    # TO DO - p values together, new colors, fehér háttér, facet_grid csillagokkal
    p1 <- ggplot(ras_tmp,
                 aes(x=category_factor, y=.data[[gene]], 
                     colour = category_factor)) +
      geom_violin(trim=F) +
      geom_jitter(size=0.5) +
      geom_boxplot(width=0.2)+
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
      xlab("") + ylab(gene)+ guides(colour="none") + ylim(-6, 6)+
      scale_x_discrete(labels=c("Normal", "Stage_I", "Stage_II", 
                                "Stage_III", "Stage_IV"))
    if (rname=="tcga") p1 <- p1 + theme_tcga
    if (rname=="pcr") p1 <- p1 + theme_pcr
      
      
    
    #y.position <- 4.5
    # testtbl <- as_tibble(list("group1"="Normal", "group2"="Tumor",
    #                           "p"=stars.pval(tres$p.value), "y.position"=y.position))
    # testtbl <- tukeydf[1,]
    # colnames(testtbl)[colnames(testtbl)=="p.signif"] <- "p"
    # testtbl$y.position <- y.position
    # 
    # p2 <- p1+ stat_pvalue_manual(testtbl)
    # plotlist[[rname]][[gene]] <- p2
    
    p2 <- p1+ stat_pvalue_manual(
      stat.test[which(stat.test$gene==gene & stat.test$experiment==rname),],  
      label = "p.adj.signif", tip.length = 0.05)
    plotlist[[rname]][[gene]] <- p2
  }
}
p_tcga <- cowplot::plot_grid(plotlist=plotlist[["tcga"]], ncol=1)
p_pcr <- cowplot::plot_grid(plotlist=plotlist[["pcr"]], ncol=1)
p_all <- cowplot::plot_grid(p_tcga, p_pcr, ncol=2, 
                            rel_widths = c(1.25,1))
plot(p_all)
outfig <- file.path(figfolder, "violinplot_pcrgenes_stages_with_stars.pdf")
ggsave(outfig, plot=p_all, width=7, height=14)
outfig2 <- file.path(figfolder, "violinplot_pcrgenes_stages_with_stars.tiff")
ggsave(outfig2, plot=p_all, width=7, height=14, dpi=600)

# t-test
pvalue_list <- list()
plotlist <- list()
for (rname in names(ras_both_scaled)){
  #print(rname)
  pvalue_list[[rname]] <- list()
  plotlist[[rname]] <- list()
  for (gname in pcr_genes){
    tres <- t.test(ras_both_scaled[[rname]][[gname]][
          which(ras_both_scaled[[rname]]$sample_type_factor == "Normal")],
          ras_both_scaled[[rname]][[gname]][
            which(ras_both_scaled[[rname]]$sample_type_factor == "Tumor")])
    # print(gname)
    # print(tres)
    pvalue_list[[rname]][[gname]] <- tres$p.value
    
    # add stars to plot
    # p1 <- ggplot(data=ras_both_scaled[[rname]]) +
    #   geom_boxplot(mapping=aes(x=sample_type_factor, y=.data[[gname]]))+
    #   scale_x_discrete(labels=c("Normal","Tumor"))+
    #   theme(axis.title.x=element_blank(),
    #         axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=1),
    #         plot.margin = unit(c(0.05, 0.05, 0.3, 0.05), 
    #                            "inches"))+
    #   ylim(-6, 6)
    
    p1 <- ggplot(ras_both_scaled[[rname]],
                aes(x=sample_type_factor, y=.data[[gname]], 
                    colour = sample_type_factor)) +
      geom_violin(trim=F) +
      geom_jitter(size=0.5) +
      geom_boxplot(width=0.2)+
      scale_color_manual(values = c("Tumor" = "red", "Normal" = "black"))+
      xlab("") + ylab(gname)+ guides(colour="none") + ylim(-6, 6)+
      scale_x_discrete(labels=c("Normal","Tumor"))+
        theme(axis.title.x=element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=1),
              plot.margin = unit(c(0.05, 0.05, 0.3, 0.05), "inches"))
    
    
    ymax <- (median(ras_both_scaled[[rname]][[gname]]) + 
               2 * IQR(ras_both_scaled[[rname]][[gname]]))
    ymin <- (median(ras_both_scaled[[rname]][[gname]]) -
               2 * IQR(ras_both_scaled[[rname]][[gname]]))
    ystep <- (ymax-ymin)/10
    #y.position = c(ymax+0.1*ystep, ymax+1.1*ystep, ymax+2.1*ystep)
    #y.position <- ymax+0.1*ystep
    y.position <- 4.5
    
    testtbl <- as_tibble(list("group1"="Normal", "group2"="Tumor",
                      "p"=stars.pval(tres$p.value), "y.position"=y.position))
    testtbl$p[testtbl$p==" "] <- "ns"
                    # "p"=tres$p.value, "y.position"=y.position))
    print(rname)
    print(gname)
    print(testtbl$p)
    p2 <- p1+ stat_pvalue_manual(testtbl)
    plotlist[[rname]][[gname]] <- p2
  }
}

p_tcga <- cowplot::plot_grid(plotlist=plotlist[["tcga"]], ncol=1)
p_pcr <- cowplot::plot_grid(plotlist=plotlist[["pcr"]], ncol=1)
p_all <- cowplot::plot_grid(p_tcga, p_pcr, ncol=2)
outfig <- file.path(figfolder, "violinplot_with_stars.pdf")
ggsave(outfig, plot=p_all, width=7, height=14)


## ----correlation---------------------------------------------------------------------------------------------
# correlation plots with legend without "cor"
cfigfolder <- file.path(figfolder, "correlation")
if (!file.exists(cfigfolder)) dir.create(cfigfolder)

# TO DO - use new color, "Corr" felierat nélkül, külön ábrákon, cím nélkül (stage I is)

corr_list <- list()
plot_list <- list()
for (rname in names(ras_both_scaled)){
  cfigfolder <- file.path(figfolder, "correlation")
  if (!file.exists(cfigfolder)) dir.create(cfigfolder)
  ras_tmp <- ras_both_scaled[[rname]]
  # plot correlation among the variables
  # do it for all samples, tumor samples and normal samples
  # order always as would be for normal
  corr <- round(cor(ras_tmp[which(ras_tmp$sample_type_factor=="Normal"), pcr_genes], 
                    use="pairwise.complete.obs",
                    method = "pearson"), digits=1)
  # order by the "Normal samples clustering
  hc <- hclust(dist(corr), method="ward.D2")
  # now plot
  corr_list[[rname]] <- list()
  plot_list[[rname]] <- list()
  for (stype in c("Normal", "Tumor")){
    sel_samples <- 1:nrow(ras_tmp)
    if (stype=="Normal"){
      sel_samples <- which(ras_tmp$sample_type_factor=="Normal")
    }
    if (stype=="Tumor"){
      sel_samples <- which(ras_tmp$sample_type_factor=="Tumor")
    }
    corr <- round(cor(ras_tmp[sel_samples, pcr_genes], 
                      use="pairwise.complete.obs",
                      method = "pearson"), digits=2)
    # mean absolute correlation
    corr <- corr[hc$order, hc$order]
    figfile <- file.path(cfigfolder, 
                         paste("correlation_pcr_genes_only_", 
                               rname, "_", stype, "_hc_order.pdf", sep=""))
    p <- ggcorrplot(corr, lab=T, legend.title = "",
                    method="circle", type="upper",
                    hc.order = F, sig.level=0.1, lab_size = 3,
                    colors = c("#234387", "white", "#b81007"),
                    show.legend = F)
    ggsave(figfile, plot=p, height=3, width=2.7)
    # store in list
    corr_list[[rname]][[stype]] <- corr
    
    # and save once with legend as well
    p <- ggcorrplot(corr, lab=T, legend.title = "",
                    method="circle", type="upper",
                    hc.order = F, sig.level=0.1, lab_size = 3,
                    colors = c("#234387", "white", "#b81007"),
                    show.legend = T)
    figfile <- file.path(cfigfolder, 
                         paste("correlation_pcr_genes_only_", 
                               rname, "_", stype, "_hc_order_withlegend.pdf", sep=""))
    ggsave(figfile, plot=p, height=3, width=2.7)
    figfile <- file.path(cfigfolder, 
                         paste("correlation_pcr_genes_only_", 
                               rname, "_", stype, "_hc_order_withlegend.tiff", sep=""))
    ggsave(figfile, plot=p, height=3, width=2.7, dpi=600)
    
    # #also for the new type figure
    # corr <- round(corr, digits=2)
    # p <- ggcorrplot(corr[hc$order, hc$order], method="circle", type="upper",
    #                 hc.order = F, sig.level=0.1, lab=T, lab_size = 3)
    # plot_list[[rname]][[stype]] <- p
    # # also save each figure separately
    # outfig <- file.path(cfigfolder, 
    #                     paste("pcr_correlation_",  stype, ".pdf", sep=""))
    # ggsave(outfig, plot=p, height=3, width=2.7)
  } # loop for stype
} # loop for rname

# repeat for stages
corr_list <- list()
plot_list <- list()
for (rname in names(ras_both_scaled)){
  cfigfolder <- file.path(figfolder, "correlation")
  if (!file.exists(cfigfolder)) dir.create(cfigfolder)
  ras_tmp <- ras_both_scaled[[rname]]
  # plot correlation among the variables
  # do it for all samples, tumor samples and normal samples
  # order always as would be for normal
  corr <- round(cor(ras_tmp[which(ras_tmp$sample_type_factor=="Normal"), pcr_genes], 
                    use="pairwise.complete.obs",
                    method = "pearson"), digits=1)
  # order by the "Normal samples clustering
  hc <- hclust(dist(corr), method="ward.D2")
  # now plot
  corr_list[[rname]] <- list()
  plot_list[[rname]] <- list()
  for (ctype in c("Normal", "Tumor_Stage_I")){
    sel_samples <- 1:nrow(ras_tmp)
    if (ctype=="Normal"){
      sel_samples <- which(ras_tmp$category_factor==ctype)
    }
    if (ctype=="Tumor_Stage_I"){
      sel_samples <- which(ras_tmp$category_factor==ctype)
    }
    corr <- round(cor(ras_tmp[sel_samples, pcr_genes], 
                      use="pairwise.complete.obs",
                      method = "pearson"), digits=2)
    # mean absolute correlation
    corr <- corr[hc$order, hc$order]
    figfile <- file.path(cfigfolder, 
                         paste("correlation_pcr_genes_only_", 
                               rname, "_", ctype, "_hc_order.pdf", sep=""))
    p <- ggcorrplot(corr, lab=T, legend.title = "",
                    method="circle", type="upper",
                    hc.order = F, sig.level=0.1, lab_size = 3,
                    colors = c("#234387", "white", "#b81007"),
                    show.legend = F)
    ggsave(figfile, plot=p, height=3, width=2.7)
    # store in list
    corr_list[[rname]][[stype]] <- corr
    
    
    
    # and save once with legend as well
    p <- ggcorrplot(corr, lab=T, legend.title = "",
                    method="circle", type="upper",
                    hc.order = F, sig.level=0.1, lab_size = 3,
                    colors = c("#234387", "white", "#b81007"),
                    show.legend = T)
    figfile <- file.path(cfigfolder, 
                         paste("correlation_pcr_genes_only_", 
                               rname, "_", ctype, "_hc_order_withlegend.pdf", sep=""))
    ggsave(figfile, plot=p, height=3, width=2.7)
    figfile <- file.path(cfigfolder, 
                         paste("correlation_pcr_genes_only_", 
                               rname, "_", ctype, "_hc_order_withlegend.tiff", sep=""))
    ggsave(figfile, plot=p, height=3, width=2.7, dpi=600)
    
    # #also for the new type figure
    # corr <- round(corr, digits=2)
    # p <- ggcorrplot(corr[hc$order, hc$order], method="circle", type="upper",
    #                 hc.order = F, sig.level=0.1, lab=T, lab_size = 3)
    # plot_list[[rname]][[stype]] <- p
    # # also save each figure separately
    # outfig <- file.path(cfigfolder, 
    #                     paste("pcr_correlation_",  stype, ".pdf", sep=""))
    # ggsave(outfig, plot=p, height=3, width=2.7)
  } # loop for stype
} # loop for rname


## ----LDA normal or tumor-------------------------------------------------------------------------------------

rname <- "tcga"
ras_tmp <- ras_both_scaled[[rname]]
ri <- which(ras_tmp$category_factor %in% c("Normal", 
                                      "Tumor_Stage_I","Tumor_Stage_IV"))
ras_tmp <- ras_tmp[ri, ]
ras_tmp$category_factor <- droplevels(ras_tmp$category_factor)

# similarly for pcr data
ras_tmp2 <- ras_both_scaled[["pcr"]]
# ri <- which(ras_tmp2$category_factor %in% c("Normal", 
#                                            "Tumor_Stage_I","Tumor_Stage_IV"))
# ras_tmp2 <- ras_tmp2[ri, ]
# ras_tmp2$category_factor <- droplevels(ras_tmp2$category_factor)

indata_mx <-  as.matrix(ras_tmp[, pcr_genes])
lda_model <- MASS::lda(indata_mx, grouping=ras_tmp$category_factor)
lda_values <- predict(lda_model, newdata = ras_tmp[, pcr_genes])
all.equal(lda_values$class, ras_tmp$category_factor) # "40 string mismatches"
lda_annot <- as_tibble(lda_values$x)
lda_annot$sample <- ras_tmp$sample
lda_annot$category_factor <- ras_tmp$category_factor
lda_annot$experiment <- "tcga"

# obtain the “proportion of trace” which is the proportion of between-class 
# variance that is explained by successive discriminant functions. 
prop_trace <- prop.table(lda_model$svd^2)
nicelabel <- list()
for (ii in 1:length(prop_trace)){
  nicelabel[[ii]] <- paste('LD',
                           as.character(ii), ' (', 
                           as.character(round(100*prop_trace[ii], digits = 2)), '%)', sep="")
}


# plot results for TCGA data
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
outfig <- file.path(figfolder, rtype, "ras_lda_pcrgenes_only_normal_stageI_stageIV.pdf")
ggsave(outfig, plot=pfinal, width=7, height=5)
# outfig <- file.path(figfolder, rtype, "ras_lda_pcrgenes_only_normal_stageI_stageIV.tiff")
# ggsave(outfig, plot=pfinal, width=7, height=7, dpi=600)

# add prediction for pcr data
lda_pcr <- predict(lda_model, newdata=ras_tmp2[, pcr_genes])
lda_pcr <- as_tibble(lda_pcr$x)
lda_pcr$sample <- ras_tmp2$sample
lda_pcr$category_factor <- ras_tmp2$category_factor
lda_pcr$experiment <- "pcr"

xname <- "LD1"
yname <- "LD2"
pfinal <- ggplot(data=lda_pcr) +
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
  xlab(nicelabel[[1]]) + ylab(nicelabel[[2]])+
  theme(plot.background=element_rect(color="white"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=20),
        legend.text=element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background=element_rect(fill="white", colour = "black"))
plot(pfinal)
outfig <- file.path(figfolder, "ras_lda_tcga_space_pcr_samples.pdf")
ggsave(outfig, plot=pfinal, width=7, height=5)
outfig <- file.path(figfolder, "ras_lda_tcga_space_pcr_samples.tiff")
ggsave(outfig, plot=pfinal, width=7, height=5, dpi=600)

# plot LD1 density plot
lda_tmp <- lda_pcr

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
outfig <- file.path(figfolder,"ras_pcr_LD1_density_categories.pdf")
ggsave(outfig, plot=p, width=7, height=5)
outfig <- file.path(figfolder,"ras_pcr_LD1_density_categories.tiff")
ggsave(outfig, plot=p, width=7, height=5, dpi=600)

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
outfig <- file.path(figfolder, "ras_pcr_LD2_density_categories.pdf")
ggsave(outfig, plot=p, width=7, height=5)

# join
lda_annot <- bind_rows(lda_annot, lda_pcr)

# factors order
lda_annot$experiment <- factor(lda_annot$experiment, levels=c("tcga", "pcr"))

# obtain the “proportion of trace” which is the proportion of between-class 
# variance that is explained by successive discriminant functions. 
prop_trace <- prop.table(lda_model$svd^2)
nicelabel <- list()
for (ii in 1:length(prop_trace)){
  nicelabel[[ii]] <- paste('LD',
                           as.character(ii), ' (', 
                           as.character(round(100*prop_trace[ii], digits = 2)), '%)', sep="")
}

xname <- "LD1"
yname <- "LD2"
pfinal <- ggplot(data=lda_annot) +
  geom_point(mapping=aes(x = .data[[xname]], y = .data[[yname]], 
                         color=category_factor, shape = experiment),
             size = 3, stroke = 1, fill="white")+
  # color points as required
  scale_color_manual("category", 
                     values = c("Normal" = "black",
                                "Tumor_Stage_I"=stage1col,
                                "Tumor_Stage_IV"=stage4col),
                     labels= c("Normal" = "Normal",
                               "Tumor_Stage_I"="Stage I",
                               "Tumor_Stage_IV"="Stage IV"))+
  scale_shape_manual(values= c("tcga" = 1,
                                "pcr"=2),
                     labels= c("tcga" = "TCGA",
                               "pcr"="PCR"))+
  xlab(nicelabel[[1]]) + ylab(nicelabel[[2]])+
  theme(plot.background=element_rect(color="white"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=20),
        legend.text=element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background=element_rect(fill="white", colour = "black"))
plot(pfinal)
outfig <- file.path(figfolder, "ras_lda_tcga_plus_pcr_normal_stageI_stageIV.pdf")
ggsave(outfig, plot=pfinal, width=7, height=5)
  
  
## ----single gene balanced accuracy normal vs stageI----------------------------------------------------------
ras_tmp <- ras_data[ras_data$category2 %in% 
                      c("Normal", "Tumor_Stage_I"),
                    c(indep_vars, "category_factor")]
# also delete unused factors
ras_tmp$category_factor <- droplevels(ras_tmp$category_factor)

# calculate all measures
acc_file <- file.path(resfolder, "single_gene_model_measures_normal_vs_stageI.rds")
all_measures <- calculate_single_gene_measures(acc_file,
                                           indep_vars,
                                           ras_tmp,
                                           "category_factor",
                                           c("Normal", "Tumor_Stage_I"),
                                           "Tumor_Stage_I")

# write confusion matrix to file for selected genes
sel_genes <- c("RAF1", "PLCE1")
conffolder <- file.path(resfolder, "confusion_matrices")
if(!file.exists(conf_folder)){dir.create(conf_folder)}
for (sel_gene in sel_genes){
  it_model  <-  glm(as.formula(paste0("category_factor ~ ",
                                          paste(sel_gene,collapse="+"))),
                        data = ras_tmp, family = "binomial")
  conf_mx <- runCV_LOOCV(ras_tmp, it_model,
                  depvar="category_factor",
                  depvar_levels=c("Normal", "Tumor_Stage_I"),
                  posclass="Tumor_Stage_I",
                  numCV_inner=3,
                  return_conf_mx=T)
  outfile <- file.path(conffolder, 
                       paste("single_gene_normal_vs_stageI_", sel_gene, ".csv", sep=""))
  print_confmx(conf_mx, outfile)
}

# save to text file as well
xls_file <- gsub(".rds", ".xlsx", acc_file)
openxlsx::write.xlsx(all_measures[
  order(all_measures$balanced_accuracy, decreasing = T),], file=xls_file)

# # now start collecting the balanced accuracy values
# acc_file <- file.path(resfolder, "balanced_accuracy_normal_vs_stageI.rds")
# if (!file.exists(acc_file)){
#   acc_list <- list()
#   for (gene in indep_vars){
#     print(gene)
#     it_model  <-  glm(as.formula(paste0("category_factor ~ ",
#                                         paste(gene,collapse="+"))),
#                       data = ras_tmp, family = "binomial")
#     acc_list[[gene]] <- runCV_LOOCV(ras_tmp, it_model, verbose=F,
#                                     depvar="category_factor",
#                      depvar_levels=c("Normal", "Tumor_Stage_I"),
#                      posclass="Tumor_Stage_I",
#                           return_conf_mx=F)
#   }  
#   saveRDS(acc_list, file=acc_file)
# } else {
#   acc_list <- readRDS(acc_file)
# }

# # create a barplot
# acc_tbl <- as.data.frame(unlist(acc_list))
# colnames(acc_tbl) <- "balanced_accuracy"
# acc_tbl <- as_tibble(acc_tbl, rownames="gene")

# put stars based on t-test
pvalue_list <- list()
for (ivar in indep_vars){
  ttest_res <- (
    t.test(ras_tmp[[ivar]][ras_tmp$category_factor=="Normal"],
           ras_tmp[[ivar]][ras_tmp$category_factor=="Tumor_Stage_I"]))
  pvalue_list[[ivar]] <- ttest_res$p.value
}
pvaldf <- as.data.frame(unlist(pvalue_list))
colnames(pvaldf) <- "p.value" # calculate adjusted p value?
pvaltbl <- as_tibble(pvaldf, rownames="gene")

# join with previous results
all_measures <- left_join(all_measures, pvaltbl, by=join_by(gene))
all_measures$ttest_stars <- stars.pval(all_measures$p.value)
all_measures <- all_measures[order(all_measures$balanced_accuracy, decreasing = F),]
all_measures$gene <- factor(all_measures$gene, levels=all_measures$gene)

# add p-values and stars then plot with the accuracy values
p <- ggplot(all_measures, aes(x = balanced_accuracy, y = gene)) +
  geom_col(fill = "gray70") +
  ## add percentage labels
  geom_text(aes(label = ttest_stars), nudge_x = 0.05) +
  geom_vline(xintercept=0.9, color="red")+
  theme_minimal()+
  xlab("Balanced accuracy") + ylab("")
outfig <- file.path(figfolder, "single_gene_balanced_acc_normal_vs_stageI.pdf")
ggsave(outfig, plot=p, width=4, height=7)


## ----single gene balanced accuracy normal vs tumor-----------------------------------------------------------
ras_tmp <- ras_data[ras_data$sample_type2 %in% 
                      c("Normal", "Tumor"),
                    c(indep_vars, "sample_type_factor")]
# also delete unused factors
ras_tmp$sample_type_factor <- droplevels(ras_tmp$sample_type_factor)

# calculate all measures
acc_file <- file.path(resfolder, "single_gene_model_measures_normal_vs_tumor.rds")
all_measures <- calculate_single_gene_measures(acc_file,
                                           indep_vars,
                                           ras_tmp,
                                           "sample_type_factor",
                                           c("Normal", "Tumor"),
                                           "Tumor")

# write confusion matrix to file for selected genes
sel_genes <- c("PLCE1")
conffolder <- file.path(resfolder, "confusion_matrices")
if(!file.exists(conf_folder)){dir.create(conf_folder)}
for (sel_gene in sel_genes){
  it_model  <-  glm(as.formula(paste0("sample_type_factor ~ ",
                                          paste(sel_gene,collapse="+"))),
                        data = ras_tmp, family = "binomial")
  conf_mx <- runCV_LOOCV(ras_tmp, it_model,
                  depvar="sample_type_factor",
                  depvar_levels=c("Normal", "Tumor"),
                  posclass="Tumor",
                  numCV_inner=3,
                  return_conf_mx=T)
  outfile <- file.path(conffolder, 
                       paste("single_gene_normal_vs_tumor_", sel_gene, ".csv", sep=""))
  print_confmx(conf_mx, outfile)
}
  

# save to text file as well
xls_file <- gsub(".rds", ".xlsx", acc_file)
openxlsx::write.xlsx(all_measures[
  order(all_measures$balanced_accuracy, decreasing = T),], file=xls_file)

# put stars based on t-test
pvalue_list <- list()
for (ivar in indep_vars){
  ttest_res <- (
    t.test(ras_tmp[[ivar]][ras_tmp$sample_type_factor=="Normal"],
           ras_tmp[[ivar]][ras_tmp$sample_type_factor=="Tumor"]))
  pvalue_list[[ivar]] <- ttest_res$p.value
}
pvaldf <- as.data.frame(unlist(pvalue_list))
colnames(pvaldf) <- "p.value" # calculate adjusted p value?
pvaltbl <- as_tibble(pvaldf, rownames="gene")

all_measures <- left_join(all_measures, pvaltbl, by=join_by(gene))
all_measures$ttest_stars <- stars.pval(all_measures$p.value)
all_measures <- all_measures[order(all_measures$balanced_accuracy, decreasing = F),]
all_measures$gene <- factor(all_measures$gene, levels=all_measures$gene)

# add p-values and stars then plot with the accuracy values
p <- ggplot(all_measures, aes(x = balanced_accuracy, y = gene)) +
  geom_col(fill = "gray70") +
  ## add percentage labels
  geom_text(aes(label = ttest_stars), nudge_x = 0.05) +
  geom_vline(xintercept=0.9, color="red")+
  theme_minimal()+
  xlab("Balanced accuracy") + ylab("")
outfig <- file.path(figfolder, "single_gene_balanced_acc_normal_vs_tumor.pdf")
ggsave(outfig, plot=p, width=4, height=7)


## ----single gene balanced accuracy stage I vs stage IV-------------------------------------------------------
ras_tmp <- ras_data[ras_data$category2 %in% 
                      c("Tumor_Stage_I", "Tumor_Stage_IV"),
                    c(indep_vars, "category_factor")]
# also delete unused factors
ras_tmp$category_factor <- droplevels(ras_tmp$category_factor)

# calculate all measures
acc_file <- file.path(resfolder, "single_gene_model_measures_stageI_vs_stageVI.rds")
all_measures <- calculate_single_gene_measures(acc_file,
                                           indep_vars,
                                           ras_tmp,
                                           "category_factor",
                                           c("Tumor_Stage_I", "Tumor_Stage_IV"),
                                           "Tumor_Stage_IV")
# save to text file as well
xls_file <- gsub(".rds", ".xlsx", acc_file)
openxlsx::write.xlsx(all_measures[
  order(all_measures$balanced_accuracy, decreasing = T),], file=xls_file)

# put stars based on t-test
pvalue_list <- list()
for (ivar in indep_vars){
  ttest_res <- (
    t.test(ras_tmp[[ivar]][ras_tmp$category_factor=="Tumor_Stage_I"],
           ras_tmp[[ivar]][ras_tmp$category_factor=="Tumor_Stage_IV"]))
  pvalue_list[[ivar]] <- ttest_res$p.value
}
pvaldf <- as.data.frame(unlist(pvalue_list))
colnames(pvaldf) <- "p.value" # calculate adjusted p value?
pvaltbl <- as_tibble(pvaldf, rownames="gene")

all_measures <- left_join(all_measures, pvaltbl, by=join_by(gene))
all_measures$ttest_stars <- stars.pval(all_measures$p.value)
all_measures <- all_measures[order(all_measures$balanced_accuracy, decreasing = F),]
all_measures$gene <- factor(all_measures$gene, levels=all_measures$gene)

# add p-values and stars then plot with the accuracy values
p <- ggplot(all_measures, aes(x = balanced_accuracy, y = gene)) +
  geom_col(fill = "gray70") +
  ## add percentage labels
  geom_text(aes(label = ttest_stars), nudge_x = 0.05) +
  #geom_vline(xintercept=0.9, color="red")+
  theme_minimal()+
  xlab("Balanced accuracy") + ylab("")
outfig <- file.path(figfolder, "single_gene_balanced_acc_stageI_vs_stageIV.pdf")
ggsave(outfig, plot=p, width=4, height=7)





## ----k gene balanced accuracy accross all comparisons----------------------------------------------------------
comptype_vector <- c("normal_vs_tumor", "normal_vs_stageI", 
                     "stageI_vs_stageIV")
depvar_vector <- c("sample_type_factor", "category_factor", "category_factor")
depvar_levels_list <- list(c("Normal", "Tumor"),
                           c("Normal", "Tumor_Stage_I"),
                           c("Tumor_Stage_I", "Tumor_Stage_IV"))

for (comparison_index in c(2)){
  plotlist <- list()
  measures_list <- list()
  depvar <- depvar_vector[comparison_index]
  depvar_levels <- depvar_levels_list[[comparison_index]]
  ras_tmp <- ras_both_scaled[["pcr"]][ras_both_scaled[["pcr"]][[depvar]]
                                      %in%  depvar_levels,
                                      c(pcr_genes, depvar)]
  ras_tmp[[depvar]] <- droplevels(ras_tmp[[depvar]])
  comptype <- comptype_vector[comparison_index]
  # calculate all measures
  for (k in c(1:2)){
    acc_file <- file.path(resfolder, 
        paste("k", as.character(k), "_pcr_gene_model_measures_",
               comptype, ".rds", sep=""))
    all_measures <- calculate_k_gene_measures(acc_file,
                                              k,
                                              pcr_genes,
                                              ras_tmp,
                                              depvar,
                                              depvar_levels_list[[comparison_index]],
                                              depvar_levels_list[[comparison_index]][2])
    # change '.' to '+'
    all_measures$genes <- gsub("\\.", "+", all_measures$genes)
    all_measures <- all_measures[order(all_measures$balanced_accuracy, decreasing = T),]
    measures_list[[k]] <- all_measures
    
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
        #geom_col(aes(fill = barcolor), color="black", width = 0.8) +
        geom_col(fill = "white", color="black", width = 0.8) +
        # scale_fill_manual("category",
        #                   values = c("#d4d4d4" = "#d4d4d4",
        #                              "#CD3232"="#CD3232",
        #                              "#33ce33"="#33ce33"))+
        ## add percentage labels
        # geom_text(aes(label = genes, x=0.2), nudge_x = -0.05) +
        geom_text(aes(label = genes, x=0.02, color=barcolor), hjust=0) +
        scale_color_manual("category",
                           values = c("#d4d4d4" = "black",
                                      "#CD3232"="#CD3232",
                                      "#33ce33"="#33ce33"))+
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
        scale_x_continuous(limits=c(0,1.05), breaks= seq(0, 1, by = 0.2))+
        xlab("Balanced accuracy")
    plot(plotlist[[k]])
    } else { # TO DO for k=2: 0.95 felett; bar-on kívülre: balanced accuracy számmal
      # multicolor text: https://andrewwhitby.com/2017/09/18/multi-color-text-ggplot2/
      # stage1 vs stage4: for k=3-5: first 10 bars
      k1data <- indata
      
      #indata <- indata[1:length(indep_vars), ]
      indata <- all_measures
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
    # # plot 
    # indata <- all_measures
    # # order for plot
    # indata <- indata[order(indata$balanced_accuracy, decreasing = F),]
    # indata$genes <- factor(indata$genes, levels=indata$genes)
    # 
    # plotlist[[k]] <- ggplot(indata, aes(x = balanced_accuracy, y = genes)) +
    #   geom_col(fill = "gray70") +
    #   ## add percentage labels
    #   geom_text(aes(label = genes, x=0.2), nudge_x = -0.05) +
    #   # geom_vline(xintercept=0.95, color="red")+
    #   theme(axis.title.y=element_blank(),
    #         axis.text.y=element_blank(),
    #         axis.ticks.y=element_blank())+
    #   xlab("Balanced accuracy") + xlim(0,1)
    
  # save to text file as well
  xls_file <-  file.path(resfolder, 
                         paste("k", "_pcr_gene_model_measures_",
                               comptype, ".xls", sep=""))
  openxlsx::write.xlsx(measures_list, file=xls_file)
  # save the figures together
  p <- cowplot::plot_grid(plotlist=plotlist)
  figfile <- file.path(figfolder, paste("k", "_pcr_gene_model_measures_",
                                        comptype, ".pdf", sep=""))
  ggsave(figfile, plot=p, width=10, height=7)
  }
}

##----confusion matrices------------------------------------

comptype_vector <- c("normal_vs_tumor", "normal_vs_stageI", 
                     "stageI_vs_stageIV")
depvar_vector <- c("sample_type_factor", "category_factor", "category_factor")
depvar_levels_list <- list(c("Normal", "Tumor"),
                           c("Normal", "Tumor_Stage_I"),
                           c("Tumor_Stage_I", "Tumor_Stage_IV"))
ci <- 2
for (rname in names(ras_both_scaled)){
  ras_tmp <- ras_both_scaled[[rname]][ras_both_scaled[[rname]][[depvar]]
                                      %in%  depvar_levels,
                                      c(pcr_genes, depvar)]
  ras_tmp[[depvar]] <- droplevels(ras_tmp[[depvar]])
  ras_data$category_factor <- droplevels(ras_data$category_factor)
  it_model <-  glm(as.formula(paste0("`", depvar_vector[ci], "`" ," ~ ",
                                     paste(pcr_genes,collapse="+"))),
                   data = ras_tmp, family = "binomial")
  
  conf_mx <- runCV_LOOCV(ras_tmp, it_model,
                  depvar=depvar_vector[ci],
                  depvar_levels=depvar_levels_list[[ci]],
                  posclass=depvar_levels_list[[ci]][2],
                  return_conf_mx=T)
  print(rname)
  print(conf_mx$table)
  print("##############")
  outfile <- file.path(resfolder,
      paste(rname, "model_pcr_genes_", comptype_vector[ci], ".csv", sep=""))
  print_confmx(conf_mx, outfile)
}


## ----double gene balanced accuracy stageI vs stageIV---------------------------------------------------------
ras_tmp <- ras_data[ras_data$category2 %in% 
                      c("Tumor_Stage_I", "Tumor_Stage_IV"),
                    c(indep_vars, "category_factor")]
# also delete unused factors
ras_tmp$category_factor <- droplevels(ras_tmp$category_factor)

# calculate all measures
acc_file <- file.path(resfolder, "two_gene_model_measures_stageI_vs_stageIV.rds")
all_measures <- calculate_double_gene_measures(acc_file,
                                           indep_vars,
                                           ras_tmp,
                                           "category_factor",
                                           c("Tumor_Stage_I", "Tumor_Stage_IV"),
                                           "Tumor_Stage_IV")
# save to text file as well
xls_file <- gsub(".rds", ".xlsx", acc_file)
# change '.' to '+'
all_measures$gene_pair <- gsub("\\.", "+", all_measures$gene_pair)
all_measures$gene1 <- limma::strsplit2(all_measures$gene_pair, split="\\+")[,1]
all_measures$gene2 <- limma::strsplit2(all_measures$gene_pair, split="\\+")[,2]
openxlsx::write.xlsx(all_measures[
  order(all_measures$balanced_accuracy, decreasing = T),], file=xls_file)

# # write confusion matrix to file for selected genes
# sel_genes <- c("RGL1+PLCE1", 
# "PLCE1+MYO10", "PIK3C2B+PLCE1", "RAF1+PLCE1", "PLCE1+RADIL", "PLCE1+GRB7")
# conffolder <- file.path(resfolder, "confusion_matrices")
# if(!file.exists(conf_folder)){dir.create(conf_folder)}
# for (sel_gene in sel_genes){
#   it_model  <-  glm(as.formula(paste0("category_factor ~ ",
#                                           sel_gene)),
#                         data = ras_tmp, family = "binomial")
#   conf_mx <- runCV_LOOCV(ras_tmp, it_model,
#                   depvar="category_factor",
#                   depvar_levels=c("Tumor_Stage_I", "Tumor_Stage_IV"),
#                   posclass="Tumor_Stage_IV",
#                   numCV_inner=3,
#                   return_conf_mx=T)
#   outfile <- file.path(conffolder, 
#                 paste("double_genes_stageI_vs_stageIV_", sel_gene, ".csv", sep=""))
#   print_confmx(conf_mx, outfile)
# }
  

# plot gene pairs with accuracy above 0.60
indata <- all_measures[which(all_measures$balanced_accuracy>=0.6),]
# order for plot
indata <- indata[order(indata$balanced_accuracy, decreasing = F),]
indata$gene_pair <- factor(indata$gene_pair, levels=indata$gene_pair)

p <- ggplot(indata, aes(x = balanced_accuracy, y = gene_pair)) +
  geom_col(fill = "gray70") +
  ## add percentage labels
  geom_text(aes(label = gene_pair, x=0.2), nudge_x = -0.05) +
  # geom_vline(xintercept=0.95, color="red")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  xlab("Balanced accuracy")
outfig <- file.path(figfolder, "double_gene_balanced_acc_stageI_vs_stageIV_thr06.pdf")
ggsave(outfig, plot=p, width = 7, height=14)


## ----double gene balanced accuracy normal vs tumor-----------------------------------------------------------
ras_tmp <- ras_data[ras_data$sample_type2 %in% 
                      c("Normal", "Tumor"),
                    c(indep_vars, "sample_type_factor")]
# also delete unused factors
ras_tmp$sample_type_factor <- droplevels(ras_tmp$sample_type_factor)

# calculate all measures
acc_file <- file.path(resfolder, "two_gene_model_measures_normal_vs_tumor.rds")
all_measures <- calculate_double_gene_measures(acc_file,
                                           indep_vars,
                                           ras_tmp,
                                           "sample_type_factor",
                                           c("Normal", "Tumor"),
                                           "Tumor")
# save to text file as well
xls_file <- gsub(".rds", ".xlsx", acc_file)
# change '.' to '+'
all_measures$gene_pair <- gsub("\\.", "+", all_measures$gene_pair)
all_measures$gene1 <- limma::strsplit2(all_measures$gene_pair, split="\\+")[,1]
all_measures$gene2 <- limma::strsplit2(all_measures$gene_pair, split="\\+")[,2]
openxlsx::write.xlsx(all_measures[
  order(all_measures$balanced_accuracy, decreasing = T),], file=xls_file)

# # write confusion matrix to file for selected genes
# sel_genes <- c("RGL1+PLCE1", 
# "PLCE1+MYO10", "PIK3C2B+PLCE1", "RAF1+PLCE1", "PLCE1+RADIL", "PLCE1+GRB7")
# conffolder <- file.path(resfolder, "confusion_matrices")
# if(!file.exists(conf_folder)){dir.create(conf_folder)}
# for (sel_gene in sel_genes){
#   it_model  <-  glm(as.formula(paste0("sample_type_factor ~ ",
#                                           sel_gene)),
#                         data = ras_tmp, family = "binomial")
#   conf_mx <- runCV_LOOCV(ras_tmp, it_model,
#                   depvar="sample_type_factor",
#                   depvar_levels=c("Normal", "Tumor"),
#                   posclass="Tumor",
#                   numCV_inner=3,
#                   return_conf_mx=T)
#   outfile <- file.path(conffolder, 
#                 paste("double_genes_normal_vs_tumor_", sel_gene, ".csv", sep=""))
#   print_confmx(conf_mx, outfile)
# }
  

# plot gene pairs with accuracy above 0.95
indata <- all_measures[which(all_measures$balanced_accuracy>=0.9),]
# order for plot
indata <- indata[order(indata$balanced_accuracy, decreasing = F),]
indata$gene_pair <- factor(indata$gene_pair, levels=indata$gene_pair)

p <- ggplot(indata, aes(x = balanced_accuracy, y = gene_pair)) +
  geom_col(fill = "gray70") +
  ## add percentage labels
  geom_text(aes(label = gene_pair, x=0.2), nudge_x = -0.05) +
  # geom_vline(xintercept=0.95, color="red")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  xlab("Balanced accuracy")
outfig <- file.path(figfolder, "double_gene_balanced_acc_normal_vs_tumor_thr09.pdf")
ggsave(outfig, plot=p, width = 7, height=14)


## ----three and four genes balanced accuracy at different comparisons: results from server----------------------------------------------------------
num_genes_vector <- c(3:5)
comptype_vector <- c("normal_vs_tumor", "normal_vs_stageI", 
                     "stageI_vs_stageIV")

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
    
    # save to text file as well
    names(all_measures) <- c("AIC", "BIC", "balanced_accuracy", "gene_trio")
    all_measures <- all_measures[c("gene_trio", "AIC", "BIC", "balanced_accuracy")]
    xls_file <- gsub(".rds", ".xlsx", acc_file)
    all_measures <- all_measures[
      order(all_measures$balanced_accuracy, decreasing = T),]
    openxlsx::write.xlsx(all_measures, file=xls_file)
    
    # for the plot, reverse order is needed
    all_measures <- all_measures[
      order(all_measures$balanced_accuracy, decreasing = F),]
    all_measures$gene_trio <- factor(all_measures$gene_trio, 
                                     levels=all_measures$gene_trio)
    p <- ggplot(all_measures, aes(x = balanced_accuracy, y = gene_trio)) +
      geom_col(fill = "gray70") +
      ## add percentage labels
      geom_text(aes(label = gene_trio, x=0.1*num_genes), nudge_x = -0.05) +
      # geom_vline(xintercept=0.95, color="red")+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())+
      xlim(0,1)+
      xlab("Balanced accuracy")
    outfig <- file.path(figfolder, 
                        paste(rtype, "best_", as.character(num_genes), 
                              "gene_models_", comptype,".pdf", sep=""))
    ggsave(outfig, plot=p, width=4, height=14)
  }
}



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

