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
  
## ----k gene balanced accuracy accross all comparisons----------------------------------------------------------
comptype_vector <- c("normal_vs_tumor", "normal_vs_stageI", 
                     "stageI_vs_stageIV")
depvar_vector <- c("sample_type_factor", "category_factor", "category_factor")
depvar_levels_list <- list(c("Normal", "Tumor"),
                           c("Normal", "Tumor_Stage_I"),
                           c("Tumor_Stage_I", "Tumor_Stage_IV"))

for (comparison_index in c(1:2)){
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