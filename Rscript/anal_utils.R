# 2024.10.01.
# Anna Lovrics
# collection of utility functions

######################x
runCV <- function(indata, test_model,
                  depvar="sample_type_factor",
                  depvar_levels=c("Normal", "Tumor"),
                  posclass="Tumor",
                  numCV_outer=5,
                  numCV_inner=3,
                  verbose=F,
                  return_conf_mx=F){
  hold_rows_list <- createFolds(indata[[depvar]],
                                k=numCV_outer, list = T)
  logreg_model_best_list <- list()
  logreg_model_best_resp <- factor(rep(posclass, times=nrow(indata)), 
                                   levels=depvar_levels)
  for (cvi in 1:numCV_outer){
    if (verbose) print(paste("cvi: ", as.character(cvi), sep=""))
    test_rows <- hold_rows_list[[cvi]]
    ras_train <- indata[-test_rows,]
    ras_test <- indata[test_rows,]
    # create the sampling data for the inner groups as well
    inner_train_list <- createFolds(ras_train[[depvar]], 
                                    k=numCV_inner, list = T, returnTrain=T)
    ctrlspecs <- trainControl(method="CV",
                              index=inner_train_list, 
                              savePredictions="all",
                              classProbs=TRUE)
    
    suppressWarnings(logreg_model_best <- train(
      #eval(log_best_list[[1]]$formula),
      as.formula(paste(depvar, " ~ ",
                       as.character(test_model$formula[3]))),
      data=ras_train,
      method="glm",
      preProcess=NULL,
      family=binomial,
      trControl=ctrlspecs))
    # save model
    logreg_model_best_list[[cvi]] <- logreg_model_best
    # predict response
    logreg_model_best_resp[test_rows] <- predict(logreg_model_best, newdata=ras_test)
  }
  # calculate confusion matrix
  conf_mx <- confusionMatrix(data=logreg_model_best_resp,
                             reference=indata[[depvar]],
                             positive=posclass)
  # also print confusion matrix of requested
  if (verbose==T) print(conf_mx)
  # use balanced accuracy to measure goodness
  if (return_conf_mx==T){
    res_val <- conf_mx
  } else {
    res_val <- conf_mx$byClass[["Balanced Accuracy"]]
  }
  return(res_val)
}

gen.next.cbn <- function(cbn, n){
  ## Generates the combination that follows the one provided as input
  cbn.bin      <- rep(0, n)
  cbn.bin[cbn] <- 1
  if (tail(cbn.bin, 1) == 0){
    ind <- tail(which(cbn.bin == 1), 1)
    cbn.bin[c(ind, ind+1)] <- c(0, 1)
  }else{
    ind <- 1 + tail(which(diff(cbn.bin) == -1), 1)
    nb  <- sum(cbn.bin[-c(1:ind)] == 1)
    cbn.bin[c(ind-1, (n-nb+1):n)] <- 0
    cbn.bin[ind:(ind+nb)]         <- 1
  }
  cbn <- which(cbn.bin == 1)
}

########################
runCV_LOOCV_archive <- function(indata, test_model,
                  depvar="sample_type_factor",
                  depvar_levels=c("Normal", "Tumor"),
                  posclass="Tumor",
                  numCV_inner=3,
                  verbose=F,
                  return_conf_mx=F){
  set.seed(7)
  logreg_model_best_resp <- rep(posclass, times=nrow(indata))
  logreg_model_best_resp <- factor(logreg_model_best_resp,
                                   levels=depvar_levels)
  coefficients_list <- list()
  for (ii in 1:nrow(indata)){
    if (verbose) print(paste("ii: ", as.character(ii), 
                             "/", as.character(nrow(indata)), sep=""))
    ras_train <- indata[-ii,, drop=F]
    ras_test <- indata[ii,, drop=F]
    # # create the sampling data for the inner groups as well
    # inner_train_list <- createFolds(ras_train[[depvar]], 
    #                           k=numCV_inner, list = T, returnTrain=T)
    # ctrlspecs <- trainControl(method="CV",
    #                           index=inner_train_list, 
    #                           savePredictions=F,
    #                           classProbs=TRUE)
    # create the sampling data for the inner groups as well
    ctrlspecs <- trainControl(method="LOOCV",
                              savePredictions=F,
                              classProbs=TRUE)
    
    logreg_model_best <- train(
      #eval(log_best_list[[1]]$formula),
      as.formula(paste(depvar, " ~ ",
                       as.character(test_model$formula[3]))),
      data=ras_train,
      method="glm",
      preProcess=NULL,
      family="binomial",
      trControl=ctrlspecs)
    
    # Extract the model coefficients
    final_model <- model$finalModel  # The trained glm object
    coefficients_list[[ii]] <- coef(final_model)  # Store coefficients

    # predict response
    logreg_model_best_resp[ii] <- predict(logreg_model_best, newdata=ras_test)
    if (verbose) print(logreg_model_best_resp[ii])
  }
  if (verbose) print(coefficients_list)
  # calculate confusion matrix
  conf_mx <- confusionMatrix(data=logreg_model_best_resp,
                             reference=indata[[depvar]],
                             positive=posclass)
  # also print confusion matrix of requested
  if (verbose==T) print(conf_mx)
  # use balanced accuracy to measure goodness
  if (return_conf_mx==T){
    res_val <- conf_mx
  } else {
    res_val <- conf_mx$byClass[["Balanced Accuracy"]]
  }
  return(res_val)
}
#########################x
runCV_LOOCV <- function(indata, test_model,
                                depvar="sample_type_factor",
                                depvar_levels=c("Normal", "Tumor"),
                                posclass="Tumor",
                                numCV_inner=3,
                                verbose=F,
                                return_conf_mx=F){
  
  negclass <- setdiff(depvar_levels, posclass)
  logreg_model_best_resp <- rep(posclass, times=nrow(indata))
  logreg_model_best_resp <- factor(logreg_model_best_resp,
                                   levels=depvar_levels)
  #coefficients_list <- list()
  for (ii in 1:nrow(indata)){
    if (verbose) print(paste("ii: ", as.character(ii), 
                             "/", as.character(nrow(indata)), sep=""))
    ras_train <- indata[-ii,, drop=F]
    ras_test <- indata[ii,, drop=F]
    # logreg_model_best <- train(
    #   as.formula(paste(depvar, " ~ ",
    #                    as.character(test_model$formula[3]))),
    #   data=ras_train,
    #   method="glm",
    #   preProcess=NULL,
    #   family="binomial")
    # # Extract the model coefficients
    # final_model <- logreg_model_best$finalModel  # The trained glm object
    # coefficients_list[[ii]] <- coef(final_model)  # Store coefficients

    # predict response
    final_model <- glm(as.formula(paste(depvar, " ~ ",
                   as.character(test_model$formula[3]))), 
                   data = ras_train, family = "binomial")
    predicted_prob <- predict.glm(final_model, 
                            newdata=ras_test, type = "response")
    predicted_class <- ifelse(predicted_prob > 0.5, posclass, negclass)
    predicted_class <- factor(predicted_class, levels=depvar_levels)
    logreg_model_best_resp[ii] <- predicted_class
    if (verbose) print(logreg_model_best_resp[ii])
  }
  #if (verbose) print(coefficients_list)
  if (verbose) print("#---------------#")
  if (verbose) print(logreg_model_best_resp)
  # calculate confusion matrix
  conf_mx <- confusionMatrix(data=logreg_model_best_resp,
                             reference=indata[[depvar]],
                             positive=posclass)
  # also print confusion matrix of requested
  if (verbose==T) print(conf_mx)
  # use balanced accuracy to measure goodness
  if (return_conf_mx==T){
    res_val <- conf_mx
  } else {
    res_val <- conf_mx$byClass[["Balanced Accuracy"]]
  }
  return(res_val)
}

###########################

runCV_LOOCV_kNN <- function(indata, 
                        depvar="sample_type_factor",
                        depvar_levels=c("Normal", "Tumor"),
                        posclass="Tumor",
                        kval=5,
                        verbose=F,
                        return_conf_mx=F){
  
  negclass <- setdiff(depvar_levels, posclass)
  kNN_resp <- rep(posclass, times=nrow(indata))
  kNN_resp <- factor(kNN_resp,levels=depvar_levels)
  #coefficients_list <- list()
  for (ii in 1:nrow(indata)){
    if (verbose) print(paste("ii: ", as.character(ii), 
                             "/", as.character(nrow(indata)), sep=""))
    ras_train <- indata[-ii, 2:(ncol(indata)-1), drop=F]
    ras_test <- indata[ii, 2:(ncol(indata)-1), drop=F]
    train_class <- indata[[depvar]][-ii]
    predicted_class <- class::knn(train=ras_train, test=ras_test, 
        cl=train_class, k = kval, 
        l = 0, prob = FALSE, use.all = TRUE)
    kNN_resp[ii] <- predicted_class
    if (verbose) print(kNN_resp[ii])
  }
  if (verbose) print("#---------------#")
  if (verbose) print(kNN_resp)
  # calculate confusion matrix
  conf_mx <- confusionMatrix(data=kNN_resp,
                             reference=indata[[depvar]],
                             positive=posclass)
  # also print confusion matrix of requested
  if (verbose==T) print(conf_mx)
  # use balanced accuracy to measure goodness
  if (return_conf_mx==T){
    res_val <- conf_mx
  } else {
    res_val <- conf_mx$byClass[["Balanced Accuracy"]]
  }
  return(res_val)
}

#########################
gen.next.cbn <- function(cbn, n){
  ## Generates the combination that follows the one provided as input
  cbn.bin      <- rep(0, n)
  cbn.bin[cbn] <- 1
  if (tail(cbn.bin, 1) == 0){
    ind <- tail(which(cbn.bin == 1), 1)
    cbn.bin[c(ind, ind+1)] <- c(0, 1)
  }else{
    ind <- 1 + tail(which(diff(cbn.bin) == -1), 1)
    nb  <- sum(cbn.bin[-c(1:ind)] == 1)
    cbn.bin[c(ind-1, (n-nb+1):n)] <- 0
    cbn.bin[ind:(ind+nb)]         <- 1
  }
  cbn <- which(cbn.bin == 1)
}

######################
# now start collecting the balanced accuracy values
calculate_single_gene_measures <- function(acc_file,
                                           indep_vars,
                                           indata,
                                           depvar,
                                           depvar_levels,
                                           posclass){
  if (!file.exists(acc_file)){
    acc_list <- list()
    aic_list <- list()
    bic_list <- list()
    for (gene in indep_vars){
      print(gene)
      it_model  <-  glm(as.formula(paste0("`", depvar, "`" ," ~ ",
                                          paste(gene,collapse="+"))),
                        data = indata, family = "binomial")
      acc_list[[gene]] <- runCV_LOOCV(indata, it_model, verbose=F,
                                      depvar=depvar,
                                      depvar_levels=depvar_levels,
                                      posclass=posclass,
                                      return_conf_mx=F)
      aic_list[[gene]] <- it_model$aic
      bic_list[[gene]] <- BIC(it_model)
    } 
    all_measures <- t(rbind(as.data.frame(acc_list),
                            as.data.frame(aic_list),
                            as.data.frame(bic_list)))
    colnames(all_measures) <- c("balanced_accuracy", "aic", "bic")
    all_measures <- as_tibble(all_measures, rownames = "gene")
    saveRDS(all_measures, file=acc_file)
  } else {
    all_measures <- readRDS(acc_file)
  }
  return(all_measures)
}

##############x

calculate_double_gene_measures <- function(acc_file,
                                           indep_vars,
                                           indata,
                                           depvar,
                                           depvar_levels,
                                           posclass){
  if (!file.exists(acc_file)){
    acc_list <- list()
    aic_list <- list()
    bic_list <- list()
    for (ii in 1:(length(indep_vars)-1)){
      print(ii)
      for (jj in (ii+1):length(indep_vars)){
        test_genes <- c(indep_vars[ii], indep_vars[jj])
        lname <- paste(test_genes,collapse="+")
        it_model  <-  glm(as.formula(paste0("`", depvar, "`" ," ~ ",
                           paste(test_genes,collapse="+"))),
                           data = indata, family = "binomial")
        acc_list[[lname]] <- runCV_LOOCV(indata, it_model, verbose=F,
                                        depvar=depvar,
                                        depvar_levels=depvar_levels,
                                        posclass=posclass,
                                        return_conf_mx=F)
        aic_list[[lname]] <- it_model$aic
        bic_list[[lname]] <- BIC(it_model)
      }
    }
    all_measures <- t(rbind(as.data.frame(acc_list),
                            as.data.frame(aic_list),
                            as.data.frame(bic_list)))
    colnames(all_measures) <- c("balanced_accuracy", "aic", "bic")
    all_measures <- as_tibble(all_measures, rownames = "gene_pair")
    saveRDS(all_measures, file=acc_file)
  } else {
    all_measures <- readRDS(acc_file)
  }
  return(all_measures)
}

##############x

calculate_k_gene_measures <- function(acc_file,
                                           k,
                                           indep_vars,
                                           indata,
                                           depvar,
                                           depvar_levels,
                                           posclass){
  n <- length(indep_vars)
  if (!file.exists(acc_file)){
    acc_list <- list()
    aic_list <- list()
    bic_list <- list()
  for (i in 1:choose(n, k)){
    if (i == 1){
      cbn <- 1:k
    } else{
      cbn <- gen.next.cbn(cbn, n)
    }
    lname <- as.character(paste(indep_vars[cbn],collapse="+"))
    it_model <-  glm(as.formula(paste0("`", depvar, "`" ," ~ ",
                                       paste(indep_vars[cbn],collapse="+"))),
                     data = indata, family = "binomial")
    aic_list[[lname]] <- it_model$aic
    bic_list[[lname]] <- BIC(it_model)
    acc_list[[lname]] <- runCV_LOOCV(indata, it_model,
                          depvar=depvar,
                          depvar_levels=depvar_levels,
                          posclass=posclass,
                          return_conf_mx=F)
    }
    all_measures <- t(rbind(as.data.frame(acc_list),
                            as.data.frame(aic_list),
                            as.data.frame(bic_list)))
    colnames(all_measures) <- c("balanced_accuracy", "aic", "bic")
    all_measures <- as_tibble(all_measures, rownames = "genes")
    saveRDS(all_measures, file=acc_file)
  } else {
    all_measures <- readRDS(acc_file)
  }
  return(all_measures)
}

###############################################################################

print_confmx <- function(conf_mx, outfile, coefs=NULL){
  toappend <- F
  if (!is.null(coefs)){
    write.table(names(coefs), file=outfile, sep=",", row.names=F, col.names=F)
    cat("\n", file=outfile, append=TRUE)
    toappend <- T
  }
  write.table(as.table(conf_mx), file=outfile, sep=",", 
              row.names=F, col.names=F, append=toappend)
  cat("\n", file=outfile, append=TRUE)
  write.table(as.data.frame(as.table(conf_mx)), file=outfile, sep=",", 
              row.names = F, append=T)
  metrices <- as.data.frame(cbind(t(conf_mx$overall),t(conf_mx$byClass)))
  metrices[, !colnames(metrices) %in% c("AccuracyPValue", "McnemarPValue")] <- 
    round(metrices[, !colnames(metrices) %in% 
                     c("AccuracyPValue", "McnemarPValue")], digits=2)
  metrices <- metrices[, c(ncol(metrices), 1:(ncol(metrices)-1))]
  cat("\n", file=outfile, append=TRUE)
  write.table(metrices, outfile, sep = ",", append = T, row.names = F)
}

significant_correlation_threshold <- function(n, alpha = 0.05) {
  df <- n - 2  # Degrees of freedom
  t_critical <- qt(1 - alpha / 2, df)  # Two-tailed t critical value
  r_threshold <- sqrt(t_critical^2 / (t_critical^2 + df))  # Convert t to r
  return(r_threshold)
}

# define my own theme for the plots

theme_la <- function(
  ...,
  line,
  rect,
  text,
  title,
  aspect.ratio,
  axis.title=element_text(size=20),
  axis.title.x,
  axis.title.x.top,
  axis.title.x.bottom,
  axis.title.y,
  axis.title.y.left,
  axis.title.y.right,
  axis.text.x,
  axis.text.x.top,
  axis.text.x.bottom,
  axis.text.y,
  axis.text.y.left,
  axis.text.y.right,
  axis.text.theta,
  axis.text.r,
  axis.ticks,
  axis.ticks.x,
  axis.ticks.x.top,
  axis.ticks.x.bottom,
  axis.ticks.y,
  axis.ticks.y.left,
  axis.ticks.y.right,
  axis.ticks.theta,
  axis.ticks.r,
  axis.minor.ticks.x.top,
  axis.minor.ticks.x.bottom,
  axis.minor.ticks.y.left,
  axis.minor.ticks.y.right,
  axis.minor.ticks.theta,
  axis.minor.ticks.r,
  axis.ticks.length,
  axis.ticks.length.x,
  axis.ticks.length.x.top,
  axis.ticks.length.x.bottom,
  axis.ticks.length.y,
  axis.ticks.length.y.left,
  axis.ticks.length.y.right,
  axis.ticks.length.theta,
  axis.ticks.length.r,
  axis.minor.ticks.length,
  axis.minor.ticks.length.x,
  axis.minor.ticks.length.x.top,
  axis.minor.ticks.length.x.bottom,
  axis.minor.ticks.length.y,
  axis.minor.ticks.length.y.left,
  axis.minor.ticks.length.y.right,
  axis.minor.ticks.length.theta,
  axis.minor.ticks.length.r,
  axis.line,
  axis.line.x,
  axis.line.x.top,
  axis.line.x.bottom,
  axis.line.y,
  axis.line.y.left,
  axis.line.y.right,
  axis.line.theta,
  axis.line.r,
  legend.background,
  legend.margin,
  legend.spacing,
  legend.spacing.x,
  legend.spacing.y,
  legend.key,
  legend.key.size,
  legend.key.height,
  legend.key.width,
  legend.key.spacing,
  legend.key.spacing.x,
  legend.key.spacing.y,
  legend.frame,
  legend.ticks,
  legend.ticks.length,
  legend.axis.line,
  legend.text.position,
  legend.title,
  legend.title.position,
  legend.position,
  legend.position.inside,
  legend.direction,
  legend.byrow,
  legend.justification,
  legend.justification.top,
  legend.justification.bottom,
  legend.justification.left,
  legend.justification.right,
  legend.justification.inside,
  legend.location,
  legend.box,
  legend.box.just,
  legend.box.margin,
  legend.box.background,
  legend.box.spacing,
  panel.border,
  panel.spacing,
  panel.spacing.x,
  panel.spacing.y,
  panel.grid,
  panel.grid.major.x,
  panel.grid.major.y,
  panel.grid.minor.x,
  panel.grid.minor.y,
  panel.ontop,
  plot.title,
  plot.title.position,
  plot.subtitle,
  plot.caption,
  plot.caption.position,
  plot.tag,
  plot.tag.position,
  plot.tag.location,
  plot.margin,
  strip.background,
  strip.background.x,
  strip.background.y,
  strip.clip,
  strip.placement,
  strip.text,
  strip.text.x,
  strip.text.x.bottom,
  strip.text.x.top,
  strip.text.y,
  strip.text.y.left,
  strip.text.y.right,
  strip.switch.pad.grid,
  strip.switch.pad.wrap,
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  axis.text=element_text(size=16),
  legend.text=element_text(size=16),
  panel.background=element_rect(fill="white", colour = "black"),
  plot.background=element_rect(color="white"),
  complete = FALSE,
  validate = TRUE)
{
  theme(...,
        line,
        rect,
        text,
        title,
        aspect.ratio,
        axis.title.x,
        axis.title.x.top,
        axis.title.x.bottom,
        axis.title.y,
        axis.title.y.left,
        axis.title.y.right,
        axis.text.x,
        axis.text.x.top,
        axis.text.x.bottom,
        axis.text.y,
        axis.text.y.left,
        axis.text.y.right,
        axis.text.theta,
        axis.text.r,
        axis.ticks,
        axis.ticks.x,
        axis.ticks.x.top,
        axis.ticks.x.bottom,
        axis.ticks.y,
        axis.ticks.y.left,
        axis.ticks.y.right,
        axis.ticks.theta,
        axis.ticks.r,
        axis.minor.ticks.x.top,
        axis.minor.ticks.x.bottom,
        axis.minor.ticks.y.left,
        axis.minor.ticks.y.right,
        axis.minor.ticks.theta,
        axis.minor.ticks.r,
        axis.ticks.length,
        axis.ticks.length.x,
        axis.ticks.length.x.top,
        axis.ticks.length.x.bottom,
        axis.ticks.length.y,
        axis.ticks.length.y.left,
        axis.ticks.length.y.right,
        axis.ticks.length.theta,
        axis.ticks.length.r,
        axis.minor.ticks.length,
        axis.minor.ticks.length.x,
        axis.minor.ticks.length.x.top,
        axis.minor.ticks.length.x.bottom,
        axis.minor.ticks.length.y,
        axis.minor.ticks.length.y.left,
        axis.minor.ticks.length.y.right,
        axis.minor.ticks.length.theta,
        axis.minor.ticks.length.r,
        axis.line,
        axis.line.x,
        axis.line.x.top,
        axis.line.x.bottom,
        axis.line.y,
        axis.line.y.left,
        axis.line.y.right,
        axis.line.theta,
        axis.line.r,
        legend.background,
        legend.margin,
        legend.spacing,
        legend.spacing.x,
        legend.spacing.y,
        legend.key,
        legend.key.size,
        legend.key.height,
        legend.key.width,
        legend.key.spacing,
        legend.key.spacing.x,
        legend.key.spacing.y,
        legend.frame,
        legend.ticks,
        legend.ticks.length,
        legend.axis.line,
        legend.text.position,
        legend.title,
        legend.title.position,
        legend.position,
        legend.position.inside,
        legend.direction,
        legend.byrow,
        legend.justification,
        legend.justification.top,
        legend.justification.bottom,
        legend.justification.left,
        legend.justification.right,
        legend.justification.inside,
        legend.location,
        legend.box,
        legend.box.just,
        legend.box.margin,
        legend.box.background,
        legend.box.spacing,
        panel.border,
        panel.spacing,
        panel.spacing.x,
        panel.spacing.y,
        panel.grid,
        panel.grid.major.x,
        panel.grid.major.y,
        panel.grid.minor.x,
        panel.grid.minor.y,
        panel.ontop,
        plot.title,
        plot.title.position,
        plot.subtitle,
        plot.caption,
        plot.caption.position,
        plot.tag,
        plot.tag.position,
        plot.tag.location,
        plot.margin,
        strip.background,
        strip.background.x,
        strip.background.y,
        strip.clip,
        strip.placement,
        strip.text,
        strip.text.x,
        strip.text.x.bottom,
        strip.text.x.top,
        strip.text.y,
        strip.text.y.left,
        strip.text.y.right,
        strip.switch.pad.grid,
        strip.switch.pad.wrap,
    plot.background=plot.background,
        axis.text=axis.text,
        axis.title=axis.title,
        legend.text=legend.text,
        panel.grid.major = panel.grid.major,
        panel.grid.minor = panel.grid.minor,
        panel.background=panel.background)
}

theme_la <- function(...){
  theme() %+replace%    #replace elements we want to change
  theme(
    axis.title=element_text(size=20),,
        panel.grid.minor=element_blank(),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        panel.background=element_rect(fill="white", colour = "black"),
        plot.background=element_rect(color="white"))
}

la_replace <-  theme(
  panel.grid.minor=element_blank(),
  panel.grid.major=element_blank(),
    axis.title=element_text(size=20),
    axis.text=element_text(size=16),
    legend.text=element_text(size=16),
    panel.background=element_rect(fill="white", colour = "black"),
    plot.background=element_rect(color="white"))

# shift legend into empty facets
library(gtable)
library(cowplot)

shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}
    