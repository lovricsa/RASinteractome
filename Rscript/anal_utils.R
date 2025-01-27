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
                           data = ras_data, family = "binomial")
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


