FitLasso <- function(data, data.test, all.acc) {
  # cml.clus.serum <- ModelTrain(data)
  # CombineSplits(cml.clus.serum, metric = "error rate")
  
  cor.cols <- findCorrelation(cor(data.matrix(data)), cutoff = .9)
  cor.cols
  
  x = model.matrix(Health_State ~ ., data = data)[,-1]
  y = data$Health_State
  grid = 10^seq(10, -2, length = 100)
  
  lasso.mod <- glmnet(x, y, family = "binomial", alpha = 1, lambda = grid)
  plot(lasso.mod)
  
  set.seed(835)
  cv.out = cv.glmnet(x, y, family="binomial", alpha = 1, type.measure="class", nfolds = 82,
                     keep=TRUE)
  bestlam = cv.out$lambda.min
  bestlam
  
  lasso.prob <- cv.out$fit.preval[, cv.out$lambda == bestlam]
  lasso.pred <- lasso.prob > .5
  
  # LOOCV accuracy
  all.acc[1, 1] <- mean(lasso.pred == y)
  
  # specificity
  idx <- y == 0
  all.acc[1, 2] <- mean(lasso.pred[idx] == y[idx])
  
  # sensitivity
  idx <- y == 1
  all.acc[1, 3] <- mean(lasso.pred[idx] == y[idx])
  
  # AUC
  all.acc[1, 4] <- as.numeric(auc(y, lasso.prob))
  
  plot(cv.out)
  
  cv.coef <- coef(cv.out, s = bestlam)
  cv.coef <- as.matrix(Matrix(cv.coef, sparse = F))
  cv.coef[cv.coef != 0, ]
  
  # prediction performance on test set
  
  #remove correlated columns
  
  # data.test <- data.test[, -cor.cols]
  x.test = model.matrix(Health_State ~ ., data = data.test)[, -1]
  y.test = data.test$Health_State
  lasso.pred <- predict(lasso.mod, s = bestlam, newx = x.test, type = "class")
  lasso.pred <- as.numeric(lasso.pred[, 1])
  
  # test accuracy
  all.acc[1, 5] <- mean(lasso.pred == y.test)
  
  # specificity.test
  idx <- y.test == 0
  all.acc[1, 6] <- mean(y.test[idx] == lasso.pred[idx])
  
  # sensitivity.test
  idx <- y.test == 1
  all.acc[1, 7] <- mean(y.test[idx] == lasso.pred[idx])
  
  # AUC
  lasso.prob <- predict(lasso.mod, s = bestlam, newx = x.test, type = "response")
  all.acc[1, 8] <- as.numeric(auc(y.test, lasso.prob[, 1]))
  
  return(all.acc)
}


FitMlmodels <- function(data, data.test, all.acc) {
  
  customSummary <- function (data, lev = NULL, model = NULL){
    spec <- specificity(data[, "pred"], data[, "obs"], lev[2])
    pred <- factor(ifelse(data[, "neg"] > 0.5, "neg", "pos"))
    spec2 <- specificity(pred, data[, "obs"], "pos")
    acclist <- as.numeric()
    for(i in seq(0.01, 0.99, 0.01)){
      predi <- factor(ifelse(data[, "neg"] > i, "neg", "pos"), levels = c("neg", "pos"))
      singleacc <- mean(predi == data[, "obs"])
      acclist <- c(acclist, singleacc)
    }
    max(acclist) -> accmax
    
    out <- c(spec, spec2, accmax)
    
    names(out) <- c("Spec", "Spec2", "Accmax")
    out
  }
  
  fitControl <- trainControl(method = "LOOCV", classProbs = TRUE, savePredictions = TRUE,
                             summaryFunction = customSummary)
  models <- c("LogitBoost", 'svmRadial', 'rf', 'pls', 'xgbTree')
  
  for (j in 1:5) {
    set.seed(9)
    if (j != 2) {
      model.fit <- train(Health_State~., method=models[j], data=data, metric = "Accmax",
                         trControl = fitControl, tuneLength = 10)
    } else {
      svmGrid <- expand.grid(sigma= 2^c(-25, -20, -15,-10, -5, 0), C= 2^c(0:5))
      model.fit <- train(Health_State~., method=models[j], data=data, metric = "Accmax",
                         trControl = fitControl, tuneGrid = svmGrid)
      sigma <- log(model.fit$bestTune$sigma, base = 2)
      C <- log(model.fit$bestTune$C, base = 2)
      set.seed(9)
      svmGrid <- expand.grid(sigma= 2^seq(sigma - .5, sigma + .5, length.out = 10),
                             C= 2^seq(C - 1, C + 1, length.out = 10))
      model.fit <- train(Health_State~., method=models[j], data=data, metric = "Accmax",
                         trControl = fitControl, tuneGrid = svmGrid)
    }
    
    maxcv <- max(model.fit$results$Accmax)
    
    # Need to find the threshold
    mod.pred <- as.data.frame(model.fit$pred[, 6:ncol(model.fit$pred)])
    tune.par <- model.fit$bestTune
    fix.order <- order(names(model.fit$bestTune))
    tune.par <- tune.par[fix.order]
    fix.order <- order(colnames(mod.pred))
    mod.pred <- as.data.frame(mod.pred[, fix.order])
    
    sel <- vector(length = nrow(mod.pred))
    for (i in 1:nrow(mod.pred)) {
      sel[i] <- all(mod.pred[i, ] == tune.par)
    }
    
    acclist <- as.numeric()
    for(i in seq(0.01, 0.99, 0.001)){
      predi <- factor(ifelse(model.fit$pred[sel, "pos"] >= i, "pos", "neg"), levels = c("pos", "neg"))
      singleacc <- mean(predi == model.fit$pred$obs[sel])
      acclist <- c(acclist, singleacc)
    }
    maxcv <- max(acclist)
    thresh <- seq(0.01, 0.99, 0.001)[which.max(acclist)]
    
    y <- data$Health_State
    
    probs.train <- model.fit$pred[sel, "pos"]
    pred.train <- ifelse(model.fit$pred[sel, "pos"] >= thresh, "pos", "neg")
    
    # LOOCV accuracy
    all.acc[j+1, 1] <- mean(pred.train == y)
    
    # sensitivity
    idx <- y == "pos"
    all.acc[j+1, 2] <- mean(pred.train[idx] == y[idx])
    
    # specificity
    idx <- y == "neg"
    all.acc[j+1, 3] <- mean(pred.train[idx] == y[idx])
    
    # AUC
    y.num <- ifelse(y == "pos", 1, 0)
    all.acc[j+1, 4] <- as.numeric(auc(y.num, probs.train))
    
    # External
    y.test <- data.test$Health_State
    probs.test <- predict.train(model.fit, data.test, type = "prob")
    pred.test <- ifelse(probs.test$pos > thresh, "pos", "neg")
    
    # test accuracy
    all.acc[j+1, 5] <- mean(pred.test == y.test)
    
    # sensitivity.test
    idx <- y.test == "pos"
    all.acc[j+1, 6] <- mean(y.test[idx] == pred.test[idx])
    
    # specificity.test
    idx <- y.test == "neg"
    all.acc[j+1, 7] <- mean(y.test[idx] == pred.test[idx])
    
    # AUC
    y.num <- ifelse(y.test == "pos", 1, 0)
    all.acc[j+1, 8] <- as.numeric(auc(y.num, probs.test[, "pos"]))
  } 
  
  return(all.acc)
}

