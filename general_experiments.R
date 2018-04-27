#!/usr/bin/Rscript --vanilla
require(methods)
library(data.table)
library(parallel)
library(doMC)
library(doRNG)
library(foreach)
library(glmnet)
library(randomForest)
library(kernlab)
library(rrBLUP)
library(caret)

#------------------------------------------------------------------------------#
# General Functions #
# Generate test indices where n is the number of samples.
GenerateTestIndices <- function(sample.number, percent.split = 0.30){
	sample.size <- floor(percent.split * sample.number)
	set.seed(1357)
	test.indices <- sample(seq_len(sample.number), size = sample.size)

    return(test.indices)
}

# Get varieties for a trait that are not missing.
ReadTraitData <- function(trait){ 
  df <- data.matrix(data.frame(fread(paste0("../trait_data/", trait, ".csv"),
                                     header = TRUE),
                               row.names = 1))
  fcount <- ncol(df)  
  x <- df[, 1: fcount-1] 
  y <- df[, fcount]
  trait.data <- list(x = x, y = y)

  return(trait.data)
}

# Get general peformance metrics. 
GetPerformanceMetrics <- function(y.actual, y.predicted){
  rsquared <- 1 - (sum((y.actual - y.predicted)^2) / 
                   sum((y.actual - mean(y.actual))^2))
  mse <- sum((y.actual - y.predicted)^2) / length(y.actual)
  rmse <- sqrt(sum((y.actual - y.predicted)^2) / length(y.actual))
    
  performance <- list(rsquared = rsquared,
                      mse = mse,
                      rmse = rmse)

  return(performance)
}

#------------------------------------------------------------------------------#
# Learners
FitLasso <- function(x, y, parallel_){
  set.seed(88982)
  lasso.fit <- cv.glmnet(x, y,
                         standardize = FALSE, 
                         alpha = 1,
                         #lambda = 10^seq(10, -3, length=100),
                         parallel = parallel_)
  
  return(lasso.fit)
}

FitRidgeReg <- function(x, y, parallel_){
  set.seed(99843)
  ridge.fit <- cv.glmnet(x, y,
                         standardize = FALSE,
                         alpha = 0,
                         #lambda = 10^seq(10, -3, length=100),
                         parallel = parallel_)
	
  return(ridge.fit)
}

FitBlup <- function(x, y){
  set.seed(35246)
  trait.response <- mixed.solve(y, 
								Z = x,
								K = NULL, 
								SE = FALSE, 
								return.Hinv = FALSE)	
  e <- as.matrix(trait.response$u)
  
  return(trait.response)
}

FitRandomForest <- function(x, y, ntrees=1000, parallel_=FALSE, ncores=NULL){
  if (parallel_){
    ncores = ncores
    set.seed(98364)
    rf.model <- foreach(ntree = rep(floor(ntrees/ncores), ncores),
                        .combine = combine, 
                        .multicombine = TRUE, 
                        .packages = "randomForest") %dorng% {
                  randomForest(x, y, ntree = ntree)  
                }

    return(rf.model)
  } else {
    rf.model <- randomForest(x, y, ntree = ntrees) 
        
    return(rf.model)
  }
}

FitGbm <- function(x, y, parallel_){
  ctrl <- trainControl(method = "cv", number = 5, allowParallel = parallel_)
  gbmgrid <- expand.grid(n.trees = 1000,
                         interaction.depth = c(1, 3, 7),
					     shrinkage = c(0.001, 0.01, 0.1),
					     n.minobsinnode = 10) 
  set.seed(56732)
  gbmfit <- train(x = x,
                  y = y,
				  method = "gbm",
				  distribution = "gaussian",
				  bag.fraction = 0.5,
				  verbose = FALSE,
				  trControl = ctrl,
				  tuneGrid = gbmgrid)

  return(gbmfit)
}

FitSvm <- function(x, y, parallel_){
  ctrl <- trainControl(method = "cv", number = 5, allowParallel = parallel_)
  svmgrid <- expand.grid(sigma = 2^seq(-15, 0, by = 1),
                         C = c(0.1, 1, 10))
  set.seed(53475)
  svm.selection <- train(x = x, 
                         y = y,
                         scaled = FALSE,
                         method = "svmRadial",
                         metric = "RMSE", 
                         trControl = ctrl,
                         tuneGrid = svmgrid) 
  svm.model <- svm.selection$finalModel
    
  return(svm.model)
}

# The last column in the matrix passed is the actual value.
CalculateCvRsquareds <- function(predictions){
  fcount <- ncol(predictions)
  actual <- predictions[, fcount] 
  learner.predictions <- predictions[, 1: fcount-1]
  learners <- colnames(learner.predictions) 
  learners.perf <- NULL 
  for (learner in learners) {
    learner.prediction <- learner.predictions[, learner] 
    learner.perf <- GetPerformanceMetrics(actual, learner.prediction)
    learner.rsquared <- learner.perf$rsquared
    learners.perf <- c(learners.perf, learner.rsquared)
  }
  names(learners.perf) <- learners

  return(learners.perf)
}

#------------------------------------------------------------------------------#
# Experiment Logic
FitMultipleModels <- function(x, y, test.indices, trait, parallel_=TRUE){
  if (parallel_ && !getDoParRegistered()){ 
    registerDoMC(cores = 10) 
  }

  # Train - Test Separation 
  x.train <- x[-test.indices, ]
  y.train <- y[-test.indices]
  x.test <- x[test.indices, ]
  y.test <- y[test.indices]
      
  # Fit Models
  print(paste0(trait, "_", "blup")) 
  blup.model <- FitBlup(x.train, y.train) 
  print(paste0(trait, "_", "lasso")) 
  lasso.model <- FitLasso(x.train, y.train, parallel_ = parallel_) 
  print(paste0(trait, "_", "ridge")) 
  ridge.model <- FitRidgeReg(x.train, y.train, parallel_ = parallel_)  
  print(paste0(trait, "_", "rf")) 
  rf.model <- FitRandomForest(x.train, y.train, parallel_ = parallel_, 
                              ncores = 10)    
  registerDoMC(cores = 5) 
  print(paste0(trait, "_", "svm")) 
  svm.model <- FitSvm(x.train, y.train, parallel_ = parallel_)  
  print(paste0(trait, "_", "gbm"))  
  #registerDoMC(cores = 5) 
  gbmfit <- FitGbm(x.train, y.train, parallel_ = parallel_)     
  best.tune <- unlist(gbmfit$bestTune)
  ntrain <- round(0.7 * nrow(x.train)) 
  gbm.model <- gbm.fit(x = x.train, 
                       y = y.train,
                       distribution = "gaussian",
                       n.trees = best.tune[1],
                       interaction.depth = best.tune[2],
                       shrinkage = best.tune[3],
                       nTrain = ntrain,
                       verbose = FALSE)

  # Log best performing parameters for GBM  
  write(best.tune, 
        append = TRUE, 
        ncolumns = length(best.tune),
        file = paste0("../output/general_experiments/model_parameters/", 
                      trait, "_gbm.txt")
       )

  # Test Models
  e <- as.matrix(blup.model$u)
  trait.pred <- x.test %*% e 
  blup.pred <- (trait.pred[,1]) + blup.model$beta
  blup.perf <- GetPerformanceMetrics(y.test, blup.pred) 

  lasso.pred <- predict(lasso.model, s = lasso.model$lambda.min, 
                        newx = x.test) 
  lasso.perf <- GetPerformanceMetrics(y.test, lasso.pred)

  ridge.pred <- predict(ridge.model, s = ridge.model$lambda.min,
                        newx = x.test) 
  ridge.perf <- GetPerformanceMetrics(y.test, ridge.pred)
  
  rf.pred <- predict(rf.model, x.test)
  rf.perf <- GetPerformanceMetrics(y.test, rf.pred)
   
  best.iter <- gbm.perf(gbm.model, plot.it = FALSE, method = "test") 
  gbm.pred <- predict(gbm.model, x.test, best.iter) 
  gbm.perf <- GetPerformanceMetrics(y.test, gbm.pred) 

  svm.pred <- predict(svm.model, x.test)
  svm.perf <- GetPerformanceMetrics(y.test, svm.pred)
   
  # Compile fold R2
  results <- c(blup = blup.perf$rsquared,
               lasso = lasso.perf$rsquared,
               rr = ridge.perf$rsquared,
               rf = rf.perf$rsquared,
               gbm = gbm.perf$rsquared, 
               svm = svm.perf$rsquared)
  write(results,
        append = TRUE,
        ncolumns = length(results),
        file = "../output/general_experiments/all_rsquareds.txt", 
        )
    
  return(results)
}

#------------------------------------------------------------------------------#
# Load trait data
PerformExperiments <- function(traits){
  all.rsquareds <- NULL 
  for (trait in traits){
    trait.data <- ReadTraitData(trait) 
    x <- trait.data$x
    y <- trait.data$y
    test.indices <- GenerateTestIndices(length(y)) 
    trait.results <- FitMultipleModels(x, y, test.indices, trait)  
    all.rsquareds <- rbind(all.rsquareds, trait.results) 
  }     
  rownames(all.rsquareds) <- traits
  
  write.csv(file = "../output/general_experiments/all_rsquareds.csv", 
            all.rsquareds)
}

traits <- c("CUDI_REPRO", "CULT_REPRO", "CUNO_REPRO", "GRLT", "GRWD", 
            "GRWT100", "HDG_80HEAD", "LIGLT", "LLT", "LWD", "PLT_POST",
            "SDHT")

PerformExperiments(traits)