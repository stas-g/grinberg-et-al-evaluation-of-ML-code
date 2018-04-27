#!/usr/bin/Rscript --vanilla
require(methods)
library(data.table)
library(doMC)
library(parallel)
library(glmnet)

#-----------------------------------------------------------------------------#
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

#-----------------------------------------------------------------------------#
# Learner and experiment logic
FitEnet <- function(x, y, alpha, parallel_){
  set.seed(83982)
  enet.fit <- cv.glmnet(x, y, 
                        standardize = FALSE, 
                        alpha = alpha,
                        #lambda = 10^seq(10, -3, length=100),
                        parallel = parallel_)
  
  return(enet.fit)
}

PerformTT <- function(x, y, alphas, test.indices, trait, parallel_){
  # All variables prefixed with "mean" hold values used in the standard
  # calculation of rsquared.
  x.train <- x[-test.indices, ]
  y.train <- y[-test.indices]
  x.test <- x[test.indices, ]
  y.test <- y[test.indices]
   
  alphas.rsquareds <- NULL
  for (alpha in alphas) {    
    enet.model <- FitEnet(x.train, y.train, alpha, parallel_)
    enet.pred <- predict(enet.model, 
                           s = enet.model$lambda.min, 
                           newx = x.test)  
    alpha.rsquared <- GetPerformanceMetrics(y.test, enet.pred[, 1])$rsquared  
    alphas.rsquareds <- c(alphas.rsquareds, alpha.rsquared) 
  }
  
  names(alphas.rsquareds) <- alphas
      
  return(alphas.rsquareds)
}

#-----------------------------------------------------------------------------#
# Load trait data
PerformExperiments <- function(traits, parallel_ = TRUE){
  if (parallel_ && !getDoParRegistered()){ 
    registerDoMC(cores = detectCores()) 
  }  
    
  alphas <- seq(0, 1, by = 0.1) 
  all.rsquareds <- NULL  
  for (trait in traits){
    trait.data <- ReadTraitData(trait) 
    x <- trait.data$x
    y <- trait.data$y 
    test.indices <- GenerateTestIndices(length(y))  
    rsquareds <- PerformTT(x, y, alphas, test.indices, trait, parallel_)   
    all.rsquareds <- rbind(all.rsquareds, rsquareds) 
  }      
  rownames(all.rsquareds) <- traits
  colnames(all.rsquareds) <- alphas

  write.csv(file = "../output/enet_experiments/all_rsquareds.csv", 
            all.rsquareds)
}

traits <- c("CUDI_REPRO", "CULT_REPRO", "CUNO_REPRO", "GRLT", "GRWD", 
            "GRWT100", "HDG_80HEAD", "LIGLT", "LLT", "LWD", "PLT_POST",
            "SDHT")

PerformExperiments(traits)