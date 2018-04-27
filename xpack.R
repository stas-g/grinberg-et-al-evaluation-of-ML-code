# collection of functions for cross-validation for several CRAN packages: glmnet, gbm and randomForest
# all take arguments: x (input); y (output); folds for cross validation (10 folds are randomly generated with seed = 100 if folds, n and seed are not supplied); if err = FALSE cross-validated predictions are output, else MSE is returned; if parallel_ = TRUE calculations are parallelised (you will need to install and load doRNG and doMC packages and register an appropriate number X of cores via registerDoMC(cores = X). '...' are for arguments specific for each underlying fit function (e.g. n.trees for gbm).
#missing values in the output are allowed, only non-missing values make into folds and subsequently get predicted (when in test fold)


#GLMNET

lasso.xpred <- function(x, y, folds = NULL, n = 10, err = FALSE, parallel_ = FALSE, seed = NULL, ...){
	if(is.null(folds)) folds <- split(sample(1 : length(y)), rep(1 : n, length = length(y)))
	train.ind <- lapply(folds, FUN = function(ind) (1 : length(y))[-ind][!is.na(y[-ind])])
	test.ind <- lapply(folds, FUN = function(ind) ind[!is.na(y[ind])])
	
	if(parallel_ == FALSE){
		pred_tmp <- lapply(1 : length(folds), FUN = function(i){
					fit <- cv.glmnet(x[train.ind[[i]], ], y[train.ind[[i]]], ...)
					predict(fit, x[test.ind[[i]], ], s = "lambda.1se")
				})
					} else {
						if(is.null(seed)) seed <- 100
							pred <- foreach(i = 1 : length(folds), .options.RNG = seed)  %dorng% {
								fit <- cv.glmnet(x[train.ind[[i]], ], y[train.ind[[i]]], ...)
								predict(fit, x[test.ind[[i]], ], s = "lambda.1se")
							}
					}
	pred_tmp <- unlist(pred_tmp)[order(unlist(test.ind))]
	pred <- rep(NA, length(y))
	pred[!is.na(y)] <- pred_tmp
		if(err == FALSE){return(pred)
				}else{	
					mean((pred - y)^2, na.rm = TRUE)}
}



#GBM

gbm.xpred <- function(x, y, folds = NULL, n = 10, err = FALSE, parallel_ = FALSE, seed = NULL, ...){
	if(is.null(folds)) folds <- split(sample(1 : length(y)), rep(1 : n, length = length(y)))
	train.ind <- lapply(folds, FUN = function(ind) (1 : length(y))[-ind][!is.na(y[-ind])])
	test.ind <- lapply(folds, FUN = function(ind) ind[!is.na(y[ind])])
	
	if(parallel_ == FALSE){
		pred_tmp <- lapply(1 : length(folds), FUN = function(i){
					fit <- gbm.fit(x[train.ind[[i]], ], y[train.ind[[i]]], ...)
					predict(fit, x[test.ind[[i]], ])
						})
							}else{
								if(is.null(seed)) seed <- 100
								pred_tmp <- foreach(i = 1 : length(folds), .options.RNG = seed) %dorng% {
									fit <- gbm.fit(x[train.ind[[i]], ], y[train.ind[[i]]], ...)
									predict(fit, x[test.ind[[i]], ])
								}
							}			
	pred_tmp <- unlist(pred_tmp)[order(unlist(test.ind))]
	pred <- rep(NA, length(y))
	pred[!is.na(y)] <- pred_tmp
		if(err == FALSE){return(pred)
				}else{	
					mean((pred - y)^2, na.rm = TRUE)}
}


#RANDOMFOREST

forest.xpred <- function(x, y, folds = NULL, n = 10, err = FALSE, parallel_ = FALSE, seed = NULL, ...){
	if(is.null(folds)) folds <- split(sample(1 : length(y)), rep(1 : n, length = length(y)))
	train.ind <- lapply(folds, FUN = function(ind) (1 : length(y))[-ind][!is.na(y[-ind])])
	test.ind <- lapply(folds, FUN = function(ind) ind[!is.na(y[ind])])
	
	if(parallel_ == FALSE){
		pred_tmp <- lapply(1 : length(folds), FUN = function(i){
						fit <- randomForest(x[train.ind[[i]], ], y[train.ind[[i]]], xtest = x[test.ind[[i]], ], ytest = y[test.ind[[i]]], ...)
						fit$test$predicted
			})
				}else{
						if(is.null(seed)) seed <- 100
						pred_tmp <- foreach(i = 1 : length(folds), .options.RNG = seed) %dorng% {
								fit <- randomForest(x[train.ind[[i]], ], y[train.ind[[i]]], xtest = x[test.ind[[i]], ], ytest = y[test.ind[[i]]], ...)
								fit$test$predicted
				}
					}
					
	pred_tmp <- unlist(pred_tmp)[order(unlist(test.ind))]
	pred <- rep(NA, length(y))
	pred[!is.na(y)] <- pred_tmp
		if(err == FALSE){return(pred)
				}else{	
					mean((pred - y)^2, na.rm = TRUE)}
}

