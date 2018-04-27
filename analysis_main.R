#---------------------------------------------------------------------------------------------------------------------------------------#
#                                                 CODE FOR MAIN ANALYSIS (RESULTS IN TABLE 1)                                           #
#---------------------------------------------------------------------------------------------------------------------------------------#

geno <- readRDS('geno.rds')
load("yeast_pheno.RData")
source("functions.r")
source('xpack.r')


library(doRNG)
library(doMC)
registerDoMC(cores = 46)


#--------------------------------------------------------------------------------------
#random froest
library(randomForest)

#CV
rf.cvpred <- foreach(i = 1 : ncol(pheno2), .options.RNG = 100, .combine = cbind) %dorng% {
	pred <- forest.xpred(geno, pheno2[, i], yfolds, err = FALSE, importance = TRUE, ntree = 700, do.trace = 100)
	pred
		}
dimnames(rf.cvpred) <- dimnames(pheno2)

#---------------------------------
#test/train 

forest.test <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
	ytrain <- traits.m2[[i]][-test.ind[[i]]]
	ytest <- traits.m2[[i]][test.ind[[i]]]
	mod <- randomForest(x = geno[names(ytrain), ], y = ytrain, importance = TRUE, xtest = geno[names(ytest), ], ytest = ytest, ntree = 700, do.trace = 100)
	pred <- mod$test$predicted
	return(list(mod = mod, pred = pred))
		}  

forest.test.mod <- lapply(forest.test, '[[', 1)
names(forest.test.mod) <- names(traits.m)

forest.test.pred <- lapply(forest.test, '[[', 2)
names(forest.test.pred) <- names(traits.m)

rf.test.rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], forest.test.pred[[i]], rsq = TRUE))

#=========================================================================================================================
#GBM
library(gbm)

#CV
gbm.cvpred <- foreach(i = 1 : ncol(pheno2), .options.RNG = 100, .combine = cbind) %dorng% {
	nTrain_ <- sapply(yfolds, FUN = function(ind) sum(!is.na(pheno2[-ind, i])))
	nTrain_ <- floor(3/4 * mean(nTrain_))
	pred <- gbm.xpred(geno, pheno2[, i], yfolds, err = FALSE, distribution = 'gaussian', n.trees = 1000, interaction.depth = 5, n.minobsinnode = 10, keep.data = FALSE, verbose = TRUE, shrinkage = 0.01, bag.fraction = 0.5, nTrain = nTrain_)
	pred
	}
dimnames(gbm.cvpred) <- dimnames(pheno2)

#---------------------------------
#train/test

gbm.test.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
	ytrain <-traits.m2[[i]][-test.ind[[i]]]
	nTrain_ <- 3/4 * length(ytrain)
	mod <- gbm.fit(geno[names(ytrain), ], ytrain, distribution = 'gaussian', n.trees = 1000, interaction.depth = 5, n.minobsinnode = 10, shrinkage = 0.01, bag.fraction = 0.5, nTrain = nTrain_, keep.data = FALSE, verbose = TRUE)
	mod
	}	
names(gbm.test.mod) <- names(traits.m)

gbm.test.pred <- sapply(1 : 46, FUN = function(i) predict(gbm.test.mod[[i]], geno[names(traits.m2[[i]][test.ind[[i]]]), ]))
names(gbm.test.pred) <- names(traits.m)
gbm.test.rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], gbm.test.pred[[i]], rsq = TRUE))


#--------------------------------------------------------------------------------------
#LASSO
library(glmnet)

#CV
lasso.cvpred <- foreach(i = 1 : ncol(pheno2), .options.RNG = 100, .combine = cbind) %dorng% {
				pred <- lasso.xpred(geno, pheno2[, i], yfolds, err = FALSE, family = 'gaussian', alpha = 1, standardize = FALSE, nfolds = 10, parallel = F)
				pred
			}

#---------------------------------
#train/test

lasso.test.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
		y <-traits.m2[[i]][-test.ind[[i]]]
		mod <- cv.glmnet(geno[names(y), ], y, family = 'gaussian', alpha = 1, standardize = FALSE, nfolds = 10, parallel = F)
		mod
}

names(lasso.test.mod) <- names(traits.m2)

lasso.test.pred <- sapply(1 : 46, FUN = function(i) predict(lasso.test.mod[[i]], geno[names(traits.m2[[i]][test.ind[[i]]]), ], s = 'lambda.1se'))
sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], lasso.test.pred[[i]], rsq = TRUE))
lasso_test_rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], lasso.test.pred[[i]], rsq = TRUE))

#--------------------------------------------------------------------------------------
RIDGE
library(glmnet)

#CV
ridge.cvpred <- foreach(i = 1 : ncol(pheno2), .options.RNG = 100, .combine = cbind) %dorng% {
				pred <- lasso.xpred(geno, pheno2[, i], yfolds, err = FALSE, family = 'gaussian', alpha = 0, standardize = FALSE, nfolds = 10, parallel = FALSE)
				pred
			}

#---------------------------------
#train/test

ridge.test.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
		y <-traits.m2[[i]][-test.ind[[i]]]
		mod <- cv.glmnet(geno[names(y), ], y, family = 'gaussian', alpha = 0, standardize = FALSE, nfolds = 10, parallel = FALSE)
		mod
}

names(ridge.test.mod) <- names(traits.m2)

ridge.test.pred <- sapply(1 : 46, FUN = function(i) predict(ridge.test.mod[[i]], geno[names(traits.m2[[i]][test.ind[[i]]]), ], s = 'lambda.1se'))
names(ridge.test.pred) <- names(traits.m)
ridge_test_rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], ridge.test.pred[[i]], rsq = TRUE))



#------------------------------------------------------------------------------------------------------------
#first element of xx.res is cvR^2, second - mean R^2 accross folds, last - corresponding sd of R^2 across folds. the first and the last figures for each traits are those found in table 1

rf.res <- fold.mse(obs = pheno2, pred = rf.cvpred, folds = yfolds, rsq = TRUE)
lasso.res <- fold.mse(obs = pheno2, pred = lasso.cvpred, folds = yfolds, rsq = TRUE)
ridge.res <- fold.mse(obs = pheno2, pred = ridge.cvpred, folds = yfolds, rsq = TRUE)
gbm.res <- fold.mse(obs = pheno2, pred = gbm.cvpred5, folds = yfolds, rsq = TRUE)

full_rsq <- data.frame(rf = rf.res[[1]], lasso = lasso.res[[1]], ridge = ridge.res[[1]], gbm = gbm.res[[1]])
full_sd <- data.frame(rf = rf.res[[3]], lasso = lasso.res[[3]], ridge = ridge.res[[3]], gbm = gbm.res[[3]])

test_rsq <- data.frame(lasso = lasso_test_rsq, ridge = ridge_test_rsq, gbm = gbm_test_rsq, rf = rf_test_rsq)

#saving R^2 results for CV scheme (main analysis) and for train-test scheme (for data degradation experiments later)
saveRDS(full_rsq, file = 'full_cv_rsq.rds')
saveRDS(full_sd, file = 'full_cv_sd.rds')
saveRDS(test_rsq, file = 'full_test_rsq.rds')
#------------------------------------------------------------------------------------------------------------
















































	
	
	
	
	





