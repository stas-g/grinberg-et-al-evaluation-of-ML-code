#---------------------------------------------------------------------------------------------------------------------------------------#
#                          CODE FOR ANALYSIS ON AN ALTERNATIVE DATA SET: geno.fuse (6064 BINARY ATTRIBUTES)                       #
#---------------------------------------------------------------------------------------------------------------------------------------#

geno.fuse <- readRDS('geno_fuse.rds')
load("yeast_pheno.RData")
source("functions.r")
source('xpack.r')

library(doRNG)
library(doMC)
registerDoMC(cores = 46)


#=========================================================================================================================
#RANDOM FOREST
library(randomForest)

write(paste('RF full, CV, start:', Sys.time()), file = log_file, append = TRUE)

rf.cvpred <- foreach(i = 1 : ncol(pheno2), .options.RNG = 100, .combine = cbind) %dorng% {
	pred <- forest.xpred(geno.fuse, pheno2[, i], yfolds, err = FALSE, importance = TRUE, ntree = 700, do.trace = 100)
	write(paste('RF full,', colnames(pheno2)[i], ', finish:', Sys.time()), file = log_file, append = TRUE)
	pred
		}
dimnames(rf.cvpred) <- dimnames(pheno2)


#=========================================================================================================================
#GBM
library(gbm)

gbm.cvpred <- foreach(i = 1 : ncol(pheno2), .options.RNG = 100, .combine = cbind) %dorng% {
	nTrain_ <- sapply(yfolds, FUN = function(ind) sum(!is.na(pheno2[-ind, i])))
	nTrain_ <- floor(3/4 * mean(nTrain_))
	pred <- gbm.xpred(geno.fuse, pheno2[, i], yfolds, err = FALSE, distribution = 'gaussian', n.trees = 1000, interaction.depth = 5, n.minobsinnode = 10, keep.data = FALSE, verbose = TRUE, shrinkage = 0.01, bag.fraction = 0.5, nTrain = nTrain_)
	pred
	}
dimnames(gbm.cvpred) <- dimnames(pheno2)


#=========================================================================================================================
#LASSO
library(glmnet)

lasso.cvpred <- foreach(i = 1 : ncol(pheno2), .options.RNG = 100, .combine = cbind) %dorng% {
				pred <- lasso.xpred(geno.fuse, pheno2[, i], yfolds, err = FALSE, family = 'gaussian', alpha = 1, standardize = FALSE, nfolds = 10, parallel = F)
				saveRDS(pred, file = file_)
				pred
			}
dimnames(elnet.cvpred) <- dimnames(pheno2)


#=========================================================================================================================
#RIDGE
library(glmnet)

ridge.cvpred <- foreach(i = 1 : ncol(pheno2), .options.RNG = 100, .combine = cbind) %dorng% {
				pred <- ridge.xpred(geno.fuse, pheno2[, i], yfolds, err = FALSE, family = 'gaussian', alpha = 0, standardize = FALSE, nfolds = 10, parallel = F)
				pred
			}
dimnames(ridge.cvpred) <- dimnames(pheno2)


#=========================================================================================================================

rf_rsq <- calc.mse(obs = pheno2, pred = rf.cvpred, rsq = TRUE)
gbm_rsq <- calc.mse(obs = pheno2, pred = gbm.cvpred, rsq = TRUE)
lasso_rsq <- calc.mse(obs = pheno2, pred = lasso.cvpred, rsq = TRUE)
ridge_rsq <- calc.mse(obs = pheno2, pred = ridge.cvpred, rsq = TRUE)

fuse_rsq <- cbind(rf = rf_rsq, gbm = gbm_rsq, lasso = lasso_rsq, ridge = ridge_rsq) 

saveRDS(fuse_rsq, file = 'fuse_rsq.rds')
	
	
	
	
	
	





