#!/usr/bin/Rscript --vanilla
# setwd("C:/Users/nastasya/Documents/manchester/wheat")
# library(rrBLUP)
library(randomForest)
library(gbm)
# library(glmnet)
library(e1071)
library(caret)
library(magrittr)

load('ready_data.r')
source('xpack.r')
source('yeast_fun.r')

library(doRNG)
library(doMC)
registerDoMC(cores = 30)

log_file <- 'wheat.log'
file.create(log_file)
write(paste('script started:', Sys.time()), file = log_file, append = TRUE)
write(" ", append = TRUE)

# log_file <- 'model_output/cvscript.log'

#-----------------------------------------------
#PREP

# wpheno_raw <- read.csv("C:/Users/nastasya/Documents/manchester/wheat/pheno.csv", sep = ",", stringsAsFactors = FALSE)

# wgeno_raw <- read.csv("C:/Users/nastasya/Documents/manchester/wheat/tpg12-06-0006-dataset-s2/DatasetS2_29SAWSN_gbs.csv", sep = ",", stringsAsFactors = FALSE)
# wgeno <- wgeno_raw[!duplicated(wgeno_raw$rs), ]
# wgeno <- wgeno[, -c(1 : 4)]
# wgeno <- t(wgeno)

# wgeno[wgeno == 'N'] <- 'H'
# wgeno <- apply(wgeno, 2, FUN = function(x){
				# tab <- table(x)[names(table(x)) != 'H']
				# # x[x == names(which.max(tab))] <- 1
				# # x[x == names(which.min(tab))] <- -1
				# x[x == names((tab))[1]] <- 1
				# x[x == names((tab))[2]] <- -1
				# x
	# })
# wgeno[wgeno == "H"] <- 0
# class(wgeno) <- 'numeric'

# wpheno_raw <- wpheno_raw[, -c(1 : 4)]
# wpheno <- scale(wpheno_raw)

# set.seed(100)
# wfolds <- cv.folds(nrow(wpheno), folds = 5)

# save(list = c('wpheno', 'wpheno_raw', 'wgeno', 'wfolds'), file = 'ready_data.r')
#---------------------------------------------------------------------------------------

# #BLUP
# write(paste('blup start:', Sys.time()), file = log_file, append = TRUE)

# blup_mcmc_pred <- foreach(i = 1 : 10, .options.RNG = 100) %dorng% {
				# set.seed(i)
				# wfolds <- cv.folds(nrow(wpheno), folds = 5)
				
				# blup_pred <- apply(wpheno, 2, FUN = function(x) blup.cv(x, wgeno, folds = wfolds))
				# print(i)
				# write(paste('blup MCMC run', i, 'done'), file = log_file, append = TRUE)
				# blup_pred
# }


# saveRDS(blup_mcmc_pred, file = 'model_output/blup/blup_mcmc.rds')

# blup_rsq <- sapply(blup_mcmc_pred, FUN = function(x) calc.mse(wpheno, x, rsq = TRUE))
# saveRDS(blup_rsq, file = 'model_output/blup/blup_mcmc_rsq.rds')

# write(paste('blup MCMC done:', Sys.time()), file = log_file, append = TRUE)
# write(' ', file = log_file, append = TRUE)

# # #----------------
# # #GLMNET
# write(paste('ELNET start:', Sys.time()), file = log_file, append = TRUE)

# glmnet_mcmc_pred <- lapply(seq(0, 1, by = 0.1), FUN = function(a){
	# glmnet_pred <- foreach(i = 1 : 10, .options.RNG = 100) %dorng% {
					# set.seed(i)
					# wfolds <- cv.folds(nrow(wpheno), folds = 5)
					
					# pred <- apply(wpheno, 2, FUN = function(x) lasso.xpred(wgeno, x, folds = wfolds, err = FALSE, family = 'gaussian', alpha = a, standardize = FALSE, parallel = FALSE))
					# print(i)
					# write(paste('elnet a =', a, 'MCMC run', i, 'done'), file = log_file, append = TRUE)
					# pred
	# }
		# print(a)
	# glmnet_pred
# })

# saveRDS(glmnet_mcmc_pred, file = 'model_output/elnet/glmnet_mcmc_pred.rds')

# glmnet_mcmc_rsq <- lapply(glmnet_mcmc, FUN = function(z) sapply(z, FUN = function(x) calc.mse(wpheno, x, rsq = TRUE)))
# saveRDS(glmnet_mama_rsq, file = 'model_output/elnet/glmnet_mcmc_rsq.rds')

# write(paste('elnet MCMC done:', Sys.time()), file = log_file, append = TRUE)
# write('', file = log_file, append = TRUE)

# #----------------
# #RF (with mtry tuning)

# write(paste('rf start:', Sys.time()), file = log_file, append = TRUE)

# rf_mcmc_pred <- foreach(i = 1 : 10, .options.RNG = 100) %dorng% {
					# set.seed(i)
					# wfolds <- cv.folds(nrow(wpheno), folds = 5)
					
					# rf_pred <- apply(wpheno, 2, FUN = function(x) forest.xpred.tune(wgeno, x, folds = wfolds, err = FALSE, importance = TRUE, ntree = 700, do.trace = 100))
					# write(paste('RF MCMC run', i, 'done'), file = log_file, append = TRUE)
					# rf_pred
# }

# saveRDS(rf_mcmc_pred, file = 'model_output/rf/rf_mcmc_pred.rds')

# rf_rsq <- sapply(rf_mcmc_pred, FUN = function(x) calc.mse(wpheno, x, rsq = TRUE))

# saveRDS(rf_rsq, file = 'model_output/rf/rf_mcmc_rsq.rds')

# write(paste('RF MCMC done:', Sys.time()), file = log_file, append = TRUE)
# write(' ', file = log_file, append = TRUE)

#----------------
#GBM 

gbm_mcmc_out <- lapply(1 : 10, FUN = function(k) {
	set.seed(k)
	wfolds <- cv.folds(nrow(wpheno), folds = 5)

	gbm_out <- lapply(1 : ncol(wpheno), FUN = function(m){
		y <- wpheno[, m]
		nTrain_ <- floor(3/4 * (nrow(wpheno) - length(wfolds[[1]])))
		out <- gbm.xpred.tune(wgeno, y, wfolds, nTrain = nTrain_, parallel_ = TRUE, file_ = sprintf('model_output/gbm/cv_tmp/run%s_%s_', k, colnames(wpheno)[m]))
		write(paste('GBM,', colnames(wpheno)[m], ', finish:', Sys.time()), file = log_file, append = TRUE)
		message(sprintf("run%s: trait %s", k, colnames(wpheno)[m]))
		out
		})
	gbm_pred <- lapply(gbm_out, "[[", 1)
	gbm_bestTune <- lapply(gbm_out, "[[", 2)
	
	list(pred = gbm_pred, bestTune = gbm_bestTune)
	}	
)

saveRDS(gbm_mcmc_out, file = 'model_output/gbm/gbm_mcmc_out.rds')
gbm_bestTune <- lapply(gbm_mcmc_out, "[[", 2)
gbm_pred <- lapply(gbm_mcmc_out, FUN = function(x) do.call(cbind, x[[1]]))

gbm_rsq <- sapply(gbm_pred, FUN = function(x) calc.mse(wpheno, x, rsq = TRUE))

saveRDS(gbm_rsq, file = 'model_output/gbm/gbm_mcmc_rsq.rds')
saveRDS(gbm_bestTune, file = 'model_output/gbm/gbm_bestTune.rds')
saveRDS(gbm_pred, file = 'model_output/gbm/gbm_mcmc_pred.rds')

#----------------
#SVM (with tuning)

write(paste('svm start:', Sys.time()), file = log_file, append = TRUE)

#3

svm_mcmc_pred <- foreach(k = 1 : 5, .options.RNG = 100) %dorng% {
					set.seed(k)
					wfolds <- cv.folds(nrow(wpheno), folds = 5)
				
				svm_pred <- lapply(1 : ncol(wpheno), FUN = function(m){
					y <- wpheno[, m]
					ans <- svm.xpred.tune(wgeno, y, wfolds, nTrain = nTrain_, parallel_ = TRUE, file_ = sprintf('model_output/svm/cv_tmp/run%s_%s_', k, colnames(wpheno)[m]))
					write(paste0(colnames(wpheno)[m], 'done'), file = log_file, append = TRUE)
					ans
						})
				svm_pred <- do.call(cbind, svm_pred)
				saveRDS(svm_pred, file = sprintf('model_output/svm/svm_mcmc_pred_run%s.rds', k))
				svm_pred
}

saveRDS(svm_mcmc_pred, file = 'model_output/svm/svm_mcmc_pred.rds')

svm_rsq <- sapply(svm_mcmc_pred, FUN = function(x) calc.mse(wpheno, x, rsq = TRUE))
saveRDS(svm_rsq, file = 'model_output/svm/svm_mcmc_rsq.rds')


write(paste('SVM MCMC done:', Sys.time()), file = log_file, append = TRUE)
write(' ', file = log_file, append = TRUE)

write('done!', file = log_file, append = TRUE)


