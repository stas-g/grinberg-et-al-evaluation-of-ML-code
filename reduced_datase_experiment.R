#===============================================================================================#
#       MCMC procedure for the data size experiment (deleting some % of sample data)			#
#===============================================================================================#

geno <- readRDS('geno')
load("data/yeast_pheno.RData")
source("functions.r")
source('xpack.r')

library(doRNG)
library(doMC)
registerDoMC(cores = 23)

library(gbm)
library(randomForest)
library(glmnet)

#===============================================================================================
#vector of proportion of data to be used (i.e. if use 0.15 of the data, then 0.85 is removed)
avec <- c(0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9)
anames <- c('15', '20', '30', '40', '50', '75', '90')

#matrices to be filled with results in the for loops
gbm.rsq <- rf.rsq <- lasso.rsq <- ridge.rsq <- matrix(0, nrow = 46, ncol = 7, dimnames = list(names(traits.m2), anames))
gbm.mcmc.rsq <- rf.mcmc.rsq <- lasso.mcmc.rsq <- ridge.mcmc.rsq <- array(0, c(46, 7, 10), dimnames = list(names(traits.m2), anames, 1 : 10))


#===============================================================================================
#GBM 

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 7){
		a <- avec[k]
		aname <- anames[k]
		
		#crating new phenotypic output by randomly selecting a*100% of data; new random seed is set at teach MCMC iteration
		traits.m2_n <- lapply(traits.m, FUN = function(x){
				n <- length(x)
				set.seed(u * 100)
				ind <- sort(sample(1 : n, round(n * a))) 
				x <- x[ind]
				out <- as.numeric(scale(x))
				names(out) <- names(x)
				out
			})
		
		#creating new test-train folds (random 30% - 70% split for each trait)
		test.ind <- lapply(1 : 46, FUN = function(i){
			set.seed(100)
			n <- length(traits.m2_n[[i]])
			sort(sample(1 : n, ceiling(n * 0.3), replace = FALSE))
				})

		#training model on training test
		gbm.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			ytrain <- traits.m2_n[[i]][-test.ind[[i]]]
			nTrain_ <- 3/4 * length(ytrain)
			mod <- gbm.fit(geno[names(ytrain), ], ytrain, distribution = 'gaussian', n.trees = 1000, interaction.depth = 5, n.minobsinnode = 10, shrinkage = 0.01, bag.fraction = 0.5, nTrain = nTrain_, keep.data = FALSE, verbose = TRUE)
			mod
		}

		#predictions on prediction portion of data
		gbm.pred <- sapply(1 : 46, FUN = function(i) predict(gbm.mod[[i]], geno[names(traits.m2_n[[i]][test.ind[[i]]]), ]))
		names(gbm.pred) <- names(traits.m2_n)
		
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2_n[[i]][test.ind[[i]]], gbm.pred[[i]], rsq = TRUE))
		gbm.rsq[, k] <- rsq
	}

	gbm.mcmc.rsq[, , u] <- gbm.rsq
	gbm.rsq <- matrix(0, nrow = 46, ncol = 7, dimnames = list(names(traits.m2), anames))
}


#===============================================================================================
#RF

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 7){
		a <- avec[k]
		aname <- anames[k]

		#crating new phenotypic output by randomly selecting a*100% of data; new random seed is set at teach MCMC iteration
		traits.m2_n <- lapply(traits.m, FUN = function(x){
					n <- length(x)
					set.seed(u * 100)
					ind <- sort(sample(1 : n, round(n * a))) 
					x <- x[ind]
					out <- as.numeric(scale(x))
					names(out) <- names(x)
					out
				})
		
		#creating new test-train folds (random 30% - 70% split for each trait)		
		test.ind <- lapply(1 : 46, FUN = function(i){
			set.seed(100)
			n <- length(traits.m2_n[[i]])
			sort(sample(1 : n, ceiling(n * 0.3), replace = FALSE))
					})

		#training model on training test
		forest.test <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			ytrain <- traits.m2_n[[i]][-test.ind[[i]]]
			ytest <- traits.m2_n[[i]][test.ind[[i]]]
			mod <- randomForest(x = geno[names(ytrain), ], y = ytrain, importance = TRUE, xtest = geno[names(ytest), ], ytest = ytest, ntree = 700, do.trace = 100)
			pred <- mod$test$predicted
			return(list(mod = mod, pred = pred))
				}  
				
		#predictions on prediction portion of data
		forest.pred <- lapply(forest.test, '[[', 2)
		names(forest.pred) <- names(traits.m2_n)
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2_n[[i]][test.ind[[i]]], forest.pred[[i]], rsq = TRUE))
		
		rf.rsq[, k] <- rsq
	}

	rf.mcmc.rsq[, , u] <- rf.rsq
	rf.rsq <- matrix(0, nrow = 46, ncol = 7, dimnames = list(names(traits.m2), anames))
}


#===============================================================================================
#LASSO

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 7){
		a <- avec[k]
		aname <- anames[k]

		#crating new phenotypic output by randomly selecting a*100% of data; new random seed is set at teach MCMC iteration
		traits.m2_n <- lapply(traits.m, FUN = function(x){
					n <- length(x)
					set.seed(u * 100)
					ind <- sort(sample(1 : n, round(n * a))) 
					x <- x[ind]
					out <- as.numeric(scale(x))
					names(out) <- names(x)
					out
				})
		#creating new test-train folds (random 30% - 70% split for each trait)			
		test.ind <- lapply(1 : 46, FUN = function(i){
			set.seed(100)
			n <- length(traits.m2_n[[i]])
			sort(sample(1 : n, ceiling(n * 0.3), replace = FALSE))
				})
				
	#training model on training test
	lasso.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			y <- traits.m2_n[[i]][-test.ind[[i]]]
			mod <- cv.glmnet(geno[names(y), ], y, family = 'gaussian', alpha = 1, standardize = FALSE, nfolds = 10, parallel = F)
			mod
	}


		lasso.pred <- sapply(1 : 46, FUN = function(i) predict(lasso.mod[[i]], geno[names(traits.m2_n[[i]][test.ind[[i]]]), ]))
		names(lasso.pred) <- names(traits.m2_n)
		
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2_n[[i]][test.ind[[i]]], lasso.pred[[i]], rsq = TRUE))
		lasso.rsq[, k] <- rsq
	}

		lasso.mcmc.rsq[, , u] <- lasso.rsq
		lasso.rsq <- matrix(0, nrow = 46, ncol = 7, dimnames = list(names(traits.m2), anames))
}



#===============================================================================================
#RIDGE

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 7){
		a <- avec[k]
		aname <- anames[k]

		#crating new phenotypic output by randomly selecting a*100% of data; new random seed is set at teach MCMC iteration
		traits.m2_n <- lapply(traits.m, FUN = function(x){
					n <- length(x)
					set.seed(u * 100)
					ind <- sort(sample(1 : n, round(n * a))) 
					x <- x[ind]
					out <- as.numeric(scale(x))
					names(out) <- names(x)
					out
				})
		#creating new test-train folds (random 30% - 70% split for each trait)			
		test.ind <- lapply(1 : 46, FUN = function(i){
			set.seed(100)
			n <- length(traits.m2_n[[i]])
			sort(sample(1 : n, ceiling(n * 0.3), replace = FALSE))
				})
				
	#training model on training test
	ridge.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			y <- traits.m2_n[[i]][-test.ind[[i]]]
			mod <- cv.glmnet(geno[names(y), ], y, family = 'gaussian', alpha = 0, standardize = FALSE, nfolds = 10, parallel = F)
			mod
	}


		ridge.pred <- sapply(1 : 46, FUN = function(i) predict(ridge.mod[[i]], geno[names(traits.m2_n[[i]][test.ind[[i]]]), ]))
		names(ridge.pred) <- names(traits.m2_n)
		
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2_n[[i]][test.ind[[i]]], ridge.pred[[i]], rsq = TRUE))
		ridge.rsq[, k] <- rsq
	}

		ridge.mcmc.rsq[, , u] <- ridge.rsq
		ridge.rsq <- matrix(0, nrow = 46, ncol = 7, dimnames = list(names(traits.m2), anames))
}


#------------------------------------------------------------

lasso.sample.rsq <- apply(lasso.mcmc.rsq, c(1, 2), mean)
gbm.sample.rsq <- apply(gbm.mcmc.rsq, c(1, 2), mean)
ridge.sample.rsq <- apply(ridge.mcmc.rsq, c(1, 2), mean)
rf.sample.rsq <- apply(rf.mcmc.rsq, c(1, 2), mean)

save(list = c('lasso.sample.rsq', 'gbm.sample.rsq', 'ridge.sample.rsq', 'rf.sample.rsq'), file = 'sample_mcmc_averaged_rsq.rds')
















