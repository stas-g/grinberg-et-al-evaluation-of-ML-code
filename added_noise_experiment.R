#===============================================================================================#
#       MCMC procedure for the noise experiment (adding sample noise to x% of the data)			#
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
#vector of noise levels (i.e. proportion of samples to be contaminated)
noise_lev <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9)
nnames <- c('01', '05', '10', '20', '30', '40', '50', '75', '90')

#matrices to be filled with results in the for loops
gbm.noise.rsq <- rf.noise.rsq <- lasso.noise.rsq <- matrix(0, nrow = 46, ncol = 9, dimnames = list(names(traits.m2), nnames))
gbm.mcmc.rsq <- rf.mcmc.rsq <- lasso.mcmc.rsq <- ridge.mcmc.rsq <- array(0, c(46, 9, 10), dimnames = list(names(traits.m2), nnames, 1 : 10))

#===============================================================================================
#GBM

#outer MCMC loop
for(u in 1 : 10){
	#inner loop running through various noise levels
	for(k in 1 : 9){
		noise <- noise_lev[k]
		nname <- nnames[k]
		
		traits.m2_n <- lapply(traits.m, FUN = function(x){
				n <- length(x)
				set.seed(100 * u)
				ind <- sample(1 : n, round(n * noise)) 
				set.seed(100 * u)
				x[ind] <- x[ind] + 2 * sd(x) * sample(c(-1, 1), round(n * noise), replace = TRUE)
				out <- scale(x)
				names(out) <- names(x)
				out
			})

	gbm.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			ytrain <- traits.m2_n[[i]][-test.ind[[i]]]
			nTrain_ <- 3/4 * length(ytrain)
			mod <- gbm.fit(geno[names(ytrain), ], ytrain, distribution = 'gaussian', n.trees = 1000, interaction.depth = 5, n.minobsinnode = 10, shrinkage = 0.01, bag.fraction = 0.5, nTrain = nTrain_, keep.data = FALSE, verbose = TRUE)
			mod
		}

		gbm.pred <- sapply(1 : 46, FUN = function(i) predict(gbm.mod[[i]], geno[names(traits.m2_n[[i]][test.ind[[i]]]), ]))
		names(gbm.pred) <- names(traits.m2_n)
		
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2_n[[i]][test.ind[[i]]], gbm.pred[[i]], rsq = TRUE))
		gbm.noise.rsq[, k] <- rsq
	}

	gbm.mcmc.rsq[, , u] <- gbm.noise.rsq
}


#===============================================================================================
#RF

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 9){
		noise <- noise_lev[k]
		nname <- nnames[k]

		traits.m2_n <- lapply(traits.m, FUN = function(x){
				n <- length(x)
				set.seed(100 * u)
				ind <- sample(1 : n, round(n * noise)) 
				set.seed(100 * u)
				x[ind] <- x[ind] + 2 * sd(x) * sample(c(-1, 1), round(n * noise), replace = TRUE)
				out <- scale(x)
				names(out) <- names(x)
				out
			})

		forest.test <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			ytrain <- traits.m2_n[[i]][-test.ind[[i]]]
			ytest <- traits.m2_n[[i]][test.ind[[i]]]
			mod <- randomForest(x = geno[names(ytrain), ], y = ytrain, importance = TRUE, xtest = geno[names(ytest), ], ytest = ytest, ntree = 700, do.trace = 100)
			pred <- mod$test$predicted
			return(list(mod = mod, pred = pred))
				}  

		forest.mod <- lapply(forest.test, '[[', 1)
		names(forest.mod) <- names(traits.m2_n)

		forest.pred <- lapply(forest.test, '[[', 2)
		names(forest.pred) <- names(traits.m2_n)
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2_n[[i]][test.ind[[i]]], forest.pred[[i]], rsq = TRUE))
		
		rf.noise.rsq[, k] <- rsq
	}

	rf.mcmc.rsq[, , u] <- rf.noise.rsq
}


#===============================================================================================
#LASSO

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 9){
		noise <- noise_lev[k]
		nname <- nnames[k]

		traits.m2_n <- lapply(traits.m, FUN = function(x){
				n <- length(x)
				set.seed(100 * u)
				ind <- sample(1 : n, round(n * noise)) 
				set.seed(100 * u)
				x[ind] <- x[ind] + 2 * sd(x) * sample(c(-1, 1), round(n * noise), replace = TRUE)
				out <- scale(x)
				names(out) <- names(x)
				out
			})


	lasso.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			y <-traits.m2_n[[i]][-test.ind[[i]]]
			mod <- cv.glmnet(geno[names(y), ], y, family = 'gaussian', alpha = 1, standardize = FALSE, nfolds = 10, parallel = F)
			mod
	}


		lasso.pred <- sapply(1 : 46, FUN = function(i) predict(lasso.mod[[i]], geno[names(traits.m2_n[[i]][test.ind[[i]]]), ]))
		names(lasso.pred) <- names(traits.m2_n)
		
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2_n[[i]][test.ind[[i]]], lasso.pred[[i]], rsq = TRUE))
		lasso.noise.rsq[, k] <- rsq
	}

	lasso.mcmc.rsq[, , u] <- lasso.noise.rsq
}



#===============================================================================================
#RIDGE

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 9){
		noise <- noise_lev[k]
		nname <- nnames[k]

		traits.m2_n <- lapply(traits.m, FUN = function(x){
				n <- length(x)
				set.seed(100 * u)
				ind <- sample(1 : n, round(n * noise)) 
				set.seed(100 * u)
				x[ind] <- x[ind] + 2 * sd(x) * sample(c(-1, 1), round(n * noise), replace = TRUE)
				out <- scale(x)
				names(out) <- names(x)
				out
			})


	ridge.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			y <-traits.m2_n[[i]][-test.ind[[i]]]
			mod <- cv.glmnet(geno[names(y), ], y, family = 'gaussian', alpha = 0, standardize = FALSE, nfolds = 10, parallel = F)
			mod
	}


		ridge.pred <- sapply(1 : 46, FUN = function(i) predict(ridge.mod[[i]], geno[names(traits.m2_n[[i]][test.ind[[i]]]), ]))
		names(ridge.pred) <- names(traits.m2_n)
		
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2_n[[i]][test.ind[[i]]], ridge.pred[[i]], rsq = TRUE))
		ridge.noise.rsq[, k] <- rsq
	}

	ridge.mcmc.rsq[, , u] <- ridge.noise.rsq
}

#------------------------------------------------------------

lasso.noise.rsq <- apply(lasso.mcmc.rsq, c(1, 2), mean)
gbm.noise.rsq <- apply(gbm.mcmc.rsq, c(1, 2), mean)
ridge.noise.rsq <- apply(ridge.mcmc.rsq, c(1, 2), mean)
rf.noise.rsq <- apply(rf.mcmc.rsq, c(1, 2), mean)

save(list = c('lasso.noise.rsq', 'gbm.noise.rsq', 'ridge.noise.rsq', 'rf.noise.rsq'), file = 'noise_mcmc_averaged_rsq.rds')


