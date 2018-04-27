#===============================================================================================#
#       MCMC procedure for the attribute number experiment (deleting some % of attributes)		#
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
#vector of proportion of markers to be used (i.e. if use 0.01 of the data, then 0.99 is removed)
avec <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9)
anames <- c('01', '05', '10', '20', '30', '40', '50', '75', '90')

#matrices to be filled with results in the for loops
gbm.subrsq <- rf.subrsq <- lasso.subrsq <- ridge.subrsq <- matrix(0, nrow = 46, ncol = 9, dimnames = list(names(traits.m2), anames))
gbm.mcmc.rsq <- rf.mcmc.rsq <- lasso.mcmc.rsq <- ridge.mcmc.rsq <- array(0, c(46, 9, 10), dimnames = list(names(traits.m2), anames, 1 : 10))


#===============================================================================================
#GBM 

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 9){
		a <- avec[k]
		aname <- anames[k]

		#crating new genotypic output by randomly selecting a*100% of data; new random seed is set at teach MCMC iteration
		set.seed(100 * u)
		ind <- sample(11623, round(11623*a))
		ind <- sort(ind)
		geno.r <- geno[, ind]
			
		#training model on training test
		gbm.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			ytrain <- traits.m2[[i]][-test.ind[[i]]]
			nTrain_ <- 3/4 * length(ytrain)
			mod <- gbm.fit(geno.r[names(ytrain), ], ytrain, distribution = 'gaussian', n.trees = 1000, interaction.depth = 5, n.minobsinnode = 10, shrinkage = 0.01, bag.fraction = 0.5, nTrain = nTrain_, keep.data = FALSE, verbose = TRUE)
			mod
		}

		gbm.pred <- sapply(1 : 46, FUN = function(i) predict(gbm.mod[[i]], geno.r[names(traits.m2[[i]][test.ind[[i]]]), ]))
		names(gbm.pred) <- names(traits.m)
		
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], gbm.pred[[i]], rsq = TRUE))
		gbm.subrsq[, k] <- rsq
	}

	gbm.mcmc.rsq[, , u] <- gbm.subrsq
}

	
#=============================================================================================== 
#RF

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 9){
		a <- avec[k]
		aname <- anames[k]

		set.seed(100 * u)
		ind <- sample(11623, round(11623*a))
		ind <- sort(ind)
		geno.r <- geno[, ind]

		forest.test <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			ytrain <- traits.m2[[i]][-test.ind[[i]]]
			ytest <- traits.m2[[i]][test.ind[[i]]]
			mod <- randomForest(x = geno.r[names(ytrain), ], y = ytrain, importance = TRUE, xtest = geno.r[names(ytest), ], ytest = ytest, ntree = 700, do.trace = 100)
			pred <- mod$test$predicted
			return(list(mod = mod, pred = pred))
				}  

		forest.mod <- lapply(forest.test, '[[', 1)
		names(forest.mod) <- names(traits.m)

		forest.pred <- lapply(forest.test, '[[', 2)
		names(forest.pred) <- names(traits.m)
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], forest.pred[[i]], rsq = TRUE))
		
		rf.subrsq[, k] <- rsq
	}

	rf.mcmc.rsq[, , u] <- rf.subrsq
}


#=============================================================================================== 
#LASSO

#outer MCMC loop
for(u in 1 : 10){

	for(k in 1 : 9){
		a <- avec[k]
		aname <- anames[k]

		set.seed(100 * u)
		ind <- sample(11623, round(11623*a))
		ind <- sort(ind)
		geno.r <- geno[, ind]

		lasso.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			y <-traits.m2[[i]][-test.ind[[i]]]
			mod <- cv.glmnet(geno.r[names(y), ], y, family = 'gaussian', alpha = 1, standardize = FALSE, nfolds = 10, parallel = F)
			mod
	}

		lasso.pred <- sapply(1 : 46, FUN = function(i) predict(lasso.mod[[i]], geno.r[names(traits.m2[[i]][test.ind[[i]]]), ]))
		names(lasso.pred) <- names(traits.m)
		
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], lasso.pred[[i]], rsq = TRUE))
		lasso.subrsq[, k] <- rsq
	}
	lasso.mcmc.rsq[, , u] <- lasso.subrsq
}


#=============================================================================
#RIDGE
library(glmnet)
for(u in 1 : 10){

	for(k in 1 : 9){
		a <- avec[k]
		aname <- anames[k]

		set.seed(100 * u)
		ind <- sample(11623, round(11623*a))
		ind <- sort(ind)
		geno.r <- geno[, ind]

	ridge.mod <- foreach(i = 1 : 46, .options.RNG = 100) %dorng% {
			y <-traits.m2[[i]][-test.ind[[i]]]
			mod <- cv.glmnet(geno.r[names(y), ], y, family = 'gaussian', alpha = 0, standardize = FALSE, nfolds = 10, parallel = F)
			mod
	}

		ridge.pred <- sapply(1 : 46, FUN = function(i) predict(ridge.mod[[i]], geno.r[names(traits.m2[[i]][test.ind[[i]]]), ]))
		names(ridge.pred) <- names(traits.m)
		
		rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], ridge.pred[[i]], rsq = TRUE))
		ridge.subrsq[, k] <- rsq
	}

	ridge.mcmc.rsq[, , u] <- ridge.subrsq
}


#------------------------------------------------------------

lasso.attr.rsq <- apply(lasso.mcmc.rsq, c(1, 2), mean)
gbm.attr.rsq <- apply(gbm.mcmc.rsq, c(1, 2), mean)
ridge.attr.rsq <- apply(ridge.mcmc.rsq, c(1, 2), mean)
rf.attr.rsq <- apply(rf.mcmc.rsq, c(1, 2), mean)

save(list = c('lasso.attr.rsq', 'gbm.attr.rsq', 'ridge.attr.rsq', 'rf.attr.rsq'), file = 'attr_mcmc_averaged_rsq.rds')






















# #----------------------------------------------------------------------------------------------------------
# #cycle through mcmc runs of a procedure and read in results for varying values of a parameter, given
# #path = location of runs, run = how the runs are called; default - 'run', n = number of runs
# library(abind)

# read_mcmc <- function(path, n, run = 'run'){
	# lapply(paste(path, '/', run , 1 : 10, sep = ''), FUN = function(folder){
		# files_ <- paste(folder, sort(list.files(folder)), sep = '/') #files in run i
		# xx <- lapply(files_, FUN = function(file_) readRDS(file_)) #results for different parameter values
		# res <- sapply(xx, FUN = function(x){
			# rsq <- sapply(1 : 46, FUN = function(i) calc.mse(traits.m2[[i]][test.ind[[i]]], x[[i]], rsq = TRUE))
			# rsq
					# })
			# res
						# }
				# )
# }


# gbm.mcmc.rsq <- read_mcmc('model_output_sub/full/gbm/pred', 10)
# rf.mcmc.rsq <- read_mcmc('model_output_sub/full/forest/pred', 10)
lasso.mcmc.rsq <- read_mcmc('model_output_sub/full/lasso/pred', 10)
# ridge.mcmc.rsq <- read_mcmc('model_output_sub/full/ridge/pred', 10)

# gbm.mcmc.rsq <- abind(gbm.mcmc.rsq, along = 3)
# rf.mcmc.rsq <- abind(rf.mcmc.rsq, along = 3)
lasso.mcmc.rsq <- abind(lasso.mcmc.rsq, along = 3)
# ridge.mcmc.rsq <- abind(ridge.mcmc.rsq, along = 3)

dimnames(gbm.mcmc.rsq) <- dimnames(rf.mcmc.rsq) <- dimnames(lasso.mcmc.rsq) <- dimnames(ridge.mcmc.rsq) <- list(names(traits.m), c('01', '05', '10', '20', '30', '40', '50', '75', '90'))

# saveRDS(gbm.mcmc.rsq, file = 'model_output_sub/full/gbm/gbm_mcmc_rsq.rds')
# saveRDS(rf.mcmc.rsq, file = 'model_output_sub/full/forest/rf_mcmc_rsq.rds')
saveRDS(lasso.mcmc.rsq, file = 'model_output_sub/full/lasso/lasso_mcmc_rsq.rds')
# saveRDS(ridge.mcmc.rsq, file = 'model_output_sub/full/ridge/ridge_mcmc_rsq.rds')











#====================================================================================================================



library(scales)

avec <- c(0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1)

#Rsq results for test sets
test_rsq <- readRDS('model_output/full/full_test_rsq.rds')
gbm_mcmc_rsq <- readRDS('model_output_sample2/gbm/gbm_mcmc_rsq.rds')
gbm_mcmc_rsq <- apply(gbm_mcmc_rsq, c(1, 2), mean)


	dat <- gbm_mcmc_rsq
	dat <- cbind(dat, test_rsq[, 'gbm'])/test_rsq[, 'gbm']
	plot(c(0.15, 1), c(-0.2, 1.2), type = 'n', xlab = 'Fraction of data used', ylab = 'Fraction of variance explained relative to full data', main = titles_[m], xaxs = 'i', yaxs = 'i')
	for(i in 1 : 46) points(avec, dat[i, ], col = alpha('navyblue', 0.15), type = 'l')
	lines(avec, colMeans(dat), col =  'navyblue', lwd = 2)
	abline(h = c(0, 1), lty = 2, col = 'grey30')


plot(c(0.15, 1), c(-0.2, 1.2), type = 'n', xlab = 'Fraction of data used', ylab = 'Fraction of variance explained relative to full data', main = 'Comparison of themofds', xaxs = 'i', yaxs = 'i')
for(m in c('lasso', 'ridge', 'gbm', 'rf')){
	dat <- sample_mcmc_rsq[[m]]
	dat <- cbind(dat, test_rsq[, m])/test_rsq[, m]
	lines(avec, colMeans(dat), col = alpha(cols_[m], 0.5), lwd = 1.5)
	abline(h = c(0, 1), lty = 2, col = 'grey30')
}
legend(0.75, 0.4, c('Lasso', 'Ridge', 'GBM', 'RF'), col = cols_, lty = 1, lwd = 2, bty = 'n', cex = 1)






















