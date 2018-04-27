#17.02.2015
load("data/pheno_raw.Rdata")
load("data/cross.Rdata")
source( "R_code/QTL_mappingFx.R")

library(qtl)

library(doRNG)
library(doMC)
registerDoMC(cores = 46)

#-------------------------------------------------------------------------------
#calculate R-squared for a vector of observed and a vector of predicted values or for matrices of observed and predicted values, where each column corresponds to a different variable
calc.rsq <- function(obs, pred){
						if(!is.matrix(obs) | !is.matrix(pred)){
													obs <- as.matrix(obs)
													pred <- as.matrix(pred)
																}
						rs <- (obs - pred)^2						
						meanmat <- matrix(rep(colMeans(obs, na.rm = TRUE), nrow(obs)), nrow = nrow(obs), byrow = TRUE)
						1 - colSums(rs, na.rm = TRUE)/colSums((obs - meanmat)^2, na.rm = TRUE)
									}

#-------------------------------------------------------------------------------
#preliminaries

# pheno <- scale(cross$pheno)
geno <- do.call('cbind', lapply(cross$geno, function(x) { x$data })) - 1
rownames(geno) <- names(pheno_raw[[2]])
pheno <- cross$pheno

mindex.split <- getMarkerIndexSplit(cross)
chr.mindex.offset <- sapply(mindex.split, min) - 1

grps <- c(rep('A',101),rep('B',101),rep('C',101),rep('D',101),rep('E',101),
                  rep('F',101),rep('G',101),rep('H',101),rep('I',100),rep('J',100) )
set.seed(100) #set random seed 
sgrps <- sample(grps)
folds <- lapply(LETTERS[1 : 10], FUN = function(a) which(sgrps == a)) #create CV folds corresponding to sgrps


#------------------------------
#Calculating and saving list of significant peaks for each trait for each fold like in the Nature paper(using FDR thresholds pre-calculated using all the data)

Apeak.index.cv_bloom <- foreach(gg = unique(grps), .options.RNG = 200) %dorng% {
							print(gg)
							nAcross <- subset(cross, ind = (sgrps != gg))
							Across  <- subset(cross, ind = (sgrps == gg))

							Agdata         <-  extractGenotype(nAcross)
							Apdata.01      <-  extractScaledPhenotype(nAcross, TRUE)
							An.pheno       <-  countStrainsPerTrait(nAcross$pheno)
							ALODS.01       <-  get.LOD.by.COR(An.pheno, Apdata.01, Agdata)
							ALODS.01s      <-  LODmatrix.2.scanone(ALODS.01, nAcross)
							Apeaklist.01   <-  getChrPeaks(mindex.split, chr.mindex.offset, ALODS.01) 
							ApeakArray.01  <-  getPeakArray(Apeaklist.01, 2.69)
							print('step1')

							Apdata.02      <- getPhenoResids(Apdata.01, Agdata, ApeakArray.01) 
							ALODS.02       <- get.LOD.by.COR(An.pheno, Apdata.02, Agdata)
							ALODS.02s      <- LODmatrix.2.scanone(ALODS.02, nAcross,ALODS.01s)
							Apeaklist.02   <- getChrPeaks(mindex.split, chr.mindex.offset, ALODS.02) 
							ApeakArray.02  <- getPeakArray(Apeaklist.02, 2.88)
							print('step2')

							Apdata.03      <- getPhenoResids(Apdata.02, Agdata, ApeakArray.02) 
							ALODS.03       <- get.LOD.by.COR(An.pheno, Apdata.03, Agdata)
							ALODS.03s      <- LODmatrix.2.scanone(ALODS.03, nAcross,ALODS.01s)
							Apeaklist.03   <- getChrPeaks(mindex.split, chr.mindex.offset, ALODS.03) 
							ApeakArray.03  <- getPeakArray(Apeaklist.03, 3.75)
							print('step3')

							Apdata.04      <- getPhenoResids(Apdata.03, Agdata, ApeakArray.03) 
							ALODS.04       <- get.LOD.by.COR(An.pheno, Apdata.04, Agdata)
							ALODS.04s      <- LODmatrix.2.scanone(ALODS.04, nAcross,ALODS.01s)
							Apeaklist.04   <- getChrPeaks(mindex.split, chr.mindex.offset, ALODS.04) 
							ApeakArray.04  <- getPeakArray(Apeaklist.04, 4.63)
							print('step4')

							ApA1 <- cbind(ApeakArray.01, 'J1')
							ApA2 <- cbind(ApeakArray.02, 'J2')
							ApA3 <- cbind(ApeakArray.03, 'J3')
							if(!is.null(ApeakArray.04)){ ApA4 <- cbind(ApeakArray.04, 'J4') }

							names(ApA1)[3] <- 'jump'
							names(ApA2)[3] <- 'jump'
							names(ApA3)[3] <- 'jump'
							if(!is.null(ApeakArray.04)){ names(ApA4)[3] <- 'jump' }
							if(!is.null(ApeakArray.04)) { p.df = rbind(ApA1, ApA2, ApA3, ApA4) } else {
								p.df = rbind(ApA1, ApA2, ApA3) }

							Apeak.index <- data.frame(p.df)
							Apeak.index <- Apeak.index[order(Apeak.index$trait, Apeak.index$markerIndex),]
							Apeak.index <- split(Apeak.index, Apeak.index$trait)	
							Apeak.index
						}

saveRDS(Apeak.index.cv, 'Apeak_index_cv_bloom.rds')


#------------------------------
#Calculating and saving list of significant peaks for each trait for each fold (calculating FDR thresholds using only data in the 9 training folds at each step of CV)

Apeak.index.cv_nfg <- foreach(gg = unique(grps), .options.RNG = 200) %dorng% {
						print(gg)
						nAcross <- subset(cross, ind = (sgrps != gg))
						Across  <- subset(cross, ind = (sgrps == gg))

						Agdata         <-  extractGenotype(nAcross)
						Apdata.01      <-  extractScaledPhenotype(nAcross, TRUE)
						An.pheno       <-  countStrainsPerTrait(nAcross$pheno)
						ALODS.01       <-  get.LOD.by.COR(An.pheno, Apdata.01, Agdata)
						ALODS.01s      <-  LODmatrix.2.scanone(ALODS.01, nAcross)
						Apeaklist.01   <-  getChrPeaks(mindex.split, chr.mindex.offset, ALODS.01) 
						LODS.01.FDR    <-  getPeakFDR(Apeaklist.01$chr.peaks.lod, Apdata.01, Agdata, 1000)
						LODS.01.FDR    <-  min(as.numeric(names(LODS.01.FDR)[LODS.01.FDR <= 0.05]))
						ApeakArray.01  <-  getPeakArray(Apeaklist.01, LODS.01.FDR)
						#ApeakArray.01  <-  getPeakArray(Apeaklist.01, 2.69)
						print('step1')

						Apdata.02      <- getPhenoResids(Apdata.01, Agdata, ApeakArray.01) 
						ALODS.02       <- get.LOD.by.COR(An.pheno, Apdata.02, Agdata)
						ALODS.02s      <- LODmatrix.2.scanone(ALODS.02, nAcross,ALODS.01s)
						Apeaklist.02   <- getChrPeaks(mindex.split, chr.mindex.offset, ALODS.02)
						LODS.02.FDR    <-  getPeakFDR(Apeaklist.02$chr.peaks.lod, Apdata.02, Agdata, 1000)
						LODS.02.FDR    <-  min(as.numeric(names(LODS.02.FDR)[LODS.02.FDR <= 0.05]))
						ApeakArray.02  <-  getPeakArray(Apeaklist.02, LODS.02.FDR)
						#ApeakArray.02  <- getPeakArray(Apeaklist.02, 2.88)
						print('step2')

						Apdata.03      <- getPhenoResids(Apdata.02, Agdata, ApeakArray.02) 
						ALODS.03       <- get.LOD.by.COR(An.pheno, Apdata.03, Agdata)
						ALODS.03s      <- LODmatrix.2.scanone(ALODS.03, nAcross,ALODS.01s)
						Apeaklist.03   <- getChrPeaks(mindex.split, chr.mindex.offset, ALODS.03) 
						LODS.03.FDR    <-  getPeakFDR(Apeaklist.03$chr.peaks.lod, Apdata.03, Agdata, 1000)
						LODS.03.FDR    <-  min(as.numeric(names(LODS.03.FDR)[LODS.03.FDR <= 0.05]))
						ApeakArray.03  <-  getPeakArray(Apeaklist.03, LODS.03.FDR)
						#ApeakArray.03  <- getPeakArray(Apeaklist.03, 3.75)
						print('step3')

						Apdata.04      <- getPhenoResids(Apdata.03, Agdata, ApeakArray.03) 
						ALODS.04       <- get.LOD.by.COR(An.pheno, Apdata.04, Agdata)
						ALODS.04s      <- LODmatrix.2.scanone(ALODS.04, nAcross,ALODS.01s)
						Apeaklist.04   <- getChrPeaks(mindex.split, chr.mindex.offset, ALODS.04) 
						LODS.04.FDR    <-  getPeakFDR(Apeaklist.04$chr.peaks.lod, Apdata.04, Agdata, 1000)
						LODS.04.FDR    <-  min(as.numeric(names(LODS.04.FDR)[LODS.04.FDR <= 0.05]))
						ApeakArray.04  <-  getPeakArray(Apeaklist.04, LODS.04.FDR)
						#ApeakArray.04  <- getPeakArray(Apeaklist.04, 4.63)
						print('step4')

						ApA1 <- cbind(ApeakArray.01, 'J1')
						ApA2 <- cbind(ApeakArray.02, 'J2')
						ApA3 <- cbind(ApeakArray.03, 'J3')
						if(!is.null(ApeakArray.04)){ ApA4 <- cbind(ApeakArray.04, 'J4') }

						names(ApA1)[3] <- 'jump'
						names(ApA2)[3] <- 'jump'
						names(ApA3)[3] <- 'jump'
						if(!is.null(ApeakArray.04)){ names(ApA4)[3] <- 'jump' }
						if(!is.null(ApeakArray.04)) { p.df = rbind(ApA1, ApA2, ApA3, ApA4) } else {
							p.df = rbind(ApA1, ApA2, ApA3) }

						Apeak.index <- data.frame(p.df)
						Apeak.index <- Apeak.index[order(Apeak.index$trait, Apeak.index$markerIndex),]
						Apeak.index <- split(Apeak.index, Apeak.index$trait)	
						Apeak.index
					}

saveRDS(Apeak.index.cv_nfg, 'Apeak_index_cv_nfg.rds')

#------------------------------------------------------------------------------------------------------
#Exact procedure from the paper. For each of the 10 CV steps:
#at each step peaks and QTLs are identified on the 9/10th of the data (nAcross). then the model (HK-regression) is fit and assessed on the held-out 1/10th of the data (Across).

Apeak.index.cv_bloom <- readRDS('Apeak_index_cv_bloom.rds')

cv10_bloom <- lapply(1 : 10, FUN = function(k){					
					gg <- unique(grps)[k]
					print(gg)
					Apeak.index <- Apeak.index.cv_bloom[[k]]
					nAcross <- subset(cross, ind = (sgrps != gg))
					Across  <- subset(cross, ind = (sgrps == gg))
					Agdata         <-  extractGenotype(nAcross)
					Apdata.01      <-  extractScaledPhenotype(nAcross, TRUE)
					An.pheno       <-  countStrainsPerTrait(nAcross$pheno)
					ALODS.01       <-  get.LOD.by.COR(An.pheno, Apdata.01, Agdata)
					ALODS.01s      <-  LODmatrix.2.scanone(ALODS.01, nAcross)
					
 aQTLS <- foreach(i = 1 : length(Apeak.index), .options.RNG = 200) %dorng% { doQTLModel(i, Apeak.index, nAcross, ALODS.01s, refine = FALSE) }
        names(aQTLS) = names(Apeak.index)
        afQTLs <- lapply(1 : length(Apeak.index), function(i){ fitqtl(Across, pheno.col = i, aQTLS[[i]]$rqtl, method='hk', get.ests = FALSE, dropone = FALSE) })
        names(afQTLs) = names(aQTLS)

        AVEQA <- sapply(afQTLs, function(x) { 
            y <- NA
            tryCatch( 
                {
                 y= as.numeric(x$result.full[,'%var'][1]) },
                   error= function(e){y=NA})
           return(y/100)} )
        AVEQA <- t(data.frame(t(AVEQA)))
       return(list(peak.index = Apeak.index, QTLs = aQTLS, fQTLs = afQTLs, ve = AVEQA))
	}
)

#R^2 for 10 folds for each trait
cv10_bloom_ve <- sapply(cv10_bloom, function(x){x$ve})
#average R^2 over 10 folds
cv10_bloom_m <- rowMeans(cv10_bloom_ve)
#sd of the R^2 over 10 folds
cv10_bloom_sd <- apply(cv10_bloom_ve, 1, sd )
			
			
#------------------------------------------------------------------------------------------------------
#Same procedure as in the paper (above) but using simple multivariate regression instead of Haley-Knott regression. For each of the 10 CV steps:
#identify peaks and QTLs on the 9/10th (nAcross) of the data. Then fit and assess the model (lm) on the held-out 1/10th of the data (Across).

cv10_bloom_lm <- lapply(1 : 10, FUN = function(k){					
					gg <- unique(grps)[k]
					print(gg)
					Apeak.index <- Apeak.index.cv_bloom[[k]]
					nAcross <- subset(cross, ind = (sgrps != gg))
					Across  <- subset(cross, ind = (sgrps == gg))
					pheno.test <- extractScaledPhenotype(Across, TRUE)	
	fits <- foreach(i = 1 : 46) %dorng% {
				qtlMarkers <- unique(Apeak.index[[i]][, 'markerIndex'])
				geno.test <- geno[sgrps == gg, qtlMarkers]
				fit <- lm(pheno.test[, i] ~ geno.test, na.action = na.exclude)
				return(fit)
		}
			names(fits) <- names(pheno.test)
			fits
	})

	
cv10_bloom_lm_ve <- sapply(cv10_bloom_lm, FUN = function(x) sapply(x, FUN = function(z) summary(z)$r.squared))
cv10_bloom_lm_m <- rowMeans(cv10_bloom_lm_ve)
cv10_bloom_lm_sd <- apply(cv10_bloom_lm_ve, 1, sd)			

#------------------------------------------------------------------------------------------------------
#Check that, given a list of significant QTLs, simple multivariate regression's performance is virtually indistinguisable from HK regression:

plot(cv10_bloom_m, cv10_bloom_lm_m, col = 'grey50', xlab = 'nature11867 - fitqtl', ylab = 'nature11867 - lm', cex = 1.2, cex.lab = 1.3)
abline(0, 1, lty = 2, col = 'red')

#hence we use lm instead of HK for our modified procedure (which makes implementing changes easier)	
	
#------------------------------------------------------------------------------------------------------
#Alternative (corrected) CV-scheme. For each of the 10 steps:
#identify peaks and QTLs and fit an lm model using only the 9/10th of the data. Then assess model's performance on the witheld 1/10th portion of the data (Across).

Apeak.index.cv_nfg <- readRDS('Apeak_index_cv_nfg.rds')

cv10_lmxpred <- lapply(1 : 10, FUN = function(k){					
					gg <- unique(grps)[k]
					print(gg)
					Apeak.index <- Apeak.index.cv_nfg[[k]]
					nAcross <- subset(cross, ind = (sgrps != gg))
					Across  <- subset(cross, ind = (sgrps == gg))
					pheno.test <- extractScaledPhenotype(Across, TRUE)	
					pheno.train <- extractScaledPhenotype(nAcross, TRUE)
	pred <- foreach(i = 1 : 46, .combine = cbind) %dorng% {		
						qtlMarkers <- unique(Apeak.index[[i]][, 'markerIndex'])
						geno.train <- geno[sgrps != gg, qtlMarkers]
						geno.test <- as.data.frame(geno[sgrps == gg, qtlMarkers])
		
						dat <- data.frame(pheno = pheno.train[, i], geno.train)
						colnames(dat)[-1] <- colnames(geno.train)
	
						fit <- lm(pheno ~ ., data = dat, na.action = na.exclude)
						pred <- predict(fit, geno.test)
						pred
		}
	colnames(pred) <- colnames(pheno.test)
	pred
})		


#cross-validated predictions for all traits (1008 x 46 matrix)
cv10_lmxpred <- do.call(rbind, cv10_lmxpred)
cv10_lmxpred <- cv10_lmxpred[order(unlist(folds)), ]

#calculating cross-validated R^2
cv10_lmxpred_rsq <- calc.rsq(scale(pheno), cv10_lmxpred)

#calculating R^2 for each fold for each trait separately
cv10_lmxpred_ve <- lapply(folds, FUN = function(x){
                                        obs <- scale(pheno[x, ])
                                        pred <- cv10_lmxpred[x, ]
                        calc.rsq(obs, pred) #see calc.mse function below
							}
								)
								
cv10_lmxpred_ve <- do.call(cbind, cv10_lmxpred_ve)    
# and the corresponding standard errors across folds
cv10_lmxpred_sd <- apply(cv10_lmxpred_ve, 1, sd)

#cv10_lmxpred_rsq and cv10_lmxpred_sd are the values found in Table 1 in the paper

#------------------------------------------------------------------------------------------------------
#Comparing results for original CV and modified CV

plot(cv10_lmxpred_rsq, cv10_bloom_m, col = 'grey50', xlab = 'nature11867 - corrected', ylab = 'nature11867', cex = 1.2, cex.lab = 1.3)
abline(0, 1, lty = 2, col = 'red')

save(list =c('cv10_lmxpred_sd', 'cv10_lmxpred_rsq'))































