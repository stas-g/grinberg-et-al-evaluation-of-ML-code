#---------------------------------------------------------------------------------------------------------------------------------------#
#                         						 PROCESSING AND SAVING DATA FOR FURTHER USE                      						#
#---------------------------------------------------------------------------------------------------------------------------------------#

load("cross.Rdata")
load("pheno_raw.Rdata")

#segregant names
seg_names <- names(pheno_raw[[2]])

#genotypic matrix
geno <- do.call('cbind', lapply(cross$geno, function(x) { x$data })) - 1
rownames(geno) <- seg_names

#phenotypic matrix
pheno <- cross$pheno
#phenotypic matrix standardised (mean 0, variance 1)
pheno2 <- scale(pheno)
rownames(pheno) <- rownames(pheno2) <- seg_names

#list of 46 phenotypes with missing values removed 
traits.m <- apply(pheno, 2, FUN = function(x) x[!is.na(x)])
#standardised version of traits.m (mean 0, variance 1)
traits.m2 <- apply(pheno2, 2, FUN = function(x) x[!is.na(x)])

#list of test indices (30% of data for each phenotype)
test.ind <- lapply(1 : 46, FUN = function(i){
	set.seed(100)
	n <- length(traits.m[[i]])
	sort(sample(1 : n, ceiling(n * 0.3), replace = FALSE))
})


test.ind2 <- lapply(1 : 46, FUN = function(i) which(!is.na(pheno[, i]))[test.ind[[i]]])

#--------------------------------------------------------------------------------------------------------
#creating folds for 10-fold cross-validation (same as folds used for re-running Bloom's method)

grps <- c(rep('A', 101), rep('B', 101), rep('C', 101), rep('D', 101), rep('E', 101),
                  rep('F', 101),rep('G', 101),rep('H', 101),rep('I', 100),rep('J', 100) )
set.seed(100) 
sgrps <- sample(grps)

yfolds <- lapply(LETTERS[1 : 10], FUN = function(a) which(sgrps == a))

#--------------------------------------------------------------------------------------------------------
#saving (geno.rds contains genotypic data used in most analysis, yeast_pheno contains various versions of phenotypic data (matrix/list, scaled/non-scaled), folds for cross-validation and test-train indices

saveRDS(geno, file = 'geno.rds')
save(list = c('traits.m', 'traits.m2', 'test.ind', 'pheno', 'pheno2', 'yfolds'), file = 'yeast_pheno.RData')



