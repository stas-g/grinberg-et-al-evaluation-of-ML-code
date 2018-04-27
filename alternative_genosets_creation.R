#---------------------------------------------------------------------------------------------------------------------------------------#
#                                           		CREATING ALTERNATIVE GENOTYPIC DATA SETS                                            #
#---------------------------------------------------------------------------------------------------------------------------------------#
geno <- readRDS('geno.rds')
load("cross.Rdata")
source('functions.r')

library(doRNG)
library(doMC)
registerDoMC(cores = 12)

#----------------------------------------------------------------------------------------------------------------------------------------
#NOTE: reference genome and other associated genomic information for this study can be found on http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/ . You will need to download the .zip file for the R62-1-1_20090218 edition, unpack it and let file_ below be the path to the orf_coding_all_R62-1-1_20090220.fasta file containing coding for all the yeast ORF regions

file_ <- "S288C_reference_genome_R62-1-1_20090218/orf_coding_all_R62-1-1_20090220.fasta"
genes.info <- parse.fasta(file_)[[3]]

#list of data frames containing information on all markers (like position, 
marker.info <- lapply(1 : 16, function(i){
    ans <- parse.markername(colnames(cross$geno[[i]]$data))
    rownames(ans) <- colnames(cross$geno[[i]]$data)
    ans
    } )

#a list of length 16, each element lists all genes and what markers they contain (if any) (takes time to compute)	
gene.marker <- foreach(i = 1 : 16) %dopar% {check.chrm(marker.info[[i]], genes.info[[i]], reverse = TRUE)}	
#removing 'empty' ORFs (containing no markers)
gene.marker <- lapply(gene.marker, FUN = function(x) x[!is.na(x)])

#--------------------------------------------------------------------------------------#
#   11                                                                                 #
#  111     combining intergenic regions together with ORFs into 'genes'                #
#   11	   0-1 boolean values according to whether most markers are 0s or 1s           #
#   11 																				   #
#   11                                                                                 #
#--------------------------------------------------------------------------------------#

#creating 'genes' -- combining ORFs with half of their flanking intergenic regions
new.genes <- 
	lapply(1 : 16, FUN = function(j){
		x <- genes.info[[j]]
		pos <- x$pos
		c.ind <- sapply(2 : length(pos), FUN = function(k){
								ind <- intervals(pos[[k - 1]], pos[[k]], cover = TRUE)	
								if(ind == 0) return(NA) else return(c(k - 1, k)[ind])
									} )
		c.ind <- c.ind[!is.na(c.ind)]
		pos <- c(pos[-c.ind], max(pos[[length(pos)]]) + 300)
		x <- x[-c.ind, ]
		mid <- sapply(1 : (length(pos) - 2), FUN = function(i){
			if(intervals(pos[[i]], pos[[i + 1]], bool = TRUE) == 1){
				mid <- (min(pos[[i + 2]]) - max(pos[[i]]) - 1)/2
				mid <- max(pos[[i]]) + round(mid)
				mid
			} else {
				mid <- (min(pos[[i + 1]]) - max(pos[[i]]) - 1)/2
				mid <- max(pos[[i]]) + round(mid)
				mid
			}
		} 
	)
			mid <- c(0, mid, mid[length(mid)] + 1000)
			x[, 'pos'] <- as.array(list(lapply(2 : length(mid), FUN = function(i) mid[c(i - 1, i)] + c(1, 0))))
			x
	} 
)


#a list of length 16, each element lists all genes that contain any markers and what markers they contain (so 'empty' genes are not listed)	
gene.marker2 <- foreach(i = 1 : 16) %dopar% {check.chrm(marker.info[[i]], new.genes[[i]], reverse = TRUE)}
gene.marker2 <- lapply(gene.marker2, FUN = function(x) x[!is.na(x)])

#counting 0s and 1s in each new gene and assigning the dominating value to it
geno.genes <- lapply(gene.marker2, FUN = function(x){
                do.call(cbind, lapply(x, FUN = function(z){
                            if(length(z) == 1) return(geno[, z]) else{
                            counts <- rowSums(geno[, z])
                            thr <- floor(length(z)/2)
                            counts[counts < thr] <- 0
                            counts[counts >= thr] <- 1
                            return(counts)
      }
  }
         ))
}
       )

geno.genes <- do.call(cbind, geno.genes)


#--------------------------------------------------------------------------------------#
#   222                                                                                #
#  2 22    separating intergenic regions into different variables                      #
#   22	   0-1 boolean values according to whether most markers are 0s or 1s		   #
#  22                                                                                  #
# 22222                                                                                #
#--------------------------------------------------------------------------------------#	   

#creating data.frames describing intergenic regions (i.e. regions between ORFs)
intergene <- 
	lapply(1 : 16, FUN = function(i){
		x <- genes.info[[i]]
		pos <- x$pos
		new.dat <- lapply(2 : length(pos), FUN = function(k){
				if(intervals(pos[[k - 1]], pos[[k]], bool = TRUE) == 1){
					return(NA)
					} else {
							p <- c(max(pos[[k - 1]]) + 1, min(pos[[k]]) - 1)
							n <- paste(rownames(x)[c(k - 1, k)], collapse = '.')
							return(list(name = n, pos = p))
								}
	} )
	new.dat <- new.dat[!is.na(new.dat)]
	new.dat <- data.frame(do.call(rbind, new.dat))
	rownames(new.dat) <- new.dat[, 1]; colnames(new.dat) <- c('name', 'pos')
	new.dat
})

#a list of length 16, each element lists all intergenuc regions that contain any markers and what markers they contain (so 'empty' genes are not listed)	
intergene.marker <- foreach(i = 1 : 16) %dopar% {check.chrm(marker.info[[i]], intergene[[i]], reverse = TRUE)}
intergene.marker <- lapply(intergene.marker, FUN = function(x) x[!is.na(x)])

#counting 0s and 1s in each intergenic region and assigning the dominating value to it
geno.inter <- lapply(intergene.marker, FUN = function(x){
          do.call(cbind, lapply(x, FUN = function(z){
                            if(length(z) == 1) return(geno[, z]) else{
                            counts <- rowSums(geno[, z])
                            thr <- floor(length(z)/2)
                            counts[counts < thr] <- 0
                            counts[counts >= thr] <- 1
                            return(counts)
      }
  }
         ))
}
       )

#counting 0s and 1s in eachORF and assigning the dominating value to it	   
geno.orf <- lapply(gene.marker, FUN = function(x){
          do.call(cbind, lapply(x, FUN = function(z){
                            if(length(z) == 1) return(geno[, z]) else{
                            counts <- rowSums(geno[, z])
                            thr <- floor(length(z)/2)
                            counts[counts < thr] <- 0
                            counts[counts >= thr] <- 1
                            return(counts)
      }
  }
         ))
}
       )	   

geno.orf <- do.call(cbind, geno.orf)
geno.inter <- do.call(cbind, geno.inter)
#combining ORFs and intergenic regions
geno.fuse <- cbind(geno.orf, geno.inter)

#saving...
saveRDS(geno_genes, file = 'geno_genes.rds')
saveRDS(geno_fuse, file = 'geno_fuse.rds')









