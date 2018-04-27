#---------------------------------------------------------------------------------------------------------------------------------------#
#                                           FUNCTIONS FOR PROCESSING AND MANUPULATING GENETIC DATA                                      #
#---------------------------------------------------------------------------------------------------------------------------------------#

#EXTRACT.GENOME: processes fasta file returning a list of sequences with an optional list of corresponding entry names (return.info = TRUE). file = path to fasta file (can also be a .txt file)
extract.genome <- function(file, return.info = FALSE){
    gen.seq <- readLines(file)
	ind <- grep('^>', gen.seq)
	gen.info <- substr(gen.seq[ind], 2, 1000L)
    seq.names <- regmatches(gen.info, regexpr("^(\\S+)", gen.info))	
	ind <- c(ind, length(gen.seq) + 1)
    gen.seq <- lapply(1 : (length(ind) - 1), FUN = function(i){
            paste(gen.seq[(ind[i] + 1) : (ind[i + 1] - 1)], collapse = '')
      })
	names(gen.seq) <- seq.names
    if(return.info == FALSE) { return(gen.seq)
		} else { return(list(gen.seq = gen.seq, gen.info = gen.info)) }
}

#------------------------------------------------------------------------------
#PARSE.FASTA: extracts the sequences and parses the names into a table with chrm, position, name etc. creates a separate file for all the GO (gene ontology) information. tailored for yeast genome with its 16 chromosomes but code is easily generalised and adapted.

parse.fasta <- function(fasta){
  genes <- extract.genome(fasta, TRUE)
  r <- regexec('(\\S+) (\\S+) SGDID:S[0-9]+, Chr ([A-Za-z]+) from ([0-9,-]+) (.+), (\\"(.+)\\")?', genes$gen.info)
  genes.info <- regmatches(genes$gen.info, r)
  ind <- sapply(genes.info, FUN = function(x) grepl('^[A-Z]+$', x[4]))
  genes <- genes[[1]][ind]
  genes.info <- genes.info[ind]
  
  genes.info <- as.data.frame(do.call(rbind, genes.info), stringsAsFactors = FALSE)[, c(2:6, 8)]
  genes.info[, 4] <- as.array(list(lapply(strsplit(gsub("-|,", ' ', genes.info[, 4]), split = ' '), as.numeric)))
  genes.info[, 3] <- sapply(genes.info[, 3], FUN = function(chrm) switch(chrm, I = 1, II = 2, III = 3, IV = 4, V = 5, VI = 6, VII = 7, VIII = 8, IX = 9, X = 10, XI = 11, XII = 12, XIII = 13, XIV = 14, XV = 15, XVI = 16))
  colnames(genes.info) <- c('orf_name', 'com_name', 'chrm', 'pos', 'orf', 'info')
  rownames(genes.info) <- genes.info$orf_name
  
  genes <- split(genes, as.numeric(genes.info$chrm))
  genes.go <- genes.info[, c(1, 5:6)]
  genes.go <- split(genes.go, as.numeric(genes.info$chrm))
  genes.info <- genes.info[, 1:4]
  genes.info <- split(genes.info, as.numeric(genes.info$chrm))
  
  ind <- lapply(genes.info, FUN = function(x) order(sapply(x$pos, min)))
  genes <- lapply(1 : 16, FUN = function(i) genes[[i]][ind[[i]]])
  genes.go <- lapply(1 : 16, FUN = function(i) genes.go[[i]][ind[[i]], ])
  genes.info <- lapply(1 : 16, FUN = function(i) genes.info[[i]][ind[[i]], ])
  
  return(list(seqn = genes, go = genes.go, info = genes.info))
}


#------------------------------------------------------------------------------
#PARSE.MARKERNAME: extract information from marker names (ORF names, position, etc)

parse.markername <- function(mname) {
    m <- data.frame(do.call('rbind', strsplit(mname, '_')), stringsAsFactors = FALSE)
    m[, 1] <- as.numeric(m[, 1])
    m[, 3] <- as.numeric(m[, 3])
    m[, 2] <- as.numeric(substr(m[, 2], 4, 5))
    names(m) <- c('cumpos', 'chrm', 'pos', 'ref', 'alt')
	rownames(m) <- mname
    return(m)  
}

#------------------------------------------------------------------------------
#INGENE: check whether marker m sits inside gene g. m = marker position, g = gene position (position is a vector with length at least two and must be of even length). if intron = FALSE, m sitting inside an intron (gaps within g) returns FALSE. similarly if intron = TRUE, if m sits inside g's intron, TRUE is returned. function assumes that both m and g sit on the same chromosome and it is user's responsibility to ensure so.

ingene <- function(m, g, intron = FALSE){
    if(max(g) < m | m < min(g))  return(FALSE) else{
        if(intron == FALSE){
            if(g[2] < g[1]) g <- rev(g)
            any(m >= g[c(T, F)] & m <= g[c(F, T)])
        } else{
            (m >= min(g) & m <= max(g))
        }
    }
}

#------------------------------------------------------------------------------
#INTERVALS: utility function used in CHECK.CHRM, essentially check whether any two given intervals (represented by vectors of length of at least two and even length) intersect, and (optionally) if they do, whether one is completely contained within the other
# if cover = TRUE then return the interval that is completely covered by the other (1 or 2), or 0 if no coverage
# if bool = TRUE then return 0 if no overlap (incl. coverage) and 1 otherwise
# if bool = FALSE and cover = FALSE, then return 1 if overlap, 2 if coverage and 0 otherwise

intervals <- function(x, y, bool = FALSE, cover = FALSE){
	if(min(x) <= min(y)){int1 <- x; int2 <- y; int.ord <- c(1, 2)
		} else{int1 <- y; int2 <- x; int.ord <- c(2, 1)}
		
		if(min(int1) < min(int2) & max(int1) >= min(int2) & max(int1) < max(int2)){ans <- 1
			}else{
				if((min(int1) < min(int2) & max(int1) >= max(int2)) | (min(int1) == min(int2))){
					ind <- c(max(int1) >= max(int2), max(int2) > max(int1)) 
						ans <- 2
						} else ans <- 0 
							}
			if(cover == TRUE){
					if(ans != 2){return(0)} else{return(int.ord[!ind])}
				} else{
			if(bool == FALSE) return(ans) else{
				ifelse(ans[[1]] == 0, 0, 1)}				
			}
				}

#------------------------------------------------------------------------------
#CHECK.CHRM: given data frames of markers and 'genes' (ORFs, intergenic regions etc) returns a list of all markers and genes that they sit in (if any). if revers = TRUE, then returns a list of all genes and all the markers (if any) that they contain.
#markers.df = marker.info, genes.df = genes.info; markers.df must contian position ('pos'), genes.df must contain position ('pos') and the orf name ('orf_name')

check.chrm <- function(markers.df, genes.df, reverse = FALSE, ...){
	if(reverse == FALSE){
		ans <- lapply(markers.df[, 'pos'], FUN = function(m){
			ind <- sapply(genes.df[, 'pos'], FUN = function(g){
			ingene(m, g)
			})
				if(sum(ind) != 0){ 
				return(genes.df[ind, 'orf_name'])
		} else return(NA)   
	})
	names(ans) <- rownames(markers.df)
	ans
	} else{
		  ans <- lapply(genes.df[, 'pos'], FUN = function(g){
		  ind <- sapply(markers.df[, 'pos'], FUN = function(m){
          ingene(m, g, ...)
  })
      if(sum(ind) != 0){ 
      return(rownames(markers.df)[ind])
   } else return(NA)   
 })
names(ans) <- rownames(genes.df)
ans
}
	}

	
	
#---------------------------------------------------------------------------------------------------------------------------------------#
#                                                        GENERAL PURPOSE FUNCTIONS                                                      #
#---------------------------------------------------------------------------------------------------------------------------------------#

#CALC.MSE: calculate MSE for a vector of observed and a vector of predicted values or for matrices of observed and predicted values, where each column corresponds to a different variable. if rsq = TRUE, calculate R^2 instead

calc.mse <- function(obs, pred, rsq = FALSE){
						if(!is.matrix(obs) | !is.matrix(pred)){
								obs <- as.matrix(obs)
								pred <- as.matrix(pred)
									}
						rs <- (obs - pred)^2						
						if(rsq == FALSE){colMeans(rs, na.rm = TRUE)} else {
							meanmat <- matrix(rep(colMeans(obs, na.rm = TRUE), nrow(obs)), nrow = nrow(obs), byrow = TRUE)
							1 - colSums(rs, na.rm = TRUE)/colSums((obs - meanmat)^2, na.rm = TRUE)
						}
							}

#------------------------------------------------------------------------------
#FOLD.MSE: a version of calc.mse for cross-validated predictions. takes vectors (or matrices with each column being a variable) of observed and predicted (cross-validated) values and a list of folds used for cross-validation. returns a list of overall MSE (R^2 if rsq = TRUE) calculated over all values (i.e. same as calc.mse), mean MSE (R^2) calculated as average of MSE/R^2 values across folds and the corresponding standard deviation

fold.mse <- function(obs, pred, folds, rsq){
		rsqs <- lapply(folds, FUN = function(x){
				if(is.vector(obs)){
					obs <- obs[x]
					pred <- pred[x]
						}else{
							obs <- obs[x, ]
							pred <- pred[x, ]
					}
			calc.mse(obs, pred, rsq = rsq)})
	rsqs <- do.call(rbind, rsqs)	
	sdev <- apply(rsqs, 2, sd)
	m <- colMeans(rsqs)
	rsq <- calc.mse(obs, pred, rsq = rsq)
	list(rsq = rsq, mean = m, sdev = sdev)
}




