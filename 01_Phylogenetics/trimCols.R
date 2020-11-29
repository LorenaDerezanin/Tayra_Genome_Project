# The following function trims alingment columns (or full codons) with less than a given proportion of data. Assumes that the alignment is of class DNAbin.

trimCols <- function(al, prop, codon = T){
	 mat <- as.character(as.matrix(al))
	 ntax <- nrow(mat)
	 propthres <- 1-prop
	 compliantSites <- apply(mat, 2, function(x){
	 		x <- as.character(x)
			compl <- (length(which(x %in% c("N", "n", "?", "-", "O", "o", "X", "x"))) / ntax) < propthres
			return(compl)
			})
			
	 if(codon){
			codIDs <- rep(1:(length(compliantSites)/3), each = 3)
			codsToKeep <- rep(F, length(compliantSites))
			for(i in 1:max(codIDs)){
			      if(all(compliantSites[which(codIDs == i)])) codsToKeep[which(codIDs == i)] <- T
			}
			al <- al[, as.logical(codsToKeep)]
			return(al)
	 } else {
			al <- al[, as.logical(compliantSites)]
			retunr(al)
	 }
}