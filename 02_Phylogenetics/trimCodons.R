# This function trims codons with heterozigosity above a given threshold.

trimCodons <- function(al, threshet = 0.5){
	   aaal <- trans(al)
	   ncaaal <- ncol(aaal)
	   toremove <- which(sapply(1:ncaaal, function(x) max(table(as.character(aaal[,x]))) / nrow(aaal)) < threshet)*3
	   if(length(toremove) > 0) al <- al[,-c(toremove, toremove-1, toremove-2)]
	   return(al)
}