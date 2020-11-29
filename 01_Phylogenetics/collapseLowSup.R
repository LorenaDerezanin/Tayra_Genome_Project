collapseLowSup <- function(phy, thres = 50){
	tocollapse <- which(as.numeric(phy$node.label) < thres) + Ntip(phy)
	tocollapse <- which(phy$edge[,2] %in% tocollapse)
	if(length(tocollapse) > 0){
		phy$edge.length[tocollapse] <- 0
		phy <- di2multi(phy)
		return(phy)
	} else {
		return(phy)
	}
}