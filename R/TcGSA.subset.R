TcGSA.subset <- 
	function(tcgsa, list_tcgsa_subset){
		index_subset <- which(tcgsa$GeneSets_gmt$geneset.names %in% list_tcgsa_subset)
		tcgsa_subset <- tcgsa
		tcgsa_subset$fit <- tcgsa_subset$fit[index_subset,]
		tcgsa_subset$GeneSets_gmt$genesets <- tcgsa_subset$GeneSets_gmt$genesets[index_subset]
		tcgsa_subset$GeneSets_gmt$geneset.names <- tcgsa_subset$GeneSets_gmt$geneset.names[index_subset]
		tcgsa_subset$GeneSets_gmt$geneset.descriptions <- tcgsa_subset$GeneSets_gmt$geneset.descriptions[index_subset]
		tcgsa_subset$Estimations <- tcgsa_subset$Estimations[index_subset]
		return(tcgsa_subset)
	}