# REPLICATION_STRESS_MODEL_Cell.Rep.2018_PMID.29768207
replication_stress <- function(x, input_species = "human") {
	rs_model <- replication_stress_model
	if (input_species == "mouse") {
		rs_model$EnsemblID <- babel_ensembl$mouse_ensembl_id[match(x = rs_model$EnsemblID, table = babel_ensembl$human_ensembl_id)]
	}
	genes_keep <- intersect(row.names(x), rs_model$EnsemblID)
	rs_model <- rs_model[rs_model$EnsemblID %in% genes_keep & !is.na(rs_model$EnsemblID), ]
	x <- x[as.character(rs_model$EnsemblID), ]
	
	x <- t(scale(t(x)))
	x[is.na(x)] <- 0
	
	replication_stress_score <- apply(X = x, MARGIN = 2, FUN = function(sample_i) {
		sum(rs_model$Coefficient * sample_i)
	})
	
	return(list("REPLICATION_STRESS_MODEL_Cell.Rep.2018_PMID.29768207" = replication_stress_score))
}