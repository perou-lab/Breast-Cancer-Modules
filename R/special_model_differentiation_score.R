# Differentiation.Score_Model_BCR.2010_PMID.20813035
differentiation_score <- function(x, input_species = "human") {
	centroid_MS <- differentiation_centroid_MS
	centroid_mL <- differentiation_centroid_mL
	
	genes <- centroid_MS$entrez_id
	ensembl_ids <- gene_map_human$ensembl_id[match(x = genes, table = gene_map_human$entrez_id)]
	
	if (input_species == "mouse") {
		ensembl_ids <- babel_ensembl$mouse_ensembl_id[match(x = ensembl_ids, table = babel_ensembl$human_ensembl_id)]
	}
	centroid_MS$ensembl_id <- ensembl_ids
	centroid_mL$ensembl_id <- ensembl_ids
	
	centroid_MS <- centroid_MS[!is.na(centroid_MS$ensembl_id), ]
	centroid_mL <- centroid_mL[!is.na(centroid_mL$ensembl_id), ]
	
	centroid_MS <- centroid_MS[centroid_MS$ensembl_id %in% row.names(x), ]
	centroid_mL <- centroid_mL[centroid_mL$ensembl_id %in% row.names(x), ]
	
	x <- x[centroid_MS$ensembl_id, ]
	
	centroid_MS$centroid_scale <- sign(centroid_MS$centroid) * sqrt(centroid_MS$centroid^2 / sum(centroid_MS$centroid^2))
	centroid_mL$centroid_scale <- sign(centroid_mL$centroid) * sqrt(centroid_mL$centroid^2 / sum(centroid_mL$centroid^2))
	
	x_scale <- apply(X = x, MARGIN = 2, FUN = function(x_i) { sign(x_i) * sqrt(x_i^2 / sum(x_i^2)) })
	
	msproj <- apply(X = x_scale, MARGIN = 2, FUN = function(x_i) { x_i %*% centroid_MS$centroid_scale })
	mlproj <- apply(X = x_scale, MARGIN = 2, FUN = function(x_i) { x_i %*% centroid_mL$centroid_scale })
	
	diff_score <- mlproj - msproj 
	
	return(list("Differentiation.Score_Model_BCR.2010_PMID.20813035" = diff_score))
}