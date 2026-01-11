# LUMA_HER2E_LTS_JCI.2020_PMID.32573490
luma_her2_score <- function(x, input_species = "human") {
	
	lumA_centroid <- centroids[["Scorr_LumA_Correlation_JCO.2009_PMID.19204204"]][[input_species]]
	her2_centroid <- centroids[["Scorr_Her2_Correlation_JCO.2009_PMID.19204204"]][[input_species]]
	
	genes_keep <- intersect(row.names(x), union(row.names(lumA_centroid), row.names(her2_centroid)))
	
	x <- x[genes_keep, ]
	lumA_centroid <- lumA_centroid[genes_keep, 1]
	her2_centroid <- her2_centroid[genes_keep, 1]
	
	x <- t(scale(t(x)))
	x[is.na(x)] <- 0
	
	return(list("LUMA_HER2E_LTS_JCI.2020_PMID.32573490" = apply(X = x, MARGIN = 2, FUN = function(sample_i) {
		cor(x = sample_i, y = lumA_centroid, use = "complete.obs", method = "spearman") -
			cor(x = sample_i, y = her2_centroid, use = "complete.obs", method = "spearman")
	})))
}