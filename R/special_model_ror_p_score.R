# ROR_P_Model_JCO.2009_PMID.19204204
ror_p_score <- function(x, input_species = "human") {
	
	proliferation_score <- proliferation_11_mean(x, input_species = input_species)[[1]]
	
	basal_centroid <- centroids[["Scorr_Basal_Correlation_JCO.2009_PMID.19204204"]][[input_species]]
	lumA_centroid <- centroids[["Scorr_LumA_Correlation_JCO.2009_PMID.19204204"]][[input_species]]
	lumB_centroid <- centroids[["Scorr_LumB_Correlation_JCO.2009_PMID.19204204"]][[input_species]]
	her2_centroid <- centroids[["Scorr_Her2_Correlation_JCO.2009_PMID.19204204"]][[input_species]]
	
	genes_keep <- intersect(row.names(x), row.names(lumA_centroid))
	
	x <- x[genes_keep, ]
	basal_centroid <- basal_centroid[genes_keep, 1]
	lumA_centroid <- lumA_centroid[genes_keep, 1]
	lumB_centroid <- lumB_centroid[genes_keep, 1]
	her2_centroid <- her2_centroid[genes_keep, 1]
	
	x <- t(scale(t(x)))
	x[is.na(x)] <- 0
	
	
	return(list("ROR_P_Model_JCO.2009_PMID.19204204" = apply(X = x, MARGIN = 2, FUN = function(sample_i) {
		-0.001 * cor(x = sample_i, y = basal_centroid, use = "complete.obs", method = "spearman") +
			0.069 * cor(x = sample_i, y = her2_centroid, use = "complete.obs", method = "spearman") +
			-0.095 * cor(x = sample_i, y = lumA_centroid, use = "complete.obs", method = "spearman") +
			0.049 * cor(x = sample_i, y = lumB_centroid, use = "complete.obs", method = "spearman")}) + 
			0.338 * proliferation_score
	))
}