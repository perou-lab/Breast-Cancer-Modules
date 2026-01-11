# Chemo_Endocrine_Score_CC.2.LumA_subtract.by_CC.2.Basal_CCR.2017_PMID.27903675
chemo_endocrine_score <- function(x, input_species = "human") {
	
	lumA_centroid <- centroids[["Scorr_LumA_Correlation_JCO.2009_PMID.19204204"]][[input_species]]
	basal_centroid <- centroids[["Scorr_Basal_Correlation_JCO.2009_PMID.19204204"]][[input_species]]
	
	genes_keep <- intersect(row.names(x), union(row.names(lumA_centroid), row.names(basal_centroid)))
	
	x <- x[genes_keep, ]
	lumA_centroid <- lumA_centroid[genes_keep, 1]
	basal_centroid <- basal_centroid[genes_keep, 1]
	
	x <- t(scale(t(x)))
	x[is.na(x)] <- 0
	
	return(list("Chemo_Endocrine_Score_CC.2.LumA_subtract.by_CC.2.Basal_CCR.2017_PMID.27903675" = apply(X = x, MARGIN = 2, FUN = function(sample_i) {
		cor(x = sample_i, y = lumA_centroid, use = "complete.obs", method = "spearman") -
			cor(x = sample_i, y = basal_centroid, use = "complete.obs", method = "spearman")
	})))
}