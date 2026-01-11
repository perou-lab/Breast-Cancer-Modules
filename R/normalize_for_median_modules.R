#' Performs normalization procedure for median-based modules
#' 
#' This function normalizes a data matrix
#' 
#' @param x a matrix, already log-transformed or similar
#' @param filter_genes logical, whether to filter genes prior to median module calculation (default = TRUE)
#' @param filter_gene_thresh a number between 0 and 1, fraction of samples each gene must be detected in to pass filter. This parameter does nothing if filter_genes is set to FALSE (default = 0.7)
#' @param row_median_center logical, whether to median-center rows prior to median module calculation (default = TRUE)
#' @param column_standardize logical, whether to column-standardize data matrix prior to median module calculation (default = TRUE)
#' 
#' @returns a matrix of normalized values
#' 
#' @export
normalize_for_median_modules <- function(x, 
										 filter_genes = TRUE, 
										 filter_gene_thresh = 0.7, 
										 row_median_center = TRUE, 
										 column_standardize = TRUE) {
	
	if (filter_genes) {
		genes_keep <- rowSums(x != 0, na.rm = TRUE) >= filter_gene_thresh * ncol(x)
		x <- x[genes_keep, , drop = FALSE]
	}
	
	if (row_median_center) {
		# Vectorized row median centering using matrixStats
		if (requireNamespace("matrixStats", quietly = TRUE)) {
			row_medians <- matrixStats::rowMedians(as.matrix(x), na.rm = TRUE)
		} else {
			row_medians <- apply(x, 1, median, na.rm = TRUE)
		}
		x <- x - row_medians
	}
	
	if (column_standardize) {
		# Vectorized column standardization - scale() works column-wise by default
		ensembl_ids <- rownames(x)
		x <- scale(x)
		rownames(x) <- ensembl_ids
	}
	
	return(x)
}
