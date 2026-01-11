#' Calculate Median-Based Modules
#' 
#' This function takes in a data matrix of Ensembl IDs x samples and returns module scores.
#' 
#' @param x a matrix, of Ensembl IDs (p) x samples (n), already log-transformed (if applicable)
#' @param filter_genes logical, whether to filter genes prior to median module calculation (default = TRUE)
#' @param filter_gene_thresh a number between 0 and 1, fraction of samples each gene must be detected in to pass filter. This parameter does nothing if filter_genes is set to FALSE (default = 0.7)
#' @param row_median_center logical, whether to median-center rows prior to median module calculation (default = TRUE)
#' @param column_standardize logical, whether to column-standardize data matrix prior to median module calculation (default = TRUE)
#' @param input_species character string, "human" or "mouse" (default = "human")
#' 
#' @returns a matrix of module scores (p) x samples (n)
#' 
#' @export
calculate_median_modules <- function(x,
									 filter_genes = TRUE,
									 filter_gene_thresh = 0.7,
									 row_median_center = TRUE,
									 column_standardize = TRUE,
									 custom_median_modules = NULL, 
									 input_species = "human") {
	
	x <- normalize_for_median_modules(x = x,
									  filter_genes = filter_genes,
									  filter_gene_thresh = filter_gene_thresh,
									  row_median_center = row_median_center,
									  column_standardize = column_standardize)
	
	x_median_modules <- sapply(X = median_modules, FUN = function(module_i) {
		
		module_i_set <- ifelse(input_species == "mouse", list(module_i$mouse$ensembl_id), list(module_i$human$ensembl_id))[[1]]
		
		matched_genes <- intersect(module_i_set, rownames(x))
		
		if (length(matched_genes) == 0) {
			return(rep(NA, ncol(x)))
		}
		module_i_matrix <- x[matched_genes, , drop = FALSE]
		# x is already a matrix, no need for as.matrix() conversion
		if (requireNamespace("matrixStats", quietly = TRUE)) {
			matrixStats::colMedians(module_i_matrix, na.rm = TRUE)
		} else {
			apply(module_i_matrix, 2, median, na.rm = TRUE)
		}
	}) |> t()
	
	
	if (!is.null(custom_median_modules)) {
		# TODO
		# currently assumes that the input matrix species is the same as the custom module list
		x_custom_median_modules <- sapply(X = custom_median_modules, FUN = function(module_i) {
			if (check_id_type(module_i) == "Entrez") {
				module_i <- entrez_to_ensembl(entrez_ids = module_i, species = input_species)
			} else if (check_id_type(module_i) == "Symbols") {
				module_i <- symbol_to_ensembl(gene_symbols = module_i, species = input_species)
			}
			
		matched_genes <- intersect(module_i, rownames(x))
		
		if (length(matched_genes) == 0) {
			return(rep(NA, ncol(x)))
		}
		module_i_matrix <- x[matched_genes, , drop = FALSE]
		# x is already a matrix, no need for as.matrix() conversion
		if (requireNamespace("matrixStats", quietly = TRUE)) {
			matrixStats::colMedians(module_i_matrix, na.rm = TRUE)
		} else {
			apply(module_i_matrix, 2, median, na.rm = TRUE)
		}
	}) |> t()
	
	x_median_modules <- rbind(x_median_modules, x_custom_median_modules)
	}
	
	return(x_median_modules)
}
