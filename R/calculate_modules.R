#' Calculate Breast Cancer Modules curated in the Perou Lab
#' 
#' This function takes in a data matrix of genes x samples and returns module scores.
#' 
#' @param x a matrix, preferably of Ensembl IDs (p) x samples (n)
#' @param standardize_output logical, whether to row standardize the output (default = TRUE)
#' @param filter_genes logical, whether to filter genes prior to median module calculation (default = TRUE)
#' @param filter_gene_thresh a number between 0 and 1, fraction of samples each gene must be detected in to pass filter. This parameter does nothing if filter_genes is set to FALSE (default = 0.7)
#' @param row_median_center logical, whether to median-center rows prior to median module calculation (default = TRUE)
#' @param column_standardize logical, whether to column-standardize data matrix prior to median module calculation (default = TRUE)
#' @param log_transform logical, whether to log-transform data matrix prior to median module calculation (default = TRUE)
#' @param include_median_modules logical, whether to calculate median-based modules (default = TRUE)
#' @param include_single_genes logical, whether to calculate single gene modules (default = TRUE)
#' @param include_special_models logical, whether to calculate centroid-based modules and special models (default = TRUE)
#' @param include_prognostic logical, whether to include prognostic models (default = TRUE)
#' @param custom_median_modules list of vectors (optional), e.g. list("CustomModule1" = c("GeneA", "GeneB", "GeneC"), "CustomModule2" = c("GeneD", "GeneE", "GeneF")) 
#' @param input_species character string, "human" or "mouse" (default = "human")
#' 
#' @returns a matrix of module scores (p) x samples (n)
#' 
#' @export
calculate_modules <- function(x,
							  as.is = FALSE, 
							  standardize_output = TRUE,
							  filter_genes = TRUE,
							  filter_gene_thresh = 0.7,
							  row_median_center = TRUE,
							  column_standardize = TRUE,
							  log_transform = TRUE, 
							  include_median_modules = TRUE,
							  include_single_genes = TRUE,
							  include_special_models = TRUE,
							  include_prognostic = TRUE, 
							  custom_median_modules = NULL, 
							  input_species = "human") {
	
	# TODO
	# Separate out conversion-to-Ensembl and microarray-vs-RNA-seq functions from main calculate_modules function
	# Currently buggy and somtimes guesses formatting incorrectly
	if (!as.is) {
	if (all(x < 100) & all(x >= 0)) {
		stop("Dataset appears to already be log-transformed. Please supply normalized counts.")
	} else if (any(x < 0)) {
		warning("Dataset appears to already be standardized. Skipping centering and standardizing steps")
		row_median_center <- FALSE
		column_standardize <- FALSE
		log_transform <- FALSE
	}
	}
	
	
	
	id_type <- check_id_type(row.names(x))
	if (is.na(id_type)) {
		stop("Please check row names of your data matrix. Row names don't appear to match any of the 3 acceptable input values (Entrez IDs, Ensembl IDs, and Gene Symbols)")
	}
	
	# if (id_type != "Ensembl") {
		x <- convert_to_ensembl(x = x, 
								species = input_species, 
								sum_or_mean = ifelse(log_transform, "sum", "mean"))
	# }
	
	if (log_transform & !as.is) {
		x <- log2(x + 1)
	}
	
	# Collect results in a list to avoid multiple rbind() calls
	module_list <- list()
	
	if (include_median_modules) {
		x_median_modules <- calculate_median_modules(x = x,
													 filter_genes = filter_genes,
													 filter_gene_thresh = filter_gene_thresh,
													 row_median_center = row_median_center,
													 column_standardize = column_standardize, 
													 custom_median_modules = custom_median_modules, 
													 input_species = input_species)
		module_list <- append(module_list, list(x_median_modules))
	}
	
	if (include_single_genes) {
		x_single_genes <- calculate_single_genes(x = x, input_species = input_species)
		module_list <- append(module_list, list(x_single_genes))
	}
	
	if (include_special_models) {
		x_centroids <- calculate_centroid_scores(x = x, input_species = input_species)
		x_special <- calculate_special_models(x = x, input_species = input_species)
		module_list <- append(module_list, list(x_centroids, x_special))
	}
	
	# Single rbind operation is much faster than multiple rbind calls
	if (length(module_list) > 0) {
		x_modules <- do.call(rbind, module_list)
	} else {
		x_modules <- x[FALSE, ]
	}
	
	if (!include_prognostic) {
		prognostic_models <- c("ROR_P_Model_JCO.2009_PMID.19204204", 
							   "GHI_RS_Model_NEJM.2004_PMID.15591335", 
							   "Pcorr_NKI70_Good_Correlation_Nature.2002_PMID.11823860", 
							   "TNBC_Clinically_Relevant_Good.26_BCR.2011_PMID.21978456", 
							   "TNBC_Clinically_Relevant_Poor.230_BCR.2011_PMID.21978456", 
							   "TNBC_Clinically_Relevant_Poor.26_BCR.2011_PMID.21978456")
		x_modules <- x_modules[-which(row.names(x_modules) %in% prognostic_models), ]
	}
	
	if (standardize_output) {
		# Ensure x_modules is a matrix for efficient operations
		if (!is.matrix(x_modules)) {
			x_modules <- as.matrix(x_modules)
		}
		# Vectorized row standardization
		x_modules <- t(scale(t(x_modules)))
		x_modules[is.na(x_modules)] <- 0
	}
	
	x_modules
}
