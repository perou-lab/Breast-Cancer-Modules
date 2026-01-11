#' Calculate Single Gene Modules
#' 
#' This function takes in a data matrix of Ensembl IDs x samples and returns module scores.
#' 
#' @param x a matrix, of Ensembl IDs (p) x samples (n), already log-transformed
#' @param single_genes a vector of Gene Symbols
#' @param input_species character string, "human" or "mouse" (default = "human")
#' 
#' @returns a matrix of module scores (p) x samples (n)
#' 
#' @export
calculate_single_genes <- function(x,
								   input_species = "human") {

	x_single_genes <- sapply(X = single_genes[[input_species]]$ensembl_id, FUN = function(i) {
		if (!(i %in% rownames(x))) {
			return(rep(NA, ncol(x)))
		}
		
		unlist(x[i, ])
	}) |> t()
	rownames(x_single_genes) <- single_genes[[input_species]]$module_name

	x_single_genes
}
