#' Convert Data Matrix to Ensembl Matrix
#' 
#' This function takes in a matrix of genes x samples and returns a matrix of mapped Ensembl IDs x samples.
#' 
#' @param x, a matrix of numbers, rows = genes and columns = samples
#' @param species, character string, "human" or "mouse" (default = "human")
#' @param sum_or_mean, a string, method to collapse rows (default = "sum")
#' 
#' @returns a matrix of mapped Ensembl IDs x samples
#' 
#' @export
convert_to_ensembl <- function(x, 
							   species = "human", 
							   sum_or_mean = "sum") {
	
	if (!species %in% c("human", "mouse")) {
		stop("`species` argument expects either \"human\" or \"mouse\"")
	}
	
	x_genes <- row.names(x)
	id_type <- check_id_type(x_genes)
	
	if (species == "mouse") {
		gene_map <- gene_map_mouse
	} else {
		gene_map <- gene_map_human
	}
	
	if (id_type == "Entrez") {
		x_genes_ensembl <- entrez_to_ensembl(entrez_ids = x_genes, species = species)
	} else if (id_type == "Symbols") {
		x_genes_ensembl <- symbol_to_ensembl(gene_symbols = x_genes, species = species)
	} else {
		# Strip version number
		x_genes_ensembl <- gsub(pattern = "\\..*$", replacement = "", x = x_genes)
	}
	
	x_ensembl <- collapse_duplicates(x = x, ids = x_genes_ensembl, sum_or_mean = sum_or_mean)
	
	return(x_ensembl)
	
}
