#' Convert Entrez IDs to Ensembl IDs
#' 
#' This function takes in a vector of Entrez IDs and return a vector of mapped Ensembl IDs.
#' 
#' @param entrez_ids a vector of numbers or numbers as characters
#' @param species character string, "human" or "mouse" (default = "human")
#' 
#' @returns a vector of mapped Ensembl IDs
#' 
#' @export
entrez_to_ensembl <- function(entrez_ids, species = "human") {
	
	validate_entrez(entrez_ids)
	
	if (!species %in% c("human", "mouse")) {
		stop("`species` argument expects either \"human\" or \"mouse\"")
	}
	
	if (species == "mouse") {
		gene_map <- gene_map_mouse
	} else {
		gene_map <- gene_map_human
	}
	
	entrez_ids <- modernize_entrez_ids(entrez_ids)
	
	ensembl_ids <- gene_map$ensembl_id[match(x = entrez_ids, table = gene_map$entrez_id)]
	
	if (any(is.na(ensembl_ids))) {
		warning("Some Ensembl IDs missing (may be discontinued genes)")
	}
	
	if (any(duplicated(ensembl_ids))) {
		warning("Mapped Ensembl IDs contain duplicated entries")
	}
	
	ensembl_ids
}
