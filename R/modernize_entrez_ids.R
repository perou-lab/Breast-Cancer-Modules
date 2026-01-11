#' Modernizes Entrez IDs
#' 
#' This function takes in vector of Entrez IDs and returns a vector of updated Entrez IDs based on NCBI's `gene_history.gz` file.
#' 
#' @param entrez_ids a vector of Entrez IDs
#' @param species character string, "human" or "mouse" (default = "human")
#' 
#' @returns a vector of "modernized" Entrez IDs
#' 
#' @export
modernize_entrez_ids <- function(entrez_ids, species = "human") {
	
	validate_entrez(entrez_ids)
	
	if (species == "mouse") {
		history_map <- history_map_mouse_entrez
	} else {
		history_map <- history_map_human_entrez
	}
	
	entrez_ids_modern <- history_map$new[match(x = entrez_ids, table = history_map$old)]
	entrez_ids_modern[is.na(entrez_ids_modern)] <- entrez_ids[is.na(entrez_ids_modern)]
	entrez_ids_modern <- unname(entrez_ids_modern)
	
	if (any(duplicated(entrez_ids_modern))) {
		warning("Modernizing Entrez IDs resulted in duplicated entries")
	}
	
	if (mode(entrez_ids) == "numeric") as.numeric(entrez_ids_modern)
	
	entrez_ids_modern
}
