#' Convert Gene Symbols to Ensembl IDs
#' 
#' This function takes in a vector of Gene Symbols and return a vector of mapped Ensembl IDs.
#' 
#' @param gene_symbols a vector of characters
#' @param species character string, "human" or "mouse" (default = "human")
#' 
#' @returns a vector of mapped Ensembl IDs
#' 
#' @importFrom HGNChelper checkGeneSymbols
#' 
#' @export
symbol_to_ensembl <- function(gene_symbols, species = "human") {
	
	if (!species %in% c("human", "mouse")) {
		stop("`species` argument expects either \"human\" or \"mouse\"")
	}
	
	if (species == "mouse") {
		gene_map <- gene_map_mouse
		history_map <- history_map_mouse_symbol
	} else {
		gene_map <- gene_map_human
		history_map <- history_map_human_symbol
	}
	
	ensembl_ids <- gene_map$ensembl_id[match(x = gene_symbols, table = gene_map$gene_name)]
	
	unmapped_symbols_index <- which(is.na(ensembl_ids))
	if (length(unmapped_symbols_index) > 0) {
		unmapped_symbols <- gene_symbols[unmapped_symbols_index]
		
		# Check for discontinued gene symbols
		updated_entrez <- history_map_human_symbol$new[match(x = unmapped_symbols, history_map_human_symbol$old)]
		ensembl_ids_remapped <- entrez_to_ensembl(entrez_ids = updated_entrez[!is.na(updated_entrez)], species = species)
		ensembl_ids[unmapped_symbols_index[!is.na(updated_entrez)]] <- ensembl_ids_remapped
	}
	
	# Now optionally check for "invalid" gene symbols using HGNChelper (if installed)
	unmapped_symbols_index <- which(is.na(ensembl_ids))
	if (length(unmapped_symbols_index) > 0 & requireNamespace("HGNChelper", quietly = TRUE)) {
		unmapped_symbols <- gene_symbols[unmapped_symbols_index]
		hgnchelper_symbols <- HGNChelper::checkGeneSymbols(x = gene_symbols, species = species)$Suggested.Symbol |> suppressMessages()
		
		ensembl_ids[unmapped_symbols_index] <- gene_map$ensembl_id[match(x = hgnchelper_symbols, table = gene_map$gene_name)]
	}
	

	
	if (any(is.na(ensembl_ids))) {
		warning("Some Ensembl IDs missing (may be discontinued genes)")
	}
	
	if (any(duplicated(ensembl_ids))) {
		warning("Mapped Ensembl IDs contain duplicated entries")
	}
	
	ensembl_ids
}
