validate_entrez <- function(entrez_ids) {
	input_type <- mode(entrez_ids)
	if (!input_type %in% c("numeric", "character")) {
		stop("Input is not a numeric or character vector")
	}
	
	if (any(is.na(suppressWarnings(as.numeric(entrez_ids))))) {
		if (length(entrez_ids) == 1) {
			stop("Input is not a number")
		} else {
			ind_not_entrez <- which(is.na(suppressWarnings(as.numeric(entrez_ids))))[1]
			stop(paste("Input contains non-numbers, e.g. at index:", ind_not_entrez))
		}
	}
	
	if (any(duplicated(entrez_ids))) {
		warning("Input contains duplicate Entrez IDs")
	}
	
}
