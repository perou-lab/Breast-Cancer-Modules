check_id_type <- function(x, cutoff = 0.8) {
	if (sum(grepl("^[0-9]+$", x)) / length(x) > cutoff) {
		return("Entrez")
	} else if (sum(grepl("^ENS(MUS)?G[0-9]+(\\.[0-9]+)?$", x)) / length(x) > cutoff) {
		return("Ensembl")
	} else if (sum(grepl("^[A-Za-z0-9\\-\\.]+$", x)) / length(x) > cutoff) {
		return("Symbols")
	} else {
		return(NA)
	}
}
