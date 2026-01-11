#' Calculate Special Model Modules
#' 
#' This function takes in a data matrix of Entrez IDs x samples and returns module scores.
#' 
#' @param x a matrix, of Entrez IDs (p) x samples (n), already log-transformed
#' @param special_model_functions vector of files, with each file containing the source function for one special model
#' @param input_species character string, "human" or "mouse" (default = "human")
#' 
#' @returns a matrix of module scores (p) x samples (n)
#' 
#' @export
calculate_special_models <- function(x,
									 # special_model_functions = list.files("./R/special_model_functions/", full.names = TRUE),
									 input_species = "human") {
	
	special_model_functions
	
	x_special_models <- sapply(X = special_model_functions, USE.NAMES = FALSE, FUN = function(f) {
		f(x, input_species)
	})
	
	do.call(what = rbind, x_special_models)
}

# special_model_functions <- list(
# 	chemo_endocrine_score,
# 	luma_her2_score, 
# 	MREscore, 
# 	oncotypeDX_research, 
# 	proliferation_11_mean, 
# 	replication_stress, 
# 	ror_s_score, 
# 	ror_p_score, 
# 	cd103_ratio, 
# 	macrophage_polarity, 
# 	differentiation_score, 
# 	ims_score, 
# 	rss_score
# )
