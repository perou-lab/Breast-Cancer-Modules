#' Calculate Centroid-Based Modules
#' 
#' This function takes in a data matrix of Entrez IDs x samples and returns module scores.
#' 
#' @param x a matrix, of Ensembl IDs (p) x samples (n), already log-transformed
#' @param input_species character string, "human" or "mouse" (default = "human")
#' 
#' @returns a matrix of module scores (p) x samples (n)
#' 
#' @export
calculate_centroid_scores <- function(x, 
									  input_species = "human") {
	
	# Vectorized row standardization
	x <- t(scale(t(x)))
	x[is.na(x)] <- 0
	
	centroid_names <- names(centroids)
	all_scores <- sapply(X = 1:length(centroids), FUN = function(i_centroid) {
		centroid_i <- centroids[[i_centroid]][[input_species]]
		genes_keep <- intersect(row.names(x), row.names(centroid_i))
		
		if (length(genes_keep) == 0) {
			return(rep(NA, ncol(x)))
		}
		
		# Extract matching genes from both matrices
		x_subset <- x[genes_keep, , drop = FALSE]
		centroid_vec <- centroid_i[genes_keep, 1]
		
		if (grepl(pattern = "Euclidean", x = centroid_names[i_centroid], ignore.case = TRUE)) {
			# Vectorized Euclidean Distance: sqrt(sum((x - y)^2)) for each column
			# Handle NAs in centroid by using complete pairs only (matching cor behavior)
			centroid_na <- is.na(centroid_vec)
			if (any(centroid_na)) {
				# For each sample, compute distance using only complete pairs
				correlation_scores <- sapply(1:ncol(x_subset), function(j) {
					complete <- !centroid_na & !is.na(x_subset[, j])
					if (sum(complete) == 0) return(NA)
					sqrt(sum((x_subset[complete, j] - centroid_vec[complete])^2))
				})
			} else {
				diff_sq <- (x_subset - centroid_vec)^2
				correlation_scores <- sqrt(colSums(diff_sq, na.rm = TRUE))
			}
		} else if (grepl(pattern = "Pcorr", x = centroid_names[i_centroid], ignore.case = TRUE)) {
			# Pearson correlation - must match cor(..., use = "complete.obs") exactly
			# The original uses cor() which centers each pair using only complete observations
			# We need to match this behavior exactly to avoid numerical differences
			correlation_scores <- apply(x_subset, 2, function(sample_i) {
				cor(x = sample_i, y = centroid_vec, use = "complete.obs", method = "pearson")
			})
		} else if (grepl(pattern = "Scorr", x = centroid_names[i_centroid], ignore.case = TRUE)) {
			# Spearman correlation - must match cor(..., use = "complete.obs") exactly
			# The original uses cor() which handles ranking and missing values in a specific way
			correlation_scores <- apply(x_subset, 2, function(sample_i) {
				cor(x = sample_i, y = centroid_vec, use = "complete.obs", method = "spearman")
			})
		} else {
			stop("unrecognized centroid name")
		}
		return(correlation_scores)
	}) |> t()
	
	rownames(all_scores) <- centroid_names
	
	return(all_scores)
}
