#' Collapse duplicates row
#' 
#' This function takes in matrix and collapses the rows by the `ids`. This is used for when e.g. multiple probes map to the same gene, and those values (rows) need to be aggregated.
#' 
#' @param x, a matrix of features x samples
#' @param ids, a vector of ids (same length as features) to collapse by
#' @param sum_or_mean, a string, method to collapse rows (default = "sum")
#' 
#' @returns a matrix of collapsed features x samples
#' 
#' @export
collapse_duplicates <- function(x, ids = row.names(x), sum_or_mean = "sum") {
	if (!(sum_or_mean %in% c("sum", "mean"))) {
		stop("`sum_or_mean` must be either string `sum` or `mean`")
	}
	
	if (nrow(x) != length(ids)) {
		stop("number of rows of `x` must match length of `ids`")
	}
	
	x <- x[!is.na(ids), ]
	ids <- ids[!is.na(ids)]
	
	if (sum_or_mean == "sum") {
		x <- rowsum(x, group = ids)
	} else {
		x_temp <- aggregate(x, list(ids), mean)
		row.names(x_temp) <- x_temp[, 1]
		x_temp <- x_temp[-1]
		x_temp <- as.matrix(x_temp)
		colnames(x_temp) <- colnames(x)
		x <- x_temp
	}
	
	x <- x[order(row.names(x)), ]
	
	return(x)
}
