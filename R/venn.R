#' Quick overlap checking function
#' 
#' @export
venn <- function(x, y) {
	x_and_y <- intersect(x, y)
	x_only <- setdiff(x, y)
	y_only <- setdiff(y, x)
	
	res <- c("1 only" = length(x_only),
			 "intersection" = length(x_and_y),
			 "2 only" = length(y_only))
	
	return(res)
}