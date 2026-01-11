# IMS_Score_CCR.2018_PMID.29921729
ims_score <- function(x, input_species = "human") {
	ims_score_model_i <- ims_score_model
	if (input_species == "mouse") {
		ims_score_model_i$ensembl_id <- babel_ensembl$mouse_ensembl_id[match(x = ims_score_model_i$ensembl_id, table = babel_ensembl$human_ensembl_id)]
	}
	ims_score_model_i <- ims_score_model_i[ims_score_model_i$ensembl_id %in% row.names(x) & !is.na(ims_score_model_i$ensembl_id), ]
	
	x <- x[ims_score_model_i$ensembl_id, ]
	
	# x <- t(scale(t(x)))
	# x[is.na(x)] <- 0
	
	x_ims_score <- apply(X = x, MARGIN = 2, FUN = function(x_i) sum(x_i * ims_score_model_i$Coefficient))
	
	return(list("IMS_Score_CCR.2018_PMID.29921729" = x_ims_score))
}