# GHI_RS_Model_NEJM.2004_PMID.15591335
oncotypeDX_research <- function(x, input_species = "human") {
	# Because we take in upper quartile-normalized counts, we will ignore Oncotype's normalization (redundant + only relies on 5 genes)
	
	genes <- c("GRB7", "ERBB2", "ESR1", "PGR", "BCL2", "SCUBE2", "BIRC5", "MKI67", "MYBL2", "CCNB1", "AURKA", "CTSV", "MMP11", "CD68", "GSTM1", "BAG1")
	
	human_ensembl <- gene_map_human$ensembl_id[match(x = genes, table = gene_map_human$gene_name)]
	mouse_ensembl <- babel_ensembl$mouse_ensembl_id[match(x = human_ensembl, table = babel_ensembl$human_ensembl_id)]
	
	ensembl_ids <- ifelse(input_species == "mouse", list(mouse_ensembl), list(human_ensembl))[[1]]
	
	genes <- setNames(ensembl_ids, genes)
	
	if (!all(genes %in% row.names(x))) {
		x_add <- matrix(data = rep(0, sum(!(genes %in% row.names(x))) * ncol(x)), 
						ncol = ncol(x))
		row.names(x_add) <- genes[!(genes %in% row.names(x))]
		colnames(x_add) <- colnames(x)
		x <- rbind(x, x_add)
	}
	
	
	GRB7_group_score <-
		0.9 * x[genes["GRB7"], ] + # GRB7
		0.1 * x[genes["ERBB2"], ]   # ERBB2
	GRB7_group_score[GRB7_group_score < 8] <- 8
	
	ER_group_score <-
		0.8 * x[genes["ESR1"], ] + # ESR1
		1.2 * x[genes["PGR"], ] + # PGR
		1.0 * x[genes["BCL2"], ] +  # BCL2
		1.0 * x[genes["SCUBE2"], ]  # SCUBE2
	ER_group_score <- ER_group_score / 4
	
	proliferation_group_score <-
		x[genes["BIRC5"], ] +  # BIRC5
		x[genes["MKI67"], ] + # MKI67
		x[genes["MYBL2"], ] + # MYBL2
		x[genes["CCNB1"], ] +  # CCNB1
		x[genes["AURKA"], ]   # AURKA
	proliferation_group_score <- proliferation_group_score / 5
	proliferation_group_score[proliferation_group_score < 6.5] <- 6.5
	
	invasion_group_score <-
		x[genes["CTSV"], ] + # CTSV
		x[genes["MMP11"], ]   # MMP11
	invasion_group_score <- invasion_group_score / 2
	
	recurrence_score <-
		0.47 * GRB7_group_score +
		-0.34 * ER_group_score +
		1.04 * proliferation_group_score +
		0.10 * invasion_group_score +
		0.05 * x[genes["CD68"], ] +   # CD68
		-0.08 * x[genes["GSTM1"], ] + # GSTM1
		-0.07 * x[genes["BAG1"], ]    # BAG1
	
	return(list("GHI_RS_Model_NEJM.2004_PMID.15591335" = recurrence_score))
}