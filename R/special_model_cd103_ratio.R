# CD103_Ratio_Cancer.Cell.2014_PMID.25446897
cd103_ratio <- function(x, input_species = "human") {
	cd103_pos_symbols <- c("KIT", "CCR7", "BATF3", "FLT3", "ZBTB46", "IRF8", "BTLA", "MYCL", "CLEC9A")
	cd103_neg_symbols <- c("C5AR1", "LYVE1", "ABCC3", "MRC1", "SIGLEC1", "STAB1", "C1QB", "C1QA", "TMEM37", "MERTK", "C1QC", "TMEM119", "MS4A7", "APOE", "CYP4F2", "TREM2", "TLR7")

	cd103_pos_genes <- gene_map_human$ensembl_id[match(x = cd103_pos_symbols, table = gene_map_human$gene_name)] |> as.character()
	cd103_neg_genes <- gene_map_human$ensembl_id[match(x = cd103_neg_symbols, table = gene_map_human$gene_name)] |> as.character()

	if (input_species == "mouse") {
		cd103_pos_genes <- babel_ensembl$mouse_ensembl_id[match(x = cd103_pos_genes, table = babel_ensembl$human_ensembl_id)]
		cd103_neg_genes <- babel_ensembl$mouse_ensembl_id[match(x = cd103_neg_genes, table = babel_ensembl$human_ensembl_id)]
	}

	cd103_pos_genes <- cd103_pos_genes[cd103_pos_genes %in% row.names(x)]
	cd103_neg_genes <- cd103_neg_genes[cd103_neg_genes %in% row.names(x)]

	if (length(cd103_pos_genes) == 0 || length(cd103_neg_genes) == 0) {
		return(list("CD103_Ratio_Cancer.Cell.2014_PMID.25446897" = rep(0, ncol(x))))
	}

	cd103_ratio <- colMeans(x = x[cd103_pos_genes, ], na.rm = TRUE) / colMeans(x = x[cd103_neg_genes, ], na.rm = TRUE)

	cd103_ratio[is.na(cd103_ratio)] <- 0

	return(list("CD103_Ratio_Cancer.Cell.2014_PMID.25446897" = cd103_ratio))
}
