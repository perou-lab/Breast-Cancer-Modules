# Proliferation_11_Mean_JCO.2009_PMID.19204204
proliferation_11_mean <- function(x, input_species = "human") {
	prolif_genes <- c("BIRC5", "CCNB1", "CDC20", "CEP55", "MKI67", "NDC80", "NUF2", "PTTG1", "RRM2", "TYMS", "UBE2C")
	
	human_ensembl <- gene_map_human$ensembl_id[match(x = prolif_genes, table = gene_map_human$gene_name)]
	mouse_ensembl <- babel_ensembl$mouse_ensembl_id[match(x = human_ensembl, table = babel_ensembl$human_ensembl_id)]
	
	ensembl_ids <- ifelse(input_species == "mouse", list(mouse_ensembl), list(human_ensembl))[[1]]
	
	prolif_genes_keep <- intersect(ensembl_ids, row.names(x))
	
	return(list("Proliferation_11_Mean_JCO.2009_PMID.19204204" = colMeans(x = x[prolif_genes_keep, ], na.rm = TRUE)))
}