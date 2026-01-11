# Macrophage_Polarity_CXCL9.SPP1.Ratio_Science.2023_PMID.37535729
macrophage_polarity <- function(x, input_species = "human") {
	#      entrezgene_id ensembl_gene_id hgnc_symbol
	# 2949          4283 ENSG00000138755       CXCL9
	#      entrezgene_id ensembl_gene_id hgnc_symbol
	# 4626          6696 ENSG00000118785        SPP1
	
	cxcl9 <- gene_map_human$ensembl_id[match(x = "CXCL9", table = gene_map_human$gene_name)] |> as.character()
	spp1 <- gene_map_human$ensembl_id[match(x = "SPP1", table = gene_map_human$gene_name)] |> as.character()
	
	if (input_species == "mouse") {
		cxcl9 <- babel_ensembl$mouse_ensembl_id[match(x = cxcl9, table = babel_ensembl$human_ensembl_id)]
		spp1 <- babel_ensembl$mouse_ensembl_id[match(x = spp1, table = babel_ensembl$human_ensembl_id)]
	}
	
	if (!all(c(cxcl9, spp1) %in% row.names(x)) & !is.na(cxcl9) & !is.na(spp1)) {
		return(rep(0, ncol(x)))
	}
	
	ratio <- unlist(x[cxcl9, ]) / unlist(x[spp1, ])
	
	return(list("Macrophage_Polarity_CXCL9.SPP1.Ratio_Science.2023_PMID.37535729" = ratio))
}