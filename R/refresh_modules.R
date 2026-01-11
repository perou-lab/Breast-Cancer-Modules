#' Refresh module database
#' 
#' This function regenerates the `median_modules.RDS` and `all_centroids.RDS` files to be loaded in their respective functions
#' 
#' @export
refresh_modules <- function() {
	# Published modules
	files_median <- list.files("./modules/median_modules", recursive = TRUE, full.names = TRUE)
	median_modules <- lapply(X = files_median, FUN = function(file_i) {
		x_i <- scan(file = file_i, what = character(), quiet = TRUE)
		
		if (check_id_type(x_i) == "Entrez") {
			# Update potentially outdated genes
			x_i[x_i %in% names(gene_history_map)] <- gene_history_map[x_i[x_i %in% names(gene_history_map)]]
			y_i_human <- gene_map_human[which(gene_map_human$entrez_id %in% x_i), ]
			mouse_orthologs <- babel_ensembl$mouse_ensembl_id[which(babel_ensembl$human_ensembl_id %in% y_i_human$ensembl_id)]
			y_i_mouse <- gene_map_mouse[which(gene_map_mouse$ensembl_id %in% mouse_orthologs), ]
		} else if (check_id_type(x_i) == "Ensembl") {
			y_i_human <- gene_map_human[which(gene_map_human$ensembl_id %in% x_i), ]
			mouse_orthologs <- babel_ensembl$mouse_ensembl_id[which(babel_ensembl$human_ensembl_id %in% y_i_human$ensembl_id)]
			y_i_mouse <- gene_map_mouse[which(gene_map_mouse$ensembl_id %in% mouse_orthologs), ]
		} else if (check_id_type(x_i) == "Symbols") {
			y_i_human <- gene_map_human[which(gene_map_human$gene_name %in% x_i), ]
			mouse_orthologs <- babel_ensembl$mouse_ensembl_id[which(babel_ensembl$human_ensembl_id %in% y_i_human$ensembl_id)]
			y_i_mouse <- gene_map_mouse[which(gene_map_mouse$ensembl_id %in% mouse_orthologs), ]
		}
		
		return(list("human" = y_i_human,
					"mouse" = y_i_mouse))
	})
	names(median_modules) <- basename(files_median)
	
	# TODO make package update?
	save(object = median_modules, file = "./data/median_modules.rda")
	
	message("Updated median module database.")
	
	files_single_genes <- list.files("./modules/single_genes")
	y_i_human <- gene_map_human[match(x = files_single_genes, table = gene_map_human$gene_name), ]
	y_i_human$module_name <- paste(y_i_human$gene_name, "Single_Gene", sep = "_")
	
	y_i_mouse <- y_i_human
	y_i_mouse$mouse_ensembl_id <- babel_ensembl$mouse_ensembl_id[match(x = y_i_mouse$ensembl_id, table = babel_ensembl$human_ensembl_id)]
	y_i_mouse$mouse_gene_name <- gene_map_mouse$gene_name[match(x = y_i_mouse$mouse_ensembl_id, table = gene_map_mouse$ensembl_id)]
	y_i_mouse$mouse_entrez_id <- gene_map_mouse$entrez_id[match(x = y_i_mouse$mouse_ensembl_id, table = gene_map_mouse$ensembl_id)]
	
	y_i_mouse$ensembl_id <- y_i_mouse$mouse_ensembl_id
	y_i_mouse$entrez_id <- y_i_mouse$mouse_entrez_id
	y_i_mouse$gene_name <- y_i_mouse$mouse_gene_norm
	
	single_genes <- list("human" = y_i_human,
						 "mouse" = y_i_mouse)
	
	save(object = single_genes, file = "./data/single_genes.rda")
	
	message("Updated single gene database.")
	
	files_centroids <- list.files(path = "./modules/centroids", full.names = TRUE, include.dirs = FALSE) |> grep(pattern = "(Pcorr|Scorr|Euclidean)", value = TRUE)
	centroids <- lapply(X = files_centroids, FUN = function(centroid_i) {
		x_i <- read.delim(file = centroid_i, header = TRUE, check.names = FALSE)
		x_i_genes <- x_i[, 1]
		id_type <- check_id_type(x_i_genes)
		
		if (id_type == "Entrez") {
			# Update potentially outdated genes
			x_i_genes[x_i_genes %in% names(gene_history_map)] <- gene_history_map[x_i_genes[x_i_genes %in% names(gene_history_map)]]
			y_i_human <- gene_map_human[match(x = x_i_genes, table = gene_map_human$entrez_id), ]
			mouse_orthologs <- babel_ensembl$mouse_ensembl_id[match(x = y_i_human$ensembl_id, table = babel_ensembl$human_ensembl_id)]
			y_i_mouse <- gene_map_mouse[match(x = mouse_orthologs, table = gene_map_mouse$ensembl_id), ]
		} else if (id_type == "Ensembl") {
			y_i_human <- gene_map_human[match(x = x_i_genes, table = gene_map_human$ensembl_id), ]
			mouse_orthologs <- babel_ensembl$mouse_ensembl_id[match(x = y_i_human$ensembl_id, table = babel_ensembl$human_ensembl_id)]
			y_i_mouse <- gene_map_mouse[match(x = mouse_orthologs, table = gene_map_mouse$ensembl_id), ]
		} else if (id_type == "Symbols") {
			y_i_human <- gene_map_human[match(x = x_i_genes, table = gene_map_human$gene_name), ]
			mouse_orthologs <- babel_ensembl$mouse_ensembl_id[match(x = y_i_human$ensembl_id, table = babel_ensembl$human_ensembl_id)]
			y_i_mouse <- gene_map_mouse[match(x = mouse_orthologs, table = gene_map_mouse$ensembl_id), ]
		}
		
		
		centroid_i_human <- x_i[!is.na(y_i_human$ensembl_id), 2, drop = FALSE]
		centroid_i_human <- as.matrix(centroid_i_human)
		row.names(centroid_i_human) <- y_i_human$ensembl_id[!is.na(y_i_human$ensembl_id)]
		
		centroid_i_mouse <- x_i[!is.na(y_i_mouse$ensembl_id), 2, drop = FALSE]
		centroid_i_mouse <- as.matrix(centroid_i_mouse)
		row.names(centroid_i_mouse) <- y_i_mouse$ensembl_id[!is.na(y_i_mouse$ensembl_id)]
		
		# message(paste("replaced", sum(x_i %in% names(gene_history_map)), "of", length(x_i), "genes"))
		return(list("human" = centroid_i_human,
					"mouse" = centroid_i_mouse))

	})
	names(centroids) <- basename(files_centroids)
	
	save(object = centroids, file = "./data/centroids.rda")
	
	message("Updated centroid database.")
	
	
	replication_stress_model <- read.delim(file = "./modules/special_models/replication_stress_model.txt")
	save(replication_stress_model, file = "./data/replication_stress_model.rda")
	
	differentiation_centroid_MS <- read.delim(file = "./modules/special_models/differentiationCentroids_LimDWD_MS.txt")
	differentiation_centroid_mL <- read.delim(file = "./modules/special_models/differentiationCentroids_LimDWD_mL.txt")
	
	save(differentiation_centroid_MS, file = "./data/differentiation_centroid_MS.rda")
	save(differentiation_centroid_mL, file = "./data/differentiation_centroid_mL.rda")
	
	ims_score_model <- read.delim(file = "./modules/special_models/IMS_Score_Clin.Cancer.Res.2018_PMID.29921729.txt")
	rss_score_model <- read.delim(file = "./modules/special_models/RSS_Score_Clin.Cancer.Res.2018_PMID.29921729.txt")
	
	ims_score_model$ensembl_id <- gene_map_human$ensembl_id[match(x = ims_score_model$EntrezID, table = gene_map_human$entrez_id)]
	rss_score_model$ensembl_id <- gene_map_human$ensembl_id[match(x = rss_score_model$EntrezID, table = gene_map_human$entrez_id)]
	
	save(ims_score_model, file = "./data/ims_score_model.rda")
	save(rss_score_model, file = "./data/rss_score_model.rda")
	
	special_models <- list.files(path = "./R", full.names = TRUE, pattern = "^special_model_")
	special_model_functions <- list()
	for (i in special_models) {
		i <- gsub(pattern = "\\./R/special_model_", replacement = "", x = i)
		special_model_functions <- append(x = special_model_functions, get(tools::file_path_sans_ext(i)))
	}
	save(special_model_functions, file = "./data/special_model_functions.rda")
	
	message("Updated special model database.")
	
	message("Please re-load library for changes to take effect.")
}
