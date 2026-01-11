convert_to_human_entrez <- function(x, 
									from_id_type, 
									input_species = "human") {
	
	if (!from_id_type %in% c("ensembl", "symbol", "entrez")) {
		stop("from_id_type must be one of: 'ensembl', 'symbol', 'entrez'")
	}
	
	if (!input_species %in% c("human", "mouse")) {
		stop("input_species must be either 'human' or 'mouse'")
	}
	
	if (is.null(rownames(x))) {
		stop("x must have gene identifiers as row names")
	}
	
	input_ids <- row.names(x)
	
	if (input_species == "mouse") {
		if (from_id_type != "entrez") {
			
			if (from_id_type == "ensembl") {
				input_ids <- gsub(pattern = "\\..*$", replacement = "", x = input_ids)
				input_ids <- mouse_id_table$entrez_id[match(x = input_ids, table = mouse_id_table$ensembl_id)]
			} else if (from_id_type == "symbol") {
				input_ids <- mouse_id_table$gene_symbol[match(x = input_ids, table = mouse_id_table$ensembl_id)]
			}
			
		}
		
		input_ids <- orthologs$human_entrez_id[match(x = input_ids, table = orthologs$mouse_entrez_id)]

	} else {
		
		if (from_id_type == "ensembl") {
			input_ids <- gsub(pattern = "\\..*$", replacement = "", x = input_ids)
			input_ids <- human_id_table$entrez_id[match(x = input_ids, table = human_id_table$ensembl_id)]
		} else if (from_id_type == "symbol") {
			input_ids <- human_id_table$entrez_id[match(x = input_ids, table = human_id_table$gene_symbol)]
		}
		
	}
	
	collapse_rows <- function(x, new_ids) {
		# x: numeric matrix (genes × samples)
		# new_ids: vector of Entrez IDs (same length as rownames(x))
		
		if (length(new_ids) != nrow(x)) {
			stop("new_ids must be the same length as nrow(x)")
		}
		
		# Collapse rows by summing over duplicate IDs
		collapsed <- rowsum(x, group = new_ids)
		
		return(collapsed)
	}
	
	x <- x[!is.na(input_ids), ]
	input_ids <- input_ids[!is.na(input_ids)]
	x <- collapse_rows(x, input_ids)
	
	return(x)
}
