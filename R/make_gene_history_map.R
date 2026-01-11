# wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz
# gunzip gene_history.gz
make_gene_history_map <- function() {
	gene_history <- read.delim(file = "~/lab/gene_history", check.names = FALSE)
	
	tax_id_human <- "9606"
	tax_id_mouse <- "10090"
	
	gene_history_human <- gene_history[which(gene_history$`#tax_id` == tax_id_human), -1]
	gene_history_mouse <- gene_history[which(gene_history$`#tax_id` == tax_id_mouse), -1]
	
	make_update_map <- function(old_names, new_names) {
		
		old_names <- old_names[!is.na(old_names)]
		new_names <- new_names[!is.na(new_names)]
		
		old_names <- as.character(old_names)
		new_names <- as.character(new_names)
		
		# Create lookup hash (named vector is fastest in R)
		lookup <- setNames(new_names, old_names)
		
		# Initialize: each ID points to itself initially
		terminals <- setNames(unique(old_names), unique(old_names))
		
		# Iteratively resolve chains
		max_iter <- 100  # Should be plenty for short chains
		
		for (iter in 1:max_iter) {
			# Vectorized lookup: find next step for all current terminals
			old_terminals <- terminals
			next_steps <- lookup[terminals]
			
			# Update terminals where a next step exists (not NA)
			has_next <- !is.na(next_steps)
			terminals[has_next] <- next_steps[has_next]
			
			# Check convergence - if nothing changed, we're done
			if (identical(old_terminals, terminals)) {
				break
			}
		}
		
		# Create result - only return mappings where from != terminal
		result_mask <- names(terminals) != terminals
		
		if (any(result_mask)) {
			result <- data.frame(
				old = names(terminals)[result_mask],
				new = terminals[result_mask],
				stringsAsFactors = FALSE
			)
		} else {
			result <- data.frame(old = character(0), new = character(0))
		}
		
		result <- result[result$new != "-", ]
		
		rownames(result) <- NULL
		
		return(result)
	}
	
	history_map_human_entrez <- make_update_map(old_names = gene_history_human$Discontinued_GeneID, 
												new_names = gene_history_human$GeneID)
	history_map_human_symbol <- make_update_map(old_names = gene_history_human$Discontinued_Symbol, 
												new_names = gene_history_human$GeneID)
	
	history_map_mouse_entrez <- make_update_map(old_names = gene_history_mouse$Discontinued_GeneID, 
												new_names = gene_history_mouse$GeneID)
	history_map_mouse_symbol <- make_update_map(old_names = gene_history_mouse$Discontinued_Symbol, 
												new_names = gene_history_mouse$GeneID)
	
	save(history_map_human_entrez, file = "./data/history_map_human_entrez.rda")
	save(history_map_human_symbol, file = "./data/history_map_human_symbol.rda")
	save(history_map_mouse_entrez, file = "./data/history_map_mouse_entrez.rda")
	save(history_map_mouse_symbol, file = "./data/history_map_mouse_symbol.rda")
	
}