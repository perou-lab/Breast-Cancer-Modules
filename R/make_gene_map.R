make_gene_map <- function() {
	
	# HUMAN -----
	
	# wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz
	# gunzip gene2ensembl.gz
	# awk 'NR==1 || $1=="9606"' gene2ensembl > gene2ensembl_human
	
	ncbi_gene2ensembl_human <- read.delim(file = "~/lab/gene2ensembl_human", header = TRUE, sep = "\t", check.names = FALSE)
	ncbi_gene2ensembl_human <- ncbi_gene2ensembl_human[, c(2, 3)]
	colnames(ncbi_gene2ensembl_human) <- c("entrez_id", "ensembl_id")
	ncbi_gene2ensembl_human <- ncbi_gene2ensembl_human[order(ncbi_gene2ensembl_human$ensembl_id), ]
	ncbi_gene2ensembl_human <- unique(ncbi_gene2ensembl_human)
	rownames(ncbi_gene2ensembl_human) <- NULL
	
	# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz
	# gunzip gencode.v48.annotation.gtf.gz
	
	gencode48 <- rtracklayer::import("~/lab/gencode.v48.annotation.gtf")
	gencode48 <- as.data.frame(gencode48)
	gencode48$ensembl_id <- gsub(pattern = "\\..*$", replacement = "", x = gencode48$gene_id)
	gencode48 <- gencode48[, c("ensembl_id", "gene_name", "type", "seqnames", "strand", "start", "end", "gene_type")]
	gencode48$seqnames <- gsub(pattern = "chr", replacement = "", x = as.character(gencode48$seqnames)) |> factor(levels = c(1:22, "X", "Y", "M"))
	gencode48 <- gencode48[order(gencode48$seqnames, gencode48$start, gencode48$type, gencode48$gene_type != "protein_coding", gencode48$ensembl_id), ]
	gencode48 <- gencode48[!duplicated(gencode48[, c("ensembl_id", "gene_name", "seqnames")]), ]
	rownames(gencode48) <- NULL
	
	gene_map_human <- merge(x = ncbi_gene2ensembl_human, y = gencode48, by = "ensembl_id", all = TRUE)
	gene_map_human <- gene_map_human[order(gene_map_human$seqnames, gene_map_human$start, gene_map_human$type, gene_map_human$gene_type != "protein_coding", gene_map_human$ensembl_id), ]
	rownames(gene_map_human) <- NULL
	
	save(gene_map_human, file = "./data/gene_map_human.rda")
	
	
	# MOUSE -----
	
	# wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz
	# gunzip gene2ensembl.gz
	# awk 'NR==1 || $1=="10090"' gene2ensembl > gene2ensembl_mouse
	
	ncbi_gene2ensembl_mouse <- read.delim(file = "~/lab/gene2ensembl_mouse", header = TRUE, sep = "\t", check.names = FALSE)
	ncbi_gene2ensembl_mouse <- ncbi_gene2ensembl_mouse[, c(2, 3)]
	colnames(ncbi_gene2ensembl_mouse) <- c("entrez_id", "ensembl_id")
	ncbi_gene2ensembl_mouse <- unique(ncbi_gene2ensembl_mouse)
	ncbi_gene2ensembl_mouse <- ncbi_gene2ensembl_mouse[order(ncbi_gene2ensembl_mouse$ensembl_id), ]
	rownames(ncbi_gene2ensembl_mouse) <- NULL
	
	# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.annotation.gtf.gz
	# gunzip gencode.vM37.annotation.gtf.gz
	
	gencodeM37 <- rtracklayer::import("~/lab/gencode.vM37.annotation.gtf")
	gencodeM37 <- as.data.frame(gencodeM37)
	gencodeM37$ensembl_id <- gsub(pattern = "\\..*$", replacement = "", x = gencodeM37$gene_id)
	gencodeM37 <- gencodeM37[, c("ensembl_id", "gene_name", "type", "seqnames", "strand", "start", "end", "gene_type")]
	
	gencodeM37$seqnames <- gsub(pattern = "chr", replacement = "", x = as.character(gencodeM37$seqnames)) |> factor(levels = c(1:19, "X", "Y", "M"))
	gencodeM37 <- gencodeM37[order(gencodeM37$seqnames, gencodeM37$start, gencodeM37$type, gencodeM37$gene_type != "protein_coding", gencodeM37$ensembl_id), ]
	gencodeM37 <- gencodeM37[!duplicated(gencodeM37[, c("ensembl_id", "gene_name", "seqnames")]), ]
	rownames(gencodeM37) <- NULL
	
	gene_map_mouse <- merge(x = ncbi_gene2ensembl_mouse, y = gencodeM37, by = "ensembl_id", all = TRUE)
	gene_map_mouse <- gene_map_mouse[order(gene_map_mouse$seqnames, gene_map_mouse$start, gene_map_mouse$type, gene_map_mouse$gene_type != "protein_coding", gene_map_mouse$ensembl_id), ]
	rownames(gene_map_mouse) <- NULL
	
	save(gene_map_mouse, file = "./data/gene_map_mouse.rda")
	
	# -----
	# wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene_orthologs.gz
	# gunzip gene_orthologs.gz
	# awk 'NR==1 || ($1 == "9606" && $4 == "10090") || ($1 == "10090" && $4 == "9606")' gene_orthologs > human_mouse_orthologs
	# human_mouse_orthologs <- read.delim(file = "~/lab/human_mouse_orthologs", check.names = FALSE)
	# human_orthologs <- ifelse(human_mouse_orthologs$`#tax_id` == "9606", human_mouse_orthologs$GeneID, human_mouse_orthologs$Other_GeneID)
	# mouse_orthologs <- ifelse(human_mouse_orthologs$`#tax_id` == "10090", human_mouse_orthologs$GeneID, human_mouse_orthologs$Other_GeneID)
	
	babel_human <- babelgene::orthologs(genes = gene_map_human$ensembl_id, species = "mouse", human = TRUE, top = TRUE)
	babel_mouse <- babelgene::orthologs(genes = gene_map_mouse$ensembl_id, species = "mouse", human = FALSE, top = TRUE)
	
	babel_all <- rbind(babel_human, babel_mouse)
	babel_all <- unique(babel_all)
	babel_all <- babel_all[order(babel_all$support_n, decreasing = TRUE), ]
	babel_ensembl <- babel_all[!duplicated(babel_all$human_ensembl), c("human_ensembl", "ensembl")]
	colnames(babel_ensembl) <- c("human_ensembl_id", "mouse_ensembl_id")
	babel_ensembl <- unique(babel_ensembl)
	babel_ensembl <- babel_ensembl[order(babel_ensembl$human_ensembl_id), ]
	rownames(babel_ensembl) <- NULL
	
	save(babel_ensembl, file = "./data/babel_ensembl.rda")
	
}