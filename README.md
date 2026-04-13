# Breast-Cancer-Modules

Calculates gene expression module scores, including a curated collection of approximately 1000 biologically-relevant modules compiled from the literature, which can be used in differential expression analysis and other genomic analyses. The package supports median-based calculations, centroid-based correlations, and more complex models.

## Repository Contents
1. **`data/`**: Contains previously-derived data frames for converting between Entrez IDs, Ensembl IDs, and Gene Symbols. Also contains pre-compiled objects for median-based modules and centroid modules, generated using the `refresh_modules()` function.
2. **`man/`**: Contains function descriptions accessible by e.g. `?calculate_modules`
3. **`modules/`**: Contains the source data for median-based modules, single genes, centroid modules, and special models.
4. **`R/`**: Contains the functions and scripts to be used.

### Important functions
1. **`R/calculate_modules.R`**: Primary function to calculate module scores. This function will call `calculate_median_modules.R`, `calculate_single_genes.R`, `calculate_centroid_scores.R`, and `calculate_special_models.R`, depending on what is requested in the options for `calculate_modules.R`.
2. **`R/refresh_modules.R`**: Will re-generate the `data/*.rda` files in case modules need to be added or removed. Note that the `calculate_modules` function also accepts custom modules. See documentation for details.
3. **`R/modernize_entrez_ids.R`**: Some modules are from older papers, and Entrez IDs have been discontinued or replaced with newer IDs. This function will update Entrez IDs to more recent versions using NCBI's gene_history file.
4. **`R/special_model_*`**: Special model functions are identified by this prefix. This organization is not ideal but allows for flexibility in adding/removing special model functions.

### Installation/Use
For now, we recommend cloning the GitHub repository, then using the `devtools` package and its `load_all()` function to load this R package.

### Citation
The code for this R package was formalized for the analyses in Sutcliffe et. al. "Transcriptomic profiling of mouse mammary tumors enables prognostic and predictive biomarker discovery for human breast cancers", which is currently available on bioRxiv:

https://doi.org/10.64898/2026.02.28.707759
