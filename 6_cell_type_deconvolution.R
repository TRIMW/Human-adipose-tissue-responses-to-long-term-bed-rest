
# Cell type deconvolution of bulk adipose RNA-seq
#
# Why standard scRNA-seq is inadequate as a reference for adipose:
#   Mature adipocytes are ~50-200 µm in diameter and filled with a lipid droplet.
#   They lyse during single-cell dissociation, so scRNA-seq datasets are essentially
#   devoid of adipocyte transcriptomes. This makes CIBERSORT with LM22 or any
#   standard scRNA-seq reference highly unsuitable for adipose tissue.
#
# Solution: use SINGLE-NUCLEUS RNA-seq (snRNA-seq) as the reference.
#   Nuclei survive the isolation protocol regardless of lipid content.
#   Emont et al. 2022 (Nature, DOI:10.1038/s41586-022-04518-2, GEO: GSE176171)
#   profiled human subcutaneous adipose by snRNA-seq and captured mature adipocytes.
#
# This script implements TWO approaches:
#
#   APPROACH A (immediately runnable): ssGSEA marker-gene scoring via GSVA.
#     Uses curated adipose cell-type marker gene sets to score each sample.
#     Gives relative activity scores per cell type per sample, not absolute fractions.
#     Suitable for: pre/post comparison, correlation with phenotypes.
#
#   APPROACH B (requires download): Reference-based deconvolution with MuSiC.
#     Needs the Emont 2022 processed SingleCellExperiment (see instructions below).
#     Gives estimated cell type PROPORTIONS per sample.
#
# For deconvolution purposes, rlog counts are used (ranked-based ssGSEA is
# robust to the normalization choice; MuSiC expects raw/normalized counts).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
if (!exists("counts")) source("analysis_preamble.R")

OUT_DIR <- "/Users/willtrim/Documents/projs/bedrest/outputs/6_cell_type_deconvolution"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =====================================================================
# APPROACH A: ssGSEA with curated adipose cell-type marker genes
# =====================================================================

# Marker genes from Emont et al. 2022 and prior adipose biology literature.
# Using gene symbols; converted to ENSEMBL IDs below.
cell_type_markers_sym <- list(
  Adipocytes      = c("ADIPOQ", "FABP4", "PLIN1", "LPL", "PPARG", "CEBPA",
                      "FASN", "SCD", "DGAT1", "LEP", "CFD", "PLIN4"),
  Preadipocytes   = c("PDGFRA", "DPP4", "CD34", "THY1", "SEMA3C", "ICAM1",
                      "ANPEP", "CD55"),
  Macrophages     = c("CD68", "CSF1R", "MRC1", "CD163", "ITGAM", "FOLR2",
                      "LYVE1", "S100A8", "S100A9"),
  Endothelial     = c("PECAM1", "VWF", "CLDN5", "CDH5", "KDR", "FLT1", "PLVAP"),
  SmoothMuscle    = c("ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5", "PDGFRB"),
  T_cells         = c("CD3D", "CD3E", "CD8A", "CD4", "TRAC", "IL7R"),
  B_cells         = c("CD19", "MS4A1", "CD79A", "IGHM"),
  NK_cells        = c("GNLY", "NKG7", "KLRD1", "NCAM1"),
  Mast_cells      = c("TPSAB1", "TPSB2", "CPA3", "KIT")
)

# Convert to ENSEMBL IDs using sym_map from preamble
cell_type_markers_ens <- lapply(cell_type_markers_sym, function(syms) {
  ens <- sym_map$ensembl_gene_id[sym_map$hgnc_symbol %in% syms]
  ens[ens %in% rownames(counts)]
})
cell_type_markers_ens <- Filter(function(x) length(x) >= 3, cell_type_markers_ens)
message(sprintf("Cell types with ≥3 markers found in data: %s",
                paste(names(cell_type_markers_ens), collapse=", ")))

# Run ssGSEA (requires GSVA package from Bioconductor)
if (!requireNamespace("GSVA", quietly = TRUE))
  stop("Install GSVA: BiocManager::install('GSVA')")
library(GSVA)

expr_mat <- as.matrix(counts)

# Handle both old API (GSVA < 1.48) and new API (GSVA >= 1.48)
gsva_scores <- tryCatch({
  param <- ssgseaParam(expr_mat, cell_type_markers_ens)  # new API
  gsva(param, verbose = FALSE)
}, error = function(e) {
  gsva(expr_mat, cell_type_markers_ens, method = "ssgsea", verbose = FALSE)  # old API
})

# gsva_scores: cell_types × samples matrix
gsva_df             <- as.data.frame(t(gsva_scores))
gsva_df$sample      <- rownames(gsva_df)
gsva_df$participant <- sub("_(Pre|Post)$", "", gsva_df$sample)
gsva_df$timepoint   <- ifelse(grepl("_Post$", gsva_df$sample), "Post", "Pre")

# ---- Pre vs Post comparison for each cell type ----
cell_types <- names(cell_type_markers_ens)
prepost_res <- do.call(rbind, lapply(cell_types, function(ct) {
  pre_scores  <- gsva_df[gsva_df$timepoint == "Pre",  ct]
  post_scores <- gsva_df[gsva_df$timepoint == "Post", ct]
  # Pair by participant
  pre_ids  <- gsva_df[gsva_df$timepoint == "Pre",  "participant"]
  post_ids <- gsva_df[gsva_df$timepoint == "Post", "participant"]
  common   <- intersect(pre_ids, post_ids)
  pre_v    <- setNames(pre_scores,  pre_ids)[common]
  post_v   <- setNames(post_scores, post_ids)[common]
  tt       <- wilcox.test(post_v, pre_v, paired = TRUE, exact = FALSE)
  data.frame(
    cell_type  = ct,
    median_pre  = median(pre_v,  na.rm = TRUE),
    median_post = median(post_v, na.rm = TRUE),
    delta_median = median(post_v - pre_v, na.rm = TRUE),
    p_value    = tt$p.value,
    stringsAsFactors = FALSE
  )
}))
prepost_res$fdr <- p.adjust(prepost_res$p_value, "BH")
prepost_res      <- prepost_res[order(prepost_res$p_value), ]

# ---- Correlate post scores with arterial insulin ----
pheno_post <- pheno[pheno$Time.Point == "Post Bed Rest", ]
rownames(pheno_post) <- pheno_post$Participant.ID

post_df <- gsva_df[gsva_df$timepoint == "Post", ]
rownames(post_df) <- post_df$participant
shared  <- intersect(rownames(post_df), rownames(pheno_post))
insulin <- as.numeric(pheno_post[shared, "Art..Insulin..mIU.L."])

ct_insulin_res <- do.call(rbind, lapply(cell_types, function(ct) {
  x   <- post_df[shared, ct]
  ok  <- !is.na(x) & !is.na(insulin)
  res <- cor.test(x[ok], insulin[ok], method = "spearman", exact = FALSE)
  data.frame(cell_type = ct, r = res$estimate, p_value = res$p.value)
}))
ct_insulin_res$fdr <- p.adjust(ct_insulin_res$p_value, "BH")

# ---- Save ----
fwrite(as.data.table(gsva_df),        file.path(OUT_DIR, "A_ssgsea_scores.csv"))
fwrite(as.data.table(prepost_res),    file.path(OUT_DIR, "A_prepost_comparison.csv"))
fwrite(as.data.table(ct_insulin_res), file.path(OUT_DIR, "A_celltype_vs_insulin_post.csv"))

# =====================================================================
# APPROACH B: MuSiC reference-based deconvolution (Emont 2022)
# =====================================================================
# To use this section:
#
#   1. Download GSE176171 from GEO. The processed SingleCellExperiment RDS
#      file is available from the Seurat-processed version at:
#      https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171
#      Look for the .rds or .h5ad file in the supplementary data.
#
#   2. Install MuSiC: devtools::install_github("xuranw/MuSiC")
#
#   3. Uncomment and run the block below.
#
# MuSiC works well here because it:
#   - Uses a weighted least squares approach that handles cell-type-specific
#     expression variance in the reference
#   - Can incorporate snRNA-seq references (where adipocytes ARE captured)
#   - Provides interpretable cell type proportion estimates

# library(MuSiC)
# library(SingleCellExperiment)
#
# # Load Emont 2022 snRNA-seq reference (adjust path)
# emont_sce <- readRDS("/path/to/emont2022_subcutaneous_sce.rds")
#
# # Raw counts needed for MuSiC (not rlog)
# raw_counts_mat <- as.matrix(raw_counts)
#
# # Run MuSiC deconvolution
# music_res <- music_prop(
#   bulk.mtx   = raw_counts_mat,
#   sc.sce     = emont_sce,
#   clusters   = "cell_type",   # column in colData(emont_sce) with cell type labels
#   samples    = "sample_id",   # column in colData with donor/sample labels
#   select.ct  = NULL,          # NULL = use all cell types
#   verbose    = TRUE
# )
#
# proportions <- as.data.frame(music_res$Est.prop.weighted)
# fwrite(as.data.table(proportions, keep.rownames="sample"),
#        file.path(OUT_DIR, "B_music_proportions.csv"))

message(sprintf("Approach A done. Output: %s", OUT_DIR))
message("Approach B: uncomment the MuSiC block after downloading Emont 2022 reference.")
