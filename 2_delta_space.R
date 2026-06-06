
# Delta-space analysis
#
# Computes post - pre changes for each participant (gene expression and phenotypes).
# This leverages the paired within-subject design and controls for baseline differences.
#
# Four analyses:
#   A. TCA module eigengene (PC1 of delta TCA genes) vs all delta phenotypes
#   B. All gene deltas vs delta arterial insulin  --> genes that change with insulin change
#   C. All gene deltas vs TCA module eigengene    --> genes that change with TCA activity change
#   D. Intersection of B & C                      --> genes linked to BOTH insulin and TCA

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
if (!exists("counts")) source("analysis_preamble.R")

# Fast Spearman correlation + p-value for a matrix vs a single vector.
# mat: genes x samples (rows=genes); y: named vector length=samples
spearman_all <- function(mat, y) {
  ok  <- names(y)[!is.na(y)]
  y_  <- y[ok]
  m_  <- mat[, ok, drop = FALSE]
  n   <- length(ok)
  # rank transform
  ry  <- rank(y_)
  rm_ <- t(apply(m_, 1, rank))
  r   <- drop(cor(t(rm_), ry))  # vector of length nrow(mat)
  t_  <- r * sqrt(n - 2) / sqrt(pmax(1 - r^2, 1e-12))
  p   <- 2 * pt(-abs(t_), df = n - 2)
  data.frame(r = r, p_value = p, fdr = p.adjust(p, "BH"),
             row.names = rownames(mat))
}

# ---- Analysis A: TCA module eigengene vs all delta phenotypes ----
tca_delta_mat <- t(delta_expr[tca_genes, pheno_paired, drop = FALSE])
pca_tca       <- prcomp(tca_delta_mat, scale. = TRUE)
tca_eigen     <- pca_tca$x[, 1]
# Ensure positive loading = higher TCA activity post
if (cor(tca_eigen, rowMeans(tca_delta_mat), use = "complete.obs") < 0)
  tca_eigen <- -tca_eigen

pct_var <- round(100 * pca_tca$sdev[1]^2 / sum(pca_tca$sdev^2), 1)
message(sprintf("TCA module PC1 explains %.1f%% of variance across TCA genes", pct_var))

tca_pheno_res <- do.call(rbind, lapply(pheno_vars, function(v) {
  y   <- setNames(delta_pheno[, v], rownames(delta_pheno))
  ok  <- names(y)[!is.na(y)]
  res <- cor.test(tca_eigen[ok], y[ok], method = "spearman", exact = FALSE)
  data.frame(phenotype = v, r = res$estimate, p_value = res$p.value)
}))
tca_pheno_res$fdr <- p.adjust(tca_pheno_res$p_value, "BH")
tca_pheno_res     <- tca_pheno_res[order(tca_pheno_res$p_value), ]

# ---- Analysis B: All gene deltas vs delta arterial insulin ----
message("Analysis B: all genes vs delta insulin...")
delta_insulin <- setNames(delta_pheno[, "Art..Insulin..mIU.L."], rownames(delta_pheno))
b_res         <- spearman_all(delta_expr[, pheno_paired], delta_insulin)
b_res$gene    <- rownames(b_res)
b_res$is_tca  <- rownames(b_res) %in% tca_genes
b_res         <- b_res[order(b_res$fdr), ]

# ---- Analysis C: All gene deltas vs TCA module eigengene ----
message("Analysis C: all genes vs TCA eigengene...")
c_res      <- spearman_all(delta_expr[, names(tca_eigen)], tca_eigen)
c_res$gene <- rownames(c_res)
c_res      <- c_res[order(c_res$fdr), ]

# ---- Analysis D: Intersection – genes correlated with both delta insulin AND TCA ----
combined <- merge(
  b_res[, c("gene", "r", "fdr")],
  c_res[, c("gene", "r", "fdr")],
  by = "gene", suffixes = c("_insulin", "_tca")
)
combined$same_direction <- sign(combined$r_insulin) == sign(combined$r_tca)
combined$combined_score <- abs(combined$r_insulin) * abs(combined$r_tca)
combined$is_tca         <- combined$gene %in% tca_genes
combined                <- combined[order(-combined$combined_score), ]

# ---- Save ----
out_dir <- "/Users/willtrim/Documents/projs/bedrest/outputs/2_delta_space"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(as.data.table(tca_pheno_res), file.path(out_dir, "A_tca_eigen_vs_phenotypes.csv"))
fwrite(as.data.table(b_res),         file.path(out_dir, "B_genes_vs_delta_insulin.csv"))
fwrite(as.data.table(c_res),         file.path(out_dir, "C_genes_vs_tca_eigen.csv"))
fwrite(as.data.table(combined),      file.path(out_dir, "D_genes_insulin_and_tca.csv"))

# ---- Plot: TCA eigengene vs delta insulin ----
plot_df <- data.frame(
  participant   = names(tca_eigen),
  tca_eigen     = tca_eigen,
  delta_insulin = delta_insulin[names(tca_eigen)]
)
r_val <- cor(plot_df$tca_eigen, plot_df$delta_insulin, use = "complete.obs",
             method = "spearman")
p_val <- cor.test(plot_df$tca_eigen, plot_df$delta_insulin, method = "spearman",
                  exact = FALSE)$p.value

p <- ggplot(plot_df, aes(tca_eigen, delta_insulin, label = participant)) +
  geom_point(size = 3, color = "#3b6ea8") +
  geom_smooth(method = "lm", se = TRUE, color = "#a73bf5", alpha = 0.15) +
  labs(
    x     = sprintf("TCA module eigengene Δ (PC1, %.1f%% var)", pct_var),
    y     = "Δ Arterial insulin (mIU/L)",
    title = sprintf("rho = %.2f, p = %.3f", r_val, p_val)
  ) +
  theme_bw(base_size = 14)

# ggrepel for labels if available
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(size = 3)
} else {
  p <- p + geom_text(hjust = -0.2, size = 3)
}
ggsave(file.path(out_dir, "A_tca_eigen_vs_delta_insulin.pdf"), p, width = 6, height = 5)

message(sprintf(
  "Done.\n  B: %d genes with FDR<0.1 for delta insulin\n  C: %d genes with FDR<0.1 for TCA eigengene\n  Output: %s",
  sum(b_res$fdr < 0.1, na.rm = TRUE),
  sum(c_res$fdr < 0.1, na.rm = TRUE),
  out_dir
))
