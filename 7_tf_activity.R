
# Transcription factor activity inference
#
# Uses decoupleR with DoRothEA regulons (confidence levels A, B, C).
# Method: ULM (univariate linear model), equivalent to the VIPER algorithm.
#
# Analyses:
#   1. TF activity matrix per sample (pre + post)
#   2. Differential TF activity: pre vs post (paired Wilcoxon, FDR corrected)
#   3. Correlation of TF activities with TCA module eigengene (post samples)
#   4. Correlation of TF activities with delta arterial insulin
#
# The TCA module eigengene is recomputed here from post-sample expression;
# this is the "absolute" version (vs delta in script 2), revealing which TFs
# track the TCA–insulin anticorrelation in the post state.
#
# Installation:
#   BiocManager::install("decoupleR")
#   BiocManager::install("OmnipathR")   # only needed if get_dorothea fails

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
if (!exists("counts")) source("analysis_preamble.R")

if (!requireNamespace("decoupleR", quietly = TRUE))
  stop("Install decoupleR: BiocManager::install('decoupleR')")
library(decoupleR)

OUT_DIR <- "/Users/willtrim/Documents/projs/bedrest/outputs/7_tf_activity"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =====================================================================
# 1. Download DoRothEA regulons
# =====================================================================
message("Fetching DoRothEA regulons (A/B/C confidence)...")
net <- tryCatch(
  get_dorothea(organism = "human", levels = c("A", "B", "C")),
  error = function(e) {
    message("get_dorothea failed; trying OmnipathR fallback...")
    if (!requireNamespace("OmnipathR", quietly = TRUE))
      stop("Install OmnipathR: BiocManager::install('OmnipathR')")
    OmnipathR::import_tf_target_interactions(
      resources         = "DoRothEA",
      organism          = 9606,
      dorothea_levels   = c("A", "B", "C")
    ) %>% rename(source = source_genesymbol, target = target_genesymbol,
                 mor = is_stimulation) %>%
      mutate(mor = ifelse(mor == 1, 1, -1))
  }
)
message(sprintf("Regulons loaded: %d TFs, %d TF-target pairs",
                length(unique(net$source)), nrow(net)))

# =====================================================================
# 2. Run ULM per sample (counts_sym has gene symbol rownames)
# =====================================================================
message("Running ULM TF activity inference...")
acts_long <- run_ulm(
  mat     = as.matrix(counts_sym),
  net     = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

# Pivot to TFs × samples matrix
tf_mat <- acts_long %>%
  filter(statistic == "ulm") %>%
  pivot_wider(id_cols = "source", names_from = "condition", values_from = "score") %>%
  column_to_rownames("source") %>%
  as.matrix()

fwrite(as.data.table(tf_mat, keep.rownames = "TF"),
       file.path(OUT_DIR, "tf_activity_all_samples.csv"))

# =====================================================================
# 3. Differential TF activity: Pre vs Post (paired Wilcoxon)
# =====================================================================
pre_ids  <- sub("_Pre$",  "", grep("_Pre$",  colnames(tf_mat), value = TRUE))
post_ids <- sub("_Post$", "", grep("_Post$", colnames(tf_mat), value = TRUE))
paired   <- intersect(pre_ids, post_ids)

tf_pre  <- tf_mat[, paste0(paired, "_Pre")]
tf_post <- tf_mat[, paste0(paired, "_Post")]
colnames(tf_pre) <- colnames(tf_post) <- paired

diff_res <- do.call(rbind, lapply(rownames(tf_mat), function(tf) {
  pre_v  <- tf_pre[tf,  ]
  post_v <- tf_post[tf, ]
  tt     <- wilcox.test(post_v, pre_v, paired = TRUE, exact = FALSE)
  data.frame(
    TF            = tf,
    median_pre    = median(pre_v),
    median_post   = median(post_v),
    delta_activity = median(post_v - pre_v),
    p_value       = tt$p.value,
    stringsAsFactors = FALSE
  )
}))
diff_res$fdr <- p.adjust(diff_res$p_value, "BH")
diff_res      <- diff_res[order(diff_res$p_value), ]
fwrite(as.data.table(diff_res), file.path(OUT_DIR, "differential_tf_activity.csv"))

# =====================================================================
# 4. Correlate TF activity (post) with TCA module eigengene (post)
# =====================================================================
post_tca_mat  <- t(counts[tca_genes, grep("_Post$", colnames(counts), value = TRUE)])
rownames(post_tca_mat) <- sub("_Post$", "", rownames(post_tca_mat))
pca_tca  <- prcomp(post_tca_mat, scale. = TRUE)
tca_eigen_post <- pca_tca$x[, 1]
if (cor(tca_eigen_post, rowMeans(post_tca_mat)) < 0) tca_eigen_post <- -tca_eigen_post
pct_var  <- round(100 * pca_tca$sdev[1]^2 / sum(pca_tca$sdev^2), 1)

post_tf  <- tf_mat[, paste0(rownames(post_tca_mat), "_Post"),  drop = FALSE]
colnames(post_tf) <- rownames(post_tca_mat)

tf_tca_cors <- apply(post_tf, 1, function(x)
  cor(x, tca_eigen_post[names(x)], method = "spearman", use = "complete.obs"))
tf_tca_pval <- apply(post_tf, 1, function(x) {
  y <- tca_eigen_post[names(x)]
  ok <- !is.na(x) & !is.na(y)
  cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)$p.value
})
tca_tf_res <- data.frame(
  TF          = names(tf_tca_cors),
  r_tca_eigen = tf_tca_cors,
  p_value     = tf_tca_pval,
  fdr         = p.adjust(tf_tca_pval, "BH")
)
tca_tf_res <- tca_tf_res[order(tca_tf_res$p_value), ]
fwrite(as.data.table(tca_tf_res), file.path(OUT_DIR, "tf_vs_tca_eigen_post.csv"))

# =====================================================================
# 5. Correlate TF activity (delta) with delta arterial insulin
# =====================================================================
delta_insulin <- setNames(delta_pheno[, "Art..Insulin..mIU.L."], rownames(delta_pheno))
tf_delta      <- tf_post[, pheno_paired] - tf_pre[, pheno_paired]
ok_ins        <- names(delta_insulin)[!is.na(delta_insulin)]
shared_ins    <- intersect(colnames(tf_delta), ok_ins)

tf_ins_cors <- apply(tf_delta[, shared_ins], 1, function(x)
  cor(x, delta_insulin[shared_ins], method = "spearman", use = "complete.obs"))
tf_ins_pval <- apply(tf_delta[, shared_ins], 1, function(x) {
  ok <- !is.na(x)
  cor.test(x[ok], delta_insulin[shared_ins][ok], method = "spearman", exact = FALSE)$p.value
})
ins_tf_res <- data.frame(
  TF              = names(tf_ins_cors),
  r_delta_insulin = tf_ins_cors,
  p_value         = tf_ins_pval,
  fdr             = p.adjust(tf_ins_pval, "BH")
)
ins_tf_res <- ins_tf_res[order(ins_tf_res$p_value), ]
fwrite(as.data.table(ins_tf_res), file.path(OUT_DIR, "tf_vs_delta_insulin.csv"))

# =====================================================================
# 6. Summary plot: TFs linking TCA eigengene and delta insulin
# =====================================================================
combined_tf <- merge(
  tca_tf_res[, c("TF", "r_tca_eigen", "fdr")],
  ins_tf_res[, c("TF", "r_delta_insulin", "fdr")],
  by = "TF", suffixes = c("_tca", "_ins")
)
combined_tf$combined_score <- abs(combined_tf$r_tca_eigen) * abs(combined_tf$r_delta_insulin)
combined_tf <- combined_tf[order(-combined_tf$combined_score), ]
fwrite(as.data.table(combined_tf), file.path(OUT_DIR, "tf_combined_tca_insulin.csv"))

# Volcano-style: differential activity vs TCA correlation
plot_df <- merge(diff_res[, c("TF", "delta_activity", "fdr")],
                 tca_tf_res[, c("TF", "r_tca_eigen")], by = "TF")
plot_df$label <- ifelse(plot_df$fdr < 0.1 | abs(plot_df$r_tca_eigen) > 0.5,
                        plot_df$TF, NA)

p <- ggplot(plot_df, aes(delta_activity, r_tca_eigen, color = -log10(fdr + 1e-6))) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "grey70", high = "#a73bf5", name = "-log10(FDR)") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  labs(x = "Delta TF activity (Post - Pre)",
       y = sprintf("r with TCA eigengene (post, PC1 %.1f%%)", pct_var)) +
  theme_bw(base_size = 13)

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(aes(label = label), size = 2.5, na.rm = TRUE,
                                     max.overlaps = 20)
} else {
  p <- p + geom_text(aes(label = label), size = 2.5, na.rm = TRUE, hjust = -0.1)
}
ggsave(file.path(OUT_DIR, "tf_activity_vs_tca_eigen.pdf"), p, width = 8, height = 6)

message(sprintf(
  "Done.\n  Differential TFs (FDR<0.1): %d\n  TFs correlated with TCA eigengene (FDR<0.1): %d\n  Output: %s",
  sum(diff_res$fdr < 0.1, na.rm = TRUE),
  sum(tca_tf_res$fdr < 0.1, na.rm = TRUE),
  OUT_DIR
))
