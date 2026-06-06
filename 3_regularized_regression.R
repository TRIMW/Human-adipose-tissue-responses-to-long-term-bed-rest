
# Regularized regression (elastic net) – post bed rest
#
# For each TCA gene as a fixed (unpenalized) anchor:
#   insulin ~ TCA_gene + all_other_genes
#   where TCA_gene has penalty.factor = 0 (always included)
#   and all other genes are penalized (L1+L2 via elastic net, alpha=0.5)
#
# Two outputs per TCA gene:
#   1. Single-run: genes selected at lambda.1se (cross-validated)
#   2. Stability: fraction of bootstrap subsamples in which each gene is selected
#      (subsampling 80% without replacement, 100 iterations)
#
# Pre-filtering to the 5000 most variable genes across post samples is applied
# to keep runtime feasible. Change N_HVG to include more genes at the cost of speed.

suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(data.table))
if (!exists("counts")) source("analysis_preamble.R")

N_HVG    <- 5000   # number of highly variable genes to include as candidates
N_BOOT   <- 100    # bootstrap iterations for stability selection
SUBSAMP  <- 0.80   # fraction of samples per bootstrap
ALPHA    <- 0.5    # elastic net mixing (0=ridge, 1=lasso)
LAMBDA_S <- "lambda.1se"  # use lambda.min for more variables

# ---- Prepare post-sample data ----
pheno_post <- pheno[pheno$Time.Point == "Post Bed Rest", ]
rownames(pheno_post) <- pheno_post$Participant.ID
shared_ids <- intersect(rownames(post_mat), rownames(pheno_post))
shared_ids <- shared_ids[!is.na(pheno_post[shared_ids, "Art..Insulin..mIU.L."])]
y          <- as.numeric(pheno_post[shared_ids, "Art..Insulin..mIU.L."])
names(y)   <- shared_ids
message(sprintf("Post samples with insulin data: %d", length(y)))

# ---- Select highly variable other genes ----
other_in_data <- other_genes[other_genes %in% colnames(post_mat)]
gene_vars     <- apply(post_mat[shared_ids, other_in_data], 2, var)
hvg           <- names(sort(gene_vars, decreasing = TRUE))[1:min(N_HVG, length(gene_vars))]
message(sprintf("Using top %d HVGs as candidate predictors", length(hvg)))

# ---- Run for each TCA gene ----
out_dir <- "/Users/willtrim/Documents/projs/bedrest/outputs/3_regularized_regression"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

all_stability <- list()

for (tca in tca_genes) {
  message(sprintf("  Elastic net anchor: %s", tca))

  # X matrix: TCA gene (col 1, unpenalized) + HVGs
  gene_cols <- c(tca, hvg)
  X_raw     <- post_mat[shared_ids, gene_cols]
  X_scaled  <- scale(X_raw)

  pf <- c(0, rep(1, length(hvg)))   # 0 = no penalty on TCA anchor

  # Single run with cross-validation
  set.seed(42)
  cv_fit <- cv.glmnet(X_scaled, y, alpha = ALPHA, penalty.factor = pf,
                      nfolds = min(5, length(y)), family = "gaussian")

  coefs     <- coef(cv_fit, s = LAMBDA_S)[-1, 1]   # drop intercept
  selected  <- coefs[coefs != 0]
  selected  <- selected[names(selected) != tca]     # exclude the anchor itself

  single_run <- data.frame(
    tca_gene    = tca,
    gene        = names(selected),
    coefficient = selected,
    row.names   = NULL
  )

  # Bootstrap stability selection
  n_sub     <- floor(length(y) * SUBSAMP)
  sel_count <- integer(length(hvg))
  names(sel_count) <- hvg

  for (b in seq_len(N_BOOT)) {
    idx    <- sample(seq_along(y), n_sub, replace = FALSE)
    X_b    <- X_scaled[idx, ]
    y_b    <- y[idx]
    cv_b   <- tryCatch(
      cv.glmnet(X_b, y_b, alpha = ALPHA, penalty.factor = pf,
                nfolds = min(5, length(y_b)), family = "gaussian"),
      error = function(e) NULL
    )
    if (is.null(cv_b)) next
    coefs_b <- coef(cv_b, s = LAMBDA_S)[-1, 1]
    coefs_b <- coefs_b[names(coefs_b) != tca]
    sel_count <- sel_count + (coefs_b != 0)[names(sel_count)]
  }

  stability <- data.frame(
    tca_gene         = tca,
    gene             = names(sel_count),
    stability_freq   = sel_count / N_BOOT,
    single_run_coef  = coefs[names(sel_count)],
    row.names        = NULL
  )
  stability <- stability[order(-stability$stability_freq), ]
  all_stability[[tca]] <- stability

  fwrite(single_run, file.path(out_dir, paste0(tca, "_singlerun.csv")))
  fwrite(stability,  file.path(out_dir, paste0(tca, "_stability.csv")))
}

# Combined stability table (genes selected for ≥1 TCA gene)
combined_stab <- rbindlist(all_stability)
combined_stab <- combined_stab[stability_freq > 0][order(-stability_freq)]
fwrite(combined_stab, file.path(out_dir, "all_tca_stability.csv"))

message(sprintf("Done. Stability results (freq>0.5): %d gene-TCA pairs  |  Output: %s",
                sum(combined_stab$stability_freq > 0.5), out_dir))
