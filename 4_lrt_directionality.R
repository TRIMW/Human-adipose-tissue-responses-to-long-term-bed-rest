
# LRT with directionality and FDR correction
#
# Improved version of the original nb_glm_tests.R. For each TCA gene as anchor:
#   base model: insulin ~ TCA_gene
#   test model: insulin ~ TCA_gene + other_gene
#
# Added vs original:
#   - BH FDR correction
#   - Coefficient value and direction for other_gene
#   - McFadden's delta pseudo-R² (effect size beyond the TCA gene alone)
#   - Coefficient of TCA gene before and after (does adding other_gene change TCA effect?)
#   - Parallelization via parallel::mclapply (set N_CORES below)
#   - Convergence filtering (failed/non-converged models are excluded)

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})
if (!exists("counts")) source("analysis_preamble.R")

N_CORES    <- max(1L, detectCores() - 1L)
PHENO_NAME <- "Art..Insulin..mIU.L."
OUT_DIR    <- "/Users/willtrim/Documents/projs/bedrest/outputs/4_lrt_directionality"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Build merged counts + phenotype data frame (post only) ----
pheno_post <- pheno[pheno$Time.Point == "Post Bed Rest", ]
rownames(pheno_post) <- pheno_post$Participant.ID
post_counts_t <- as.data.frame(t(counts[, grep("_Post$", colnames(counts))]))
rownames(post_counts_t) <- sub("_Post$", "", rownames(post_counts_t))
counts_pheno  <- merge(post_counts_t,
                       pheno_post[, PHENO_NAME, drop = FALSE],
                       by = "row.names")
rownames(counts_pheno) <- counts_pheno$Row.names
counts_pheno$Row.names <- NULL

# Null model log-likelihood for pseudo-R²
null_model  <- glm(reformulate("1", response = PHENO_NAME),
                   family = Gamma(link = "log"), data = counts_pheno)
loglik_null <- as.numeric(logLik(null_model))

# ---- Worker function for one TCA gene ----
run_one_tca <- function(tca) {
  base_form  <- reformulate(tca, response = PHENO_NAME)
  base_model <- tryCatch(
    glm(base_form, family = Gamma(link = "log"), data = counts_pheno),
    error = function(e) NULL
  )
  if (is.null(base_model) || !base_model$converged) return(NULL)

  loglik_base <- as.numeric(logLik(base_model))
  coef_tca_base <- coef(base_model)[tca]

  results <- mclapply(other_genes, function(other) {
    full_form  <- reformulate(c(tca, other), response = PHENO_NAME)
    full_model <- tryCatch(
      glm(full_form, family = Gamma(link = "log"), data = counts_pheno),
      error = function(e) NULL
    )
    if (is.null(full_model) || !full_model$converged) return(NULL)

    lrt         <- anova(base_model, full_model, test = "Chisq")
    coef_other  <- coef(full_model)[other]
    coef_tca_full <- coef(full_model)[tca]
    loglik_full <- as.numeric(logLik(full_model))

    data.frame(
      gene            = other,
      lrt_pvalue      = lrt$`Pr(>Chi)`[2],
      delta_aic       = AIC(base_model) - AIC(full_model),
      coef_other      = coef_other,
      direction       = ifelse(coef_other > 0, "positive", "negative"),
      coef_tca_base   = coef_tca_base,
      coef_tca_full   = coef_tca_full,
      # McFadden pseudo-R² gain from adding other_gene
      delta_pseudo_r2 = (1 - loglik_full / loglik_null) -
                        (1 - loglik_base / loglik_null),
      stringsAsFactors = FALSE
    )
  }, mc.cores = N_CORES)

  results <- rbindlist(Filter(Negate(is.null), results))
  results[, fdr := p.adjust(lrt_pvalue, method = "BH")]
  results[, tca_gene := tca]
  setcolorder(results, c("tca_gene", "gene"))
  results[order(fdr)]
}

# ---- Run for all TCA genes ----
all_results <- list()
for (tca in tca_genes) {
  message(sprintf("  Processing TCA anchor: %s", tca))
  res <- run_one_tca(tca)
  if (!is.null(res)) {
    all_results[[tca]] <- res
    fwrite(res, file.path(OUT_DIR, paste0(tca, ".csv")))
  }
}

message(sprintf("Done. Output: %s", OUT_DIR))
