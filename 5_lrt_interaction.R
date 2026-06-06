
# LRT with interaction term
#
# For each TCA gene as anchor, fits three nested Gamma GLMs:
#   M0: insulin ~ TCA_gene                           (base)
#   M1: insulin ~ TCA_gene + other_gene              (additive)
#   M2: insulin ~ TCA_gene * other_gene              (full interaction)
#
# Comparisons:
#   M1 vs M0 -> additive effect of other_gene
#   M2 vs M1 -> interaction beyond additive
#   M2 vs M0 -> combined test (use for ranking; most powerful overall)
#
# The interaction coefficient direction tells you whether other_gene
# amplifies (+) or dampens (-) the TCA gene effect on insulin.
#
# Interaction hits are more biologically specific than purely additive hits:
# they mean the TCA–insulin relationship is *modulated* by the other gene.

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})
if (!exists("counts")) source("analysis_preamble.R")

N_CORES    <- max(1L, detectCores() - 1L)
PHENO_NAME <- "Art..Insulin..mIU.L."
OUT_DIR    <- "/Users/willtrim/Documents/projs/bedrest/outputs/5_lrt_interaction"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Merged data frame (post samples) ----
pheno_post <- pheno[pheno$Time.Point == "Post Bed Rest", ]
rownames(pheno_post) <- pheno_post$Participant.ID
post_counts_t <- as.data.frame(t(counts[, grep("_Post$", colnames(counts))]))
rownames(post_counts_t) <- sub("_Post$", "", rownames(post_counts_t))
counts_pheno  <- merge(post_counts_t,
                       pheno_post[, PHENO_NAME, drop = FALSE],
                       by = "row.names")
rownames(counts_pheno) <- counts_pheno$Row.names
counts_pheno$Row.names <- NULL

run_one_tca <- function(tca) {
  m0 <- tryCatch(
    glm(reformulate(tca, response = PHENO_NAME),
        family = Gamma(link = "log"), data = counts_pheno),
    error = function(e) NULL
  )
  if (is.null(m0) || !m0$converged) return(NULL)

  results <- mclapply(other_genes, function(other) {
    m1 <- tryCatch(
      glm(reformulate(c(tca, other), response = PHENO_NAME),
          family = Gamma(link = "log"), data = counts_pheno),
      error = function(e) NULL
    )
    m2 <- tryCatch(
      glm(as.formula(paste(PHENO_NAME, "~", tca, "*", other)),
          family = Gamma(link = "log"), data = counts_pheno),
      error = function(e) NULL
    )
    if (is.null(m1) || !m1$converged) return(NULL)
    if (is.null(m2) || !m2$converged) return(NULL)

    lrt_additive    <- anova(m0, m1, test = "Chisq")
    lrt_interaction <- anova(m1, m2, test = "Chisq")
    lrt_combined    <- anova(m0, m2, test = "Chisq")

    # Interaction term name in glm output: "tca:other" or "other:tca"
    int_term <- grep(":", names(coef(m2)), value = TRUE)[1]
    int_coef <- if (!is.na(int_term)) coef(m2)[int_term] else NA_real_

    data.frame(
      gene                = other,
      p_additive          = lrt_additive$`Pr(>Chi)`[2],
      p_interaction_only  = lrt_interaction$`Pr(>Chi)`[2],
      p_combined          = lrt_combined$`Pr(>Chi)`[2],
      delta_aic_m0_m1     = AIC(m0) - AIC(m1),
      delta_aic_m1_m2     = AIC(m1) - AIC(m2),
      delta_aic_m0_m2     = AIC(m0) - AIC(m2),
      interaction_coef    = int_coef,
      interaction_dir     = ifelse(is.na(int_coef), NA_character_,
                                   ifelse(int_coef > 0, "amplifying", "dampening")),
      coef_other_m1       = coef(m1)[other],
      stringsAsFactors    = FALSE
    )
  }, mc.cores = N_CORES)

  results <- rbindlist(Filter(Negate(is.null), results))
  results[, fdr_additive   := p.adjust(p_additive,   method = "BH")]
  results[, fdr_interaction := p.adjust(p_interaction_only, method = "BH")]
  results[, fdr_combined   := p.adjust(p_combined,   method = "BH")]
  results[, tca_gene := tca]
  setcolorder(results, c("tca_gene", "gene"))
  results[order(fdr_combined)]
}

all_results <- list()
for (tca in tca_genes) {
  message(sprintf("  Interaction test anchor: %s", tca))
  res <- run_one_tca(tca)
  if (!is.null(res)) {
    all_results[[tca]] <- res
    fwrite(res, file.path(OUT_DIR, paste0(tca, ".csv")))
  }
}

message(sprintf("Done. Output: %s", OUT_DIR))
