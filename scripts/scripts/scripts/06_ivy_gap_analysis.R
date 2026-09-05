# ============================================================
# Ivy Glioblastoma Atlas Project regional analysis
#
# Manuscript:
# Integrated Transcriptomic Analysis Reveals Coordinated
# Gliovascular Remodeling Across Diffuse Glioma Progression
#
# Nine sentinel gliovascular genes were evaluated across
# anatomically annotated GBM regions using linear mixed-effects
# models accounting for repeated regional sampling within tumors.
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

library(lme4)
library(lmerTest)
library(emmeans)


# ============================================================
# 2. Load processed Ivy GAP data
# ============================================================

ivy <- read.csv(
  "processed_data/IvyGAP_sentinel_gene_expression.csv",
  stringsAsFactors = FALSE
)

sentinel_genes <- c(
  "KCNJ10",
  "SNTA1",
  "AQP4",
  "SERPINA3",
  "CCL2",
  "FN1",
  "COL4A1",
  "ANGPT2",
  "TJP1"
)


# ============================================================
# 3. Data checks
# ============================================================

ivy$Region <- factor(
  ivy$Region,
  levels = c(
    "LE",
    "IT",
    "CT",
    "MVP",
    "PAN"
  )
)

ivy$tumor_id <- factor(ivy$tumor_id)

stopifnot(
  all(sentinel_genes %in% unique(ivy$gene_symbol))
)

# Verify stored expression transformation
transformation_difference <- max(
  abs(
    ivy$expression -
      log2(ivy$FPKM + 1)
  ),
  na.rm = TRUE
)

if (transformation_difference > 1e-10) {
  stop(
    "Stored expression values do not equal log2(FPKM + 1)."
  )
}

cat("Observations:", nrow(ivy), "\n")
cat(
  "Sentinel genes:",
  length(unique(ivy$gene_symbol)),
  "\n"
)
cat(
  "Tumors:",
  length(unique(ivy$tumor_id)),
  "\n\n"
)

cat("Regional samples per gene:\n")
print(
  table(ivy$Region) /
    length(unique(ivy$gene_symbol))
)


# ============================================================
# 4. Linear mixed-effects models
#
# For each sentinel gene:
#
# expression ~ Region + (1 | tumor_id)
#
# Region is modeled as a fixed effect and tumor donor as a
# random intercept to account for repeated anatomical sampling
# from the same tumor.
# ============================================================

model_results <- list()
pairwise_results <- list()

for (gene in sentinel_genes) {

  dat <- ivy[
    ivy$gene_symbol == gene &
      complete.cases(
        ivy[, c(
          "expression",
          "Region",
          "tumor_id"
        )]
      ),
  ]

  dat$Region <- droplevels(dat$Region)

  model <- lmer(
    expression ~ Region + (1 | tumor_id),
    data = dat
  )

  # Omnibus test for anatomical region
  aov_table <- anova(model)

  model_results[[gene]] <- data.frame(
    Gene = gene,
    N = nrow(dat),
    Tumors = length(unique(dat$tumor_id)),
    F_value = aov_table["Region", "F value"],
    df_num = aov_table["Region", "NumDF"],
    df_den = aov_table["Region", "DenDF"],
    P_value = aov_table["Region", "Pr(>F)"],
    stringsAsFactors = FALSE
  )

  # Estimated marginal means with Tukey-adjusted pairwise tests
  emm <- emmeans(
    model,
    ~ Region
  )

  pairs_gene <- as.data.frame(
    pairs(
      emm,
      adjust = "tukey"
    )
  )

  pairs_gene$Gene <- gene

  pairwise_results[[gene]] <- pairs_gene
}


# ============================================================
# 5. Combine results
# ============================================================

model_results <- do.call(
  rbind,
  model_results
)

pairwise_results <- do.call(
  rbind,
  pairwise_results
)

rownames(model_results) <- NULL
rownames(pairwise_results) <- NULL


# ============================================================
# 6. Save reproduced results
# ============================================================

write.csv(
  model_results,
  "results/IvyGAP_mixed_models_reproduced.csv",
  row.names = FALSE
)

write.csv(
  pairwise_results,
  "results/IvyGAP_pairwise_Tukey_reproduced.csv",
  row.names = FALSE
)


# ============================================================
# 7. Display results
# ============================================================

cat("\n=== IVY GAP OMNIBUS MIXED MODELS ===\n")
print(model_results)

cat("\n=== IVY GAP TUKEY PAIRWISE COMPARISONS ===\n")
print(pairwise_results)
