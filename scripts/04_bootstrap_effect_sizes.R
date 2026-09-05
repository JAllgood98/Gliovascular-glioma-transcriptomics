# ============================================================
# Pathway-level grade comparison bootstrap analysis
# Manuscript:
# Integrated Transcriptomic Analysis Reveals Coordinated
# Gliovascular Remodeling Across Diffuse Glioma Progression
# ============================================================

set.seed(20260904)

B <- 2000

pathways <- c(
  "Endfoot_Organization",
  "Homeostatic_Astrocytes",
  "Reactive_Astrocytes",
  "BBB_Specialization",
  "Endothelial_Angiogenesis",
  "Basement_Membrane",
  "ECM_Remodeling",
  "Neuroinflammation"
)


# ------------------------------------------------------------
# Hedges' g
# Positive values indicate higher pathway activity in the
# later/higher-grade group.
# ------------------------------------------------------------

hedges_g <- function(x_early, x_late) {

  x_early <- x_early[is.finite(x_early)]
  x_late  <- x_late[is.finite(x_late)]

  n1 <- length(x_early)
  n2 <- length(x_late)

  pooled_sd <- sqrt(
    ((n1 - 1) * var(x_early) +
       (n2 - 1) * var(x_late)) /
      (n1 + n2 - 2)
  )

  d <- (mean(x_late) - mean(x_early)) / pooled_sd

  J <- 1 - (3 / (4 * (n1 + n2) - 9))

  g <- J * d

  return(g)
}


# ------------------------------------------------------------
# Bootstrap comparison
#
# Specimens are resampled with replacement separately within
# grades II, III, and IV.
#
# The same resampled grade III specimens are used for both
# II-vs-III and III-vs-IV comparisons within each bootstrap
# iteration.
#
# Delta_g = g(III-IV) - g(II-III)
# ------------------------------------------------------------

bootstrap_transition <- function(data, variable, B = 2000) {

  grade_II  <- data[data$Grade == "II", variable]
  grade_III <- data[data$Grade == "III", variable]
  grade_IV  <- data[data$Grade == "IV", variable]

  grade_II  <- grade_II[is.finite(grade_II)]
  grade_III <- grade_III[is.finite(grade_III)]
  grade_IV  <- grade_IV[is.finite(grade_IV)]

  n_II  <- length(grade_II)
  n_III <- length(grade_III)
  n_IV  <- length(grade_IV)

  g_II_III <- hedges_g(grade_II, grade_III)
  g_III_IV <- hedges_g(grade_III, grade_IV)

  delta_observed <- g_III_IV - g_II_III

  delta_boot <- numeric(B)

  for (b in seq_len(B)) {

    boot_II <- sample(
      grade_II,
      size = n_II,
      replace = TRUE
    )

    boot_III <- sample(
      grade_III,
      size = n_III,
      replace = TRUE
    )

    boot_IV <- sample(
      grade_IV,
      size = n_IV,
      replace = TRUE
    )

    g23 <- hedges_g(
      boot_II,
      boot_III
    )

    g34 <- hedges_g(
      boot_III,
      boot_IV
    )

    delta_boot[b] <- g34 - g23
  }

  ci <- quantile(
    delta_boot,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )

  # Two-sided sign-based bootstrap P value relative to Delta_g = 0
  p_boot <- 2 * min(
    mean(delta_boot <= 0, na.rm = TRUE),
    mean(delta_boot >= 0, na.rm = TRUE)
  )

  # Minimum reportable P value for B = 2000
  if (p_boot == 0) {
    p_boot <- 1 / B
  }

  data.frame(
    Variable = variable,
    N_II = n_II,
    N_III = n_III,
    N_IV = n_IV,
    g_II_III = g_II_III,
    g_III_IV = g_III_IV,
    Delta_g = delta_observed,
    CI_low = ci[1],
    CI_high = ci[2],
    P_boot = p_boot
  )
}


# ============================================================
# TCGA
# ============================================================

tcga <- read.csv(
  "processed_data/TCGA_analysis_cohort.csv",
  stringsAsFactors = FALSE
)

tcga_transition <- do.call(
  rbind,
  lapply(
    pathways,
    function(x) bootstrap_transition(
      data = tcga,
      variable = x,
      B = B
    )
  )
)

tcga_transition$Cohort <- "TCGA"

tcga_transition$P_BH <- p.adjust(
  tcga_transition$P_boot,
  method = "BH"
)

write.csv(
  tcga_transition,
  "results/grade_comparison_bootstrap_TCGA_reproduced.csv",
  row.names = FALSE
)


# ============================================================
# CGGA
# ============================================================

cgga <- read.csv(
  "processed_data/CGGA_analysis_cohort.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cgga_transition <- do.call(
  rbind,
  lapply(
    pathways,
    function(x) bootstrap_transition(
      data = cgga,
      variable = x,
      B = B
    )
  )
)

cgga_transition$Cohort <- "CGGA"

cgga_transition$P_BH <- p.adjust(
  cgga_transition$P_boot,
  method = "BH"
)

write.csv(
  cgga_transition,
  "results/grade_comparison_bootstrap_CGGA_reproduced.csv",
  row.names = FALSE
)


# ============================================================
# Print results
# ============================================================

cat("\n=== TCGA ===\n")
print(tcga_transition)

cat("\n=== CGGA ===\n")
print(cgga_transition)
