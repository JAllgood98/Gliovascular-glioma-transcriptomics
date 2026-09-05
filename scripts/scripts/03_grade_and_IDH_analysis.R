# ============================================================
# Grade- and IDH-associated pathway comparisons
#
# Manuscript:
# Integrated Transcriptomic Analysis Reveals Coordinated
# Gliovascular Remodeling Across Diffuse Glioma Progression
#
# This script reproduces the descriptive pathway comparisons
# across WHO grade and IDH groups in TCGA and CGGA.
# ============================================================


# ============================================================
# 1. Load processed analytical cohorts
# ============================================================

tcga <- read.csv(
  "processed_data/TCGA_analysis_cohort.csv",
  stringsAsFactors = FALSE
)

cgga <- read.csv(
  "processed_data/CGGA_analysis_cohort.csv",
  stringsAsFactors = FALSE
)


# ============================================================
# 2. Pathway names
# ============================================================

pathways <- c(
  "Homeostatic_Astrocyte",
  "Reactive_Astrocyte",
  "Neuroinflammation",
  "Endothelial_Angiogenesis",
  "Basement_Membrane",
  "Extracellular_Matrix",
  "Astrocyte_Endfoot",
  "BBB_Specialization"
)


# ============================================================
# 3. Helper: grade comparison
# ============================================================

run_grade_analysis <- function(dat, cohort_name) {

  dat$Grade <- factor(
    dat$Grade,
    levels = c("II", "III", "IV")
  )

  results <- list()

  for (pathway in pathways) {

    model <- aov(
      reformulate(
        "Grade",
        response = pathway
      ),
      data = dat
    )

    aov_tab <- summary(model)[[1]]

    tuk <- TukeyHSD(model, "Grade")$Grade

    results[[pathway]] <- list(
      omnibus = data.frame(
        Cohort = cohort_name,
        Pathway = pathway,
        F_value = aov_tab["Grade", "F value"],
        df_between = aov_tab["Grade", "Df"],
        df_within = aov_tab["Residuals", "Df"],
        P_value = aov_tab["Grade", "Pr(>F)"]
      ),
      tukey = data.frame(
        Cohort = cohort_name,
        Pathway = pathway,
        Comparison = rownames(tuk),
        Difference = tuk[, "diff"],
        CI_low = tuk[, "lwr"],
        CI_high = tuk[, "upr"],
        P_adjusted = tuk[, "p adj"],
        row.names = NULL
      )
    )
  }

  omnibus <- do.call(
    rbind,
    lapply(results, `[[`, "omnibus")
  )

  tukey <- do.call(
    rbind,
    lapply(results, `[[`, "tukey")
  )

  rownames(omnibus) <- NULL
  rownames(tukey) <- NULL

  list(
    omnibus = omnibus,
    tukey = tukey
  )
}


# ============================================================
# 4. Helper: IDH comparison
# ============================================================

run_idh_analysis <- function(dat, cohort_name) {

  dat <- dat[
    !is.na(dat$IDH_Group) &
      dat$IDH_Group != "",
  ]

  results <- lapply(
    pathways,
    function(pathway) {

      form <- reformulate(
        "IDH_Group",
        response = pathway
      )

      fit <- t.test(
        form,
        data = dat
      )

      data.frame(
        Cohort = cohort_name,
        Pathway = pathway,
        Group1 = names(fit$estimate)[1],
        Mean1 = unname(fit$estimate[1]),
        Group2 = names(fit$estimate)[2],
        Mean2 = unname(fit$estimate[2]),
        Difference = unname(
          fit$estimate[1] -
            fit$estimate[2]
        ),
        CI_low = fit$conf.int[1],
        CI_high = fit$conf.int[2],
        T_value = unname(fit$statistic),
        df = unname(fit$parameter),
        P_value = fit$p.value
      )
    }
  )

  out <- do.call(
    rbind,
    results
  )

  out$P_BH <- p.adjust(
    out$P_value,
    method = "BH"
  )

  rownames(out) <- NULL

  out
}


# ============================================================
# 5. Run TCGA analyses
# ============================================================

tcga_grade <- run_grade_analysis(
  tcga,
  "TCGA"
)

tcga_idh <- run_idh_analysis(
  tcga,
  "TCGA"
)


# ============================================================
# 6. Run CGGA analyses
# ============================================================

cgga_grade <- run_grade_analysis(
  cgga,
  "CGGA"
)

cgga_idh <- run_idh_analysis(
  cgga,
  "CGGA"
)


# ============================================================
# 7. Save reproduced results
# ============================================================

write.csv(
  tcga_grade$omnibus,
  "results/grade_omnibus_TCGA_reproduced.csv",
  row.names = FALSE
)

write.csv(
  tcga_grade$tukey,
  "results/grade_Tukey_TCGA_reproduced.csv",
  row.names = FALSE
)

write.csv(
  tcga_idh,
  "results/IDH_comparison_TCGA_reproduced.csv",
  row.names = FALSE
)

write.csv(
  cgga_grade$omnibus,
  "results/grade_omnibus_CGGA_reproduced.csv",
  row.names = FALSE
)

write.csv(
  cgga_grade$tukey,
  "results/grade_Tukey_CGGA_reproduced.csv",
  row.names = FALSE
)

write.csv(
  cgga_idh,
  "results/IDH_comparison_CGGA_reproduced.csv",
  row.names = FALSE
)


# ============================================================
# 8. Print summary
# ============================================================

cat("\n=== TCGA GRADE OMNIBUS ===\n")
print(tcga_grade$omnibus)

cat("\n=== TCGA IDH COMPARISONS ===\n")
print(tcga_idh)

cat("\n=== CGGA GRADE OMNIBUS ===\n")
print(cgga_grade$omnibus)

cat("\n=== CGGA IDH COMPARISONS ===\n")
print(cgga_idh)
