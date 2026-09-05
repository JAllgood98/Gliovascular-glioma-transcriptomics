# ============================================================
# Multivariable survival analysis and proportional-hazards
# diagnostics
#
# Manuscript:
# Integrated Transcriptomic Analysis Reveals Coordinated
# Gliovascular Remodeling Across Diffuse Glioma Progression
# ============================================================

library(survival)

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


# ============================================================
# Helper functions
# ============================================================

extract_pathway_cox <- function(fit, pathway, cohort, n_model) {

  s <- summary(fit)

  coef_table <- s$coefficients
  ci_table   <- s$conf.int

  data.frame(
    Cohort = cohort,
    Pathway = pathway,
    N = n_model,
    HR = ci_table[pathway, "exp(coef)"],
    CI_low = ci_table[pathway, "lower .95"],
    CI_high = ci_table[pathway, "upper .95"],
    P = coef_table[pathway, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}


extract_ph <- function(ph_list, cohort) {

  out <- do.call(
    rbind,
    lapply(names(ph_list), function(pathway) {

      tab <- as.data.frame(ph_list[[pathway]]$table)

      tab$Term <- rownames(tab)
      rownames(tab) <- NULL

      names(tab)[1:3] <- c(
        "chisq",
        "df",
        "P"
      )

      tab$Pathway_Model <- pathway
      tab$Cohort <- cohort

      tab[, c(
        "Cohort",
        "Pathway_Model",
        "Term",
        "chisq",
        "df",
        "P"
      )]
    })
  )

  rownames(out) <- NULL
  out
}


# ============================================================
# TCGA
#
# Model:
# Surv(OS_days, OS_status) ~ pathway + Grade + IDH_Group
# ============================================================

tcga <- read.csv(
  "processed_data/TCGA_survival_cohort.csv",
  stringsAsFactors = FALSE
)

tcga$Grade <- factor(
  tcga$Grade,
  levels = c("II", "III", "IV")
)

tcga$IDH_Group <- factor(tcga$IDH_Group)

tcga_cox <- list()
tcga_ph  <- list()
tcga_results <- list()


for (pathway in pathways) {

  vars <- c(
    "OS_days",
    "OS_status",
    "Grade",
    "IDH_Group",
    pathway
  )

  dat <- tcga[
    complete.cases(tcga[, vars]),
    vars
  ]

  form <- as.formula(
    paste0(
      "Surv(OS_days, OS_status) ~ ",
      pathway,
      " + Grade + IDH_Group"
    )
  )

  fit <- coxph(
    form,
    data = dat
  )

  ph <- cox.zph(fit)

  tcga_cox[[pathway]] <- fit
  tcga_ph[[pathway]] <- ph

  tcga_results[[pathway]] <- extract_pathway_cox(
    fit = fit,
    pathway = pathway,
    cohort = "TCGA",
    n_model = nrow(dat)
  )
}


tcga_survival_results <- do.call(
  rbind,
  tcga_results
)

rownames(tcga_survival_results) <- NULL

tcga_survival_results$P_BH <- p.adjust(
  tcga_survival_results$P,
  method = "BH"
)


tcga_ph_complete <- extract_ph(
  tcga_ph,
  "TCGA"
)


# ============================================================
# CGGA
#
# Model:
# Surv(OS, Censor) ~ pathway + Grade + IDH_Group + Age
# ============================================================

cgga <- read.csv(
  "processed_data/CGGA_analysis_cohort.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cgga$Grade <- factor(
  cgga$Grade,
  levels = c("II", "III", "IV")
)

cgga$IDH_Group <- factor(cgga$IDH_Group)

cgga_cox <- list()
cgga_ph  <- list()
cgga_results <- list()


for (pathway in pathways) {

  vars <- c(
    "OS",
    "Censor (alive=0; dead=1)",
    "Grade",
    "IDH_Group",
    "Age",
    pathway
  )

  dat <- cgga[
    complete.cases(cgga[, vars]),
    vars
  ]

  # Rename censor variable locally for formula simplicity
  names(dat)[
    names(dat) == "Censor (alive=0; dead=1)"
  ] <- "Censor"

  form <- as.formula(
    paste0(
      "Surv(OS, Censor) ~ ",
      pathway,
      " + Grade + IDH_Group + Age"
    )
  )

  fit <- coxph(
    form,
    data = dat
  )

  ph <- cox.zph(fit)

  cgga_cox[[pathway]] <- fit
  cgga_ph[[pathway]] <- ph

  cgga_results[[pathway]] <- extract_pathway_cox(
    fit = fit,
    pathway = pathway,
    cohort = "CGGA",
    n_model = nrow(dat)
  )
}


cgga_survival_results <- do.call(
  rbind,
  cgga_results
)

rownames(cgga_survival_results) <- NULL

cgga_survival_results$P_BH <- p.adjust(
  cgga_survival_results$P,
  method = "BH"
)


cgga_ph_complete <- extract_ph(
  cgga_ph,
  "CGGA"
)


# ============================================================
# Save results
# ============================================================

write.csv(
  tcga_survival_results,
  "results/survival_results_TCGA_reproduced.csv",
  row.names = FALSE
)

write.csv(
  cgga_survival_results,
  "results/survival_results_CGGA_reproduced.csv",
  row.names = FALSE
)

write.csv(
  tcga_ph_complete,
  "results/proportional_hazards_complete_TCGA_reproduced.csv",
  row.names = FALSE
)

write.csv(
  cgga_ph_complete,
  "results/proportional_hazards_complete_CGGA_reproduced.csv",
  row.names = FALSE
)


# ============================================================
# Print results
# ============================================================

cat("\n=== TCGA SURVIVAL RESULTS ===\n")
print(tcga_survival_results)

cat("\n=== CGGA SURVIVAL RESULTS ===\n")
print(cgga_survival_results)

cat("\n=== TCGA PH DIAGNOSTICS ===\n")
print(tcga_ph_complete)

cat("\n=== CGGA PH DIAGNOSTICS ===\n")
print(cgga_ph_complete)
