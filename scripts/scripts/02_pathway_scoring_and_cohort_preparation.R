# ============================================================
# Gliovascular pathway scoring
#
# Manuscript:
# Integrated Transcriptomic Analysis Reveals Coordinated
# Gliovascular Remodeling Across Diffuse Glioma Progression
#
# This script defines the pathway-scoring procedure applied
# independently within each cohort.
#
# Input expression matrices are expected to contain:
#   rows    = genes
#   columns = specimens
#
# Raw source data are not redistributed in this repository.
# ============================================================


# ============================================================
# 1. Curated pathway definitions
# ============================================================

gene_sets <- list(

  Homeostatic_Astrocytes = c(
    "ALDH1L1",
    "SOX9",
    "SLC1A2",
    "SLC1A3",
    "GLUL",
    "KCNJ10"
  ),

  Reactive_Astrocytes = c(
    "GFAP",
    "VIM",
    "GJA1",
    "S100B",
    "SERPINA3",
    "LCN2"
  ),

  Neuroinflammation = c(
    "IL1B",
    "TNF",
    "CCL2",
    "CXCL10",
    "TREM2",
    "TYROBP",
    "CSF1R",
    "AIF1",
    "CD68",
    "C1QA",
    "C1QB",
    "C1QC"
  ),

  Endothelial_Angiogenesis = c(
    "PECAM1",
    "VWF",
    "EMCN",
    "KDR",
    "ANGPT2",
    "PLVAP"
  ),

  Basement_Membrane = c(
    "COL4A1",
    "COL4A2",
    "LAMA2",
    "LAMB2",
    "LAMC1",
    "HSPG2",
    "AGRN"
  ),

  ECM_Remodeling = c(
    "VCAN",
    "CSPG4",
    "MMP2",
    "MMP9",
    "TIMP1",
    "FN1",
    "COL1A1",
    "COL1A2"
  ),

  Endfoot_Organization = c(
    "AQP4",
    "SNTA1",
    "DAG1",
    "DTNA",
    "DMD"
  ),

  BBB_Specialization = c(
    "CLDN5",
    "OCLN",
    "TJP1",
    "MFSD2A",
    "SLC2A1",
    "ABCB1"
  )
)


# ============================================================
# 2. Pathway-scoring function
#
# Gene expression is standardized across specimens within
# each cohort. Genes with zero variance are excluded.
#
# The pathway score for each specimen is the mean z-score
# across all available genes assigned to that pathway.
# ============================================================

score_pathways <- function(expression_matrix,
                           gene_sets,
                           cohort_name = "cohort",
                           minimum_genes = 3) {

  expression_matrix <- as.matrix(expression_matrix)

  if (is.null(rownames(expression_matrix))) {
    stop("Expression matrix must have gene symbols as row names.")
  }

  score_list <- list()
  availability_list <- list()

  for (pathway in names(gene_sets)) {

    requested <- gene_sets[[pathway]]

    available <- intersect(
      requested,
      rownames(expression_matrix)
    )

    missing <- setdiff(
      requested,
      rownames(expression_matrix)
    )

    # Remove zero-variance genes
    variable_genes <- available[
      vapply(
        available,
        function(g) {
          x <- expression_matrix[g, ]
          isTRUE(sd(x, na.rm = TRUE) > 0)
        },
        logical(1)
      )
    ]

    if (length(variable_genes) < minimum_genes) {
      stop(
        paste(
          "Pathway", pathway,
          "has fewer than",
          minimum_genes,
          "usable genes in",
          cohort_name
        )
      )
    }

    pathway_expr <- expression_matrix[
      variable_genes,
      ,
      drop = FALSE
    ]

    # Standardize each gene across specimens
    pathway_z <- t(
      scale(
        t(pathway_expr),
        center = TRUE,
        scale = TRUE
      )
    )

    # Mean standardized expression per specimen
    score_list[[pathway]] <- colMeans(
      pathway_z,
      na.rm = TRUE
    )

    availability_list[[pathway]] <- data.frame(
      cohort = cohort_name,
      component = pathway,
      requested_n = length(requested),
      available_n = length(variable_genes),
      percent_available =
        100 * length(variable_genes) / length(requested),
      available_genes =
        paste(variable_genes, collapse = ";"),
      missing_genes =
        paste(
          union(
            missing,
            setdiff(available, variable_genes)
          ),
          collapse = ";"
        ),
      stringsAsFactors = FALSE
    )
  }

  scores <- data.frame(
    sample_id = colnames(expression_matrix),
    score_list,
    check.names = FALSE
  )

  availability <- do.call(
    rbind,
    availability_list
  )

  rownames(availability) <- NULL

  list(
    scores = scores,
    gene_availability = availability
  )
}


# ============================================================
# 3. Expected cohort-specific gene availability
# ============================================================
#
# TCGA:
#   56 of 56 curated genes available
#
# CGGA:
#   55 of 56 curated genes available
#   PECAM1 unavailable
#
# Pathway scoring was performed independently within each
# cohort. No cross-cohort normalization was performed.
#
# Final analytical cohorts used elsewhere in this repository:
#
# TCGA:
#   n = 800 primary-tumor specimens
#   Grade II  = 214
#   Grade III = 241
#   Grade IV  = 345
#
# CGGA:
#   n = 422 primary tumors
#   Grade II  = 138
#   Grade III = 144
#   Grade IV  = 140
#
# The unfiltered pathway-score files are provided in:
#
# processed_data/TCGA_pathway_scores_all_samples.csv
# processed_data/CGGA_pathway_scores_all_samples.csv
#
# The final analysis-ready files are provided in:
#
# processed_data/TCGA_analysis_cohort.csv
# processed_data/CGGA_analysis_cohort.csv
#
# ============================================================


# ============================================================
# 4. Example usage
# ============================================================
#
# TCGA:
#
# tcga_scored <- score_pathways(
#   expression_matrix = tcga_expression,
#   gene_sets = gene_sets,
#   cohort_name = "TCGA"
# )
#
# write.csv(
#   tcga_scored$scores,
#   "processed_data/TCGA_pathway_scores_reproduced.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   tcga_scored$gene_availability,
#   "pathway_definitions/gene_availability_TCGA_reproduced.csv",
#   row.names = FALSE
# )
#
#
# CGGA:
#
# cgga_scored <- score_pathways(
#   expression_matrix = cgga_expression,
#   gene_sets = gene_sets,
#   cohort_name = "CGGA"
# )
#
# write.csv(
#   cgga_scored$scores,
#   "processed_data/CGGA_pathway_scores_reproduced.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   cgga_scored$gene_availability,
#   "pathway_definitions/gene_availability_CGGA_reproduced.csv",
#   row.names = FALSE
# )
#
# ============================================================
