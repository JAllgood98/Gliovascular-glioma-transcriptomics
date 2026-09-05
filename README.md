# Gliovascular-glioma-transcriptomics

Analysis code and processed data accompanying the manuscript:

**Integrated Transcriptomic Analysis Reveals Coordinated Gliovascular Remodeling Across Diffuse Glioma Progression**

## Overview

This repository contains processed analytical data, curated pathway definitions, statistical results, and analysis code used to characterize gliovascular-associated transcriptional remodeling across diffuse glioma.

The study integrates two independent bulk RNA-sequencing cohorts with anatomically annotated glioblastoma data:

- The Cancer Genome Atlas (TCGA-LGG and TCGA-GBM)
- Chinese Glioma Genome Atlas (CGGA; mRNAseq_693)
- Ivy Glioblastoma Atlas Project (Ivy GAP)

Eight curated transcriptional programs representing major components of gliovascular organization and remodeling were evaluated:

1. Endfoot organization
2. Homeostatic astrocytes
3. Reactive astrocytes
4. BBB specialization
5. Endothelial angiogenesis
6. Basement membrane
7. Extracellular matrix remodeling
8. Neuroinflammation

## Cohorts

### TCGA

RNA-sequencing and clinical data were obtained from the TCGA-LGG and TCGA-GBM projects through the Genomic Data Commons using TCGAbiolinks. Gene-expression data were derived from STAR-Counts FPKM-UQ values.

The final primary-tumor analytical cohort contained 800 specimens:

- WHO grade II: n = 214
- WHO grade III: n = 241
- WHO grade IV: n = 345

Analyses were performed at the specimen level. Multiple primary-tumor RNA aliquots were retained when available for the same patient.

### CGGA

The independent CGGA mRNAseq_693 cohort was used for replication. Following restriction to primary tumors and availability of required clinical variables, the analytical cohort contained 422 tumors:

- WHO grade II: n = 138
- WHO grade III: n = 144
- WHO grade IV: n = 140

### Ivy GAP

Anatomically annotated glioblastoma RNA-sequencing data from the Ivy Glioblastoma Atlas Project were used to evaluate regional variation in gliovascular-associated transcriptional programs.

## Pathway scoring

The eight gliovascular-associated programs were defined using a curated 56-gene panel.

Within each cohort, expression of each available marker gene was standardized as a z-score, and pathway activity was calculated as the mean standardized expression of the genes assigned to that pathway.

All 56 genes were available in TCGA. Fifty-five were available in CGGA because PECAM1 was unavailable in the CGGA expression matrix.

Pathway scoring was performed independently within each cohort; no cross-cohort normalization was performed.

## Repository structure

```text
pathway_definitions/
    Curated 56-gene pathway definitions and cohort-specific gene availability.

processed_data/
    Processed pathway scores and analysis-ready cohort data used in the study.

results/
    Statistical outputs, including grade-comparison bootstrap analyses,
    gene-level Loss/Gain analyses, and proportional-hazards diagnostics.

scripts/
    Analysis and figure-generation scripts used to reproduce the reported analyses.

figures/
    Figure-related outputs where applicable.
