# Data Dictionary

This document describes the processed analytical datasets distributed with this repository.

## TCGA_analysis_cohort.csv

Processed TCGA diffuse glioma analytical cohort used for pathway-level grade and IDH analyses.

| Variable | Description |
|---|---|
| Sample_ID | TCGA specimen identifier |
| Patient_ID | TCGA patient identifier |
| Sample_Type_Code | TCGA sample-type code; the analytical cohort contains primary tumor specimens (code 01) |
| Grade | Original histopathological WHO grade (II, III, or IV) |
| Grade_numeric | Numeric representation of WHO grade |
| IDH_status | Recorded IDH mutation status |
| IDH_Group | IDH grouping used for statistical analyses |
| Homeostatic_Astrocyte | Standardized homeostatic astrocyte pathway score |
| Reactive_Astrocyte | Standardized reactive astrocyte pathway score |
| Neuroinflammation | Standardized neuroinflammatory pathway score |
| Endothelial_Angiogenesis | Standardized endothelial/angiogenesis pathway score |
| Basement_Membrane | Standardized basement membrane pathway score |
| Extracellular_Matrix | Standardized extracellular matrix pathway score |
| Astrocyte_Endfoot | Standardized astrocyte endfoot pathway score |
| BBB_Specialization | Standardized blood-brain barrier specialization pathway score |

## TCGA_survival_cohort.csv

TCGA analytical dataset used for Cox proportional-hazards analyses.

| Variable | Description |
|---|---|
| Sample_ID | TCGA specimen identifier |
| Grade | Original histopathological WHO grade |
| IDH_Group | IDH grouping used in survival models |
| OS_days | Overall survival time in days |
| OS_status | Overall survival event indicator |
| Pathway score columns | Standardized scores for the eight gliovascular transcriptional programs |

## CGGA_analysis_cohort.csv

Processed primary-tumor cohort from the CGGA mRNAseq_693 dataset.

| Variable | Description |
|---|---|
| CGGA_ID | CGGA specimen identifier |
| Grade | Original histopathological WHO grade (II, III, or IV) |
| Grade_numeric | Numeric representation of WHO grade |
| IDH_mutation_status | Recorded IDH mutation status |
| IDH_Group | IDH grouping used for statistical analyses |
| Age | Age at diagnosis |
| OS | Overall survival time |
| Censor (alive=0; dead=1) | Overall survival event indicator |
| Pathway score columns | Standardized scores for the eight gliovascular transcriptional programs |

## IvyGAP_sentinel_gene_expression.csv

Processed Ivy Glioblastoma Atlas Project data used for regional sentinel-gene analyses.

| Variable | Description |
|---|---|
| gene_symbol | Gene symbol for one of the nine sentinel genes |
| rna_well_id | Ivy GAP RNA sample identifier |
| tumor_id | Tumor/donor identifier used for repeated-measures modeling |
| Region | Anatomical region: LE, IT, CT, MVP, or PAN |
| FPKM | FPKM expression value |
| expression | log2(FPKM + 1)-transformed expression |

### Ivy GAP region abbreviations

- **LE** — Leading edge
- **IT** — Infiltrating tumor
- **CT** — Cellular tumor
- **MVP** — Microvascular proliferation
- **PAN** — Pseudopalisading cells around necrosis

## Pathway scores

Pathway scores were calculated independently within each cohort. For each pathway, available marker genes were standardized across specimens using z-scores, and the pathway score was calculated as the mean standardized expression of the available genes.

The TCGA analysis contained all 56 curated genes. The CGGA analysis contained 55 of 56 genes because PECAM1 was unavailable. No cross-cohort normalization of pathway scores was performed.

See `pathway_definitions/gliovascular_pathway_definitions.csv` for the complete gene definitions.
