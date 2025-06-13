# cCRE Annotation Pipeline (Human hg19)

**Script:** `annotationScript_cre_human_hg19.R`

This R script:
1. Merges UCSC table annotations with a BED file of hg19 coordinates.
2. Builds Bioconductor GRanges objects for both UCSC cCRE regions and Illumina probes.
3. Finds overlaps between probes and cCRE regions.
4. Outputs a combined annotation table as CSV and RDS.

## How to run

```bash
Rscript annotationScript_cre_human_hg19.R

cat << 'EOF' > README.md
# cCRE Annotation Pipeline (Human hg19)

**Script:** annotationScript_cre_human_hg19.R

This R script:
1. Merges UCSC table annotations with a BED file of hg19 coordinates.
2. Builds Bioconductor GRanges objects for both UCSC cCRE regions and Illumina probes.
3. Finds overlaps between probes and cCRE regions.
4. Outputs a combined annotation table as CSV and RDS.

## How to run

Rscript annotationScript_cre_human_hg19.R

## Dependencies

Install in R:

install.packages("dplyr")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges","IRanges","readr"))

## Inputs & Outputs

- **Inputs:** UCSC table CSV, BED file, Illumina manifest  
- **Outputs:**  
  - final_merged_hg19_ucsc_file.txt & .csv  
  - Illumina_Human_hg19_Methylation_UCSC_cCRE_Overlaps.csv & .rds

