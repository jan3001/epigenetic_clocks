# Epigenetic Clocks

A single-script R pipeline to reproduce the humice methylation-aging analysis:

- **Script**: `epigenetic_clocks.R`  
- **Purpose**: merge β-values, filter by probe annotations, normalize, compute DNAm age & age acceleration, and plot results.

---

## Data

Place all inputs under a local `DATA/` directory:

### Public GEO datasets  
Download these from GEO and drop into `DATA/public/`:

- **GSE103010** → humanized-mouse blood (`humice`)  
- **GSE40799** → cord-blood samples batch 1 (`cordblood`)  
- **GSE30870** → cord-blood batch 2 & PBMC (89–103 yrs) (`cordblood`, `89_to_103_years`)  
- **GSE83334** → 5 year-old whole-blood (`5_years`)  
- **GSE65638** → 21–32 year whole-blood (`22_to_32_years`)  

### Lab-provided (private)  
Keep these under `DATA/lab/` (not public):

- **Samples.tsv**  
  Maps each sample ID to its chronological **Age** and **Condition** (`cordblood`, `humice`, `5_years`, `22_to_32_years`, `89_to_103_years`).

- **get450kProbeAnnotations_…Rda**  
  Contains vectors of probe IDs for autosomes, sex chromosomes, CpG-island/TSS regions, etc., used for targeted filtering.

---

## Requirements

* **R** (version 3.6 or higher)
* **R packages**:

  ```r
  install.packages(c("devtools", "preprocessCore", "ggplot2", "ggfortify"))
  devtools::install_github("yiluyucheng/dnaMethyAge")
  ```

## Usage

1. **Configure**
   Open `epigenetic_clocks.R` and set the path to your data folder:

   ```r
   data_dir <- "/full/path/to/DATA"
   ```

2. **Run**
   From the command line:

   ```bash
   Rscript epigenetic_clocks.R
   ```

## What it does

1. **Loads** public β-value matrices (adds an `ID_REF` column if missing)
2. **Merges** all datasets by CpG ID
3. **Filters** to autosomal probes using the lab’s annotation Rda
4. **Quantile-normalizes** the merged matrix and removes probes with missing data
5. **Reads** sample ages and conditions from `Samples.tsv`
6. **Computes** DNAm age (`mAge`) and age acceleration for published clocks:

   * HorvathS2013
   * HannumG2013
   * LevineM2018
   * ZhangQ2019
   * ShirebyG2020
   * YangZ2016
   * ZhangY2017
   * LuA2019
   * HorvathS2018
   * DunedinPACE
   * McEwenL2019
   * CBL\_specific
   * PCHorvathS2013, PCHannumG2013, PCHorvathS2018, PCPhenoAge
   * CBL\_common, Cortex\_common
   * epiTOC2
   * LuA2023p1, LuA2023p2, LuA2023p3
7. **Plots** scatterplots of epigenetic age vs. chronological age and boxplots comparing age acceleration across conditions

## Outputs

* Interactive **scatterplots** in your R session (one per clock)
* Example **boxplot** of age acceleration by condition




This repository contains only the R script to run epigenetic-age clocks:

- **epigenetic_clocks.R**: load, normalize, and compute DNAm ages.
