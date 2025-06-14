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

- R (≥ 3.6)  
- Packages:  
  ```r
  install.packages(c("devtools","preprocessCore","ggplot2","ggfortify"))
  devtools::install_github("yiluyucheng/dnaMethyAge")
````

---

## Usage

1. **Configure**
   At the top of `epigenetic_clocks.R`, set:

   ```r
   data_dir <- "/full/path/to/DATA"
   ```

2. **Run**

   ```bash
   Rscript epigenetic_clocks.R
   ```

---

## What it does

1. **Loads** public β-matrices (adds an `ID_REF` column if needed)
2. **Merges** all datasets by CpG ID
3. **Filters** to autosomal probes via the lab annotation `.Rda`
4. **Quantile-normalizes** and drops incomplete probes
5. **Reads** sample ages & conditions from `Samples.tsv`
6. **Computes** epigenetic age (mAge) & age acceleration for published clocks (HorvathS2013, HannumG2013, LevineM2018, ZhangQ2019, ShirebyG2020, YangZ2016, ZhangY2017, LuA2019, HorvathS2018, DunedinPACE, McEwenL2019, CBL\_specific, PC‐ and CBL‐ variants, epiTOC2, LuA2023p1–p3)
7. **Plots** scatterplots of mAge vs. chronological age and boxplots by condition

---

## Outputs

* Scatterplots appear in your R session (one per clock)
* Example boxplot comparing mAge distributions across conditions

---

## License

MIT © Your Name

```
```


# Epigenetic Clocks

This repository contains only the R script to run epigenetic-age clocks:

- **epigenetic_clocks.R**: load, normalize, and compute DNAm ages.
