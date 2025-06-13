# Load necessary libraries
library(dplyr)

# Read the BED file (hglft_genome)
bed_file <- read.table("/Users/anjan/Library/CloudStorage/OneDrive-JohnsHopkins/Easwaran Lab/Aging_NIH/Annotation Datasets/hglft_genome_217ad_6ae950.bed", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
bed_v2_yuba <- read.table("/Users/anjan/Library/CloudStorage/OneDrive-JohnsHopkins/Easwaran Lab/Aging_NIH/Annotation Datasets/yuba_encodeCcreCombined_4column.bed", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
  


# Assign column names to the BED file
colnames(bed_file) <- c("chrom_hg19", "start_hg19", "end_hg19", "name", "score", "strand")

# Read the UCSC human table file
ucsc_file <- read.csv("/Users/anjan/Library/CloudStorage/OneDrive-JohnsHopkins/Easwaran Lab/Aging_NIH/Annotation Datasets/ucsc_humantable_nonBED", header = TRUE, stringsAsFactors = FALSE)


# Merge the two files based on the 'name' column
merged_data <- merge(ucsc_file, bed_file[, c("name", "chrom_hg19", "start_hg19", "end_hg19")], by = "name", all.x = TRUE)

# Replace the old coordinate columns (chromStart, chromEnd) with the new coordinates from hg19 (start_hg19, end_hg19)
merged_data$chromStart <- merged_data$start_hg19
merged_data$chromEnd <- merged_data$end_hg19

# Drop the 'chrom', 'start', and 'end' columns (from UCSC)
merged_data <- merged_data %>%
  select(-chrom, -start, -end)

# Drop hg19 temporary columns
merged_data <- merged_data %>%
  select(-start_hg19, -end_hg19, -chrom_hg19)

merged_data_clean <- merged_data %>%
  filter(!is.na(chromStart) & !is.na(chromEnd))

# Write the final merged data to a new file
write.table(merged_data_clean, "/Users/anjan/Library/CloudStorage/OneDrive-JohnsHopkins/Easwaran Lab/Aging_NIH/Annotation Datasets/final_merged_hg19_ucsc_file.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(merged_data_clean, "/Users/anjan/Library/CloudStorage/OneDrive-JohnsHopkins/Easwaran Lab/Aging_NIH/Annotation Datasets/final_merged_hg19_ucsc_file.csv", sep = "\t", row.names = FALSE, quote = FALSE)

# View the first few rows to confirm
head(merged_data)



###NOW CREATING ANNOTATION SCRIPT
# Loading required libraries
library(GenomicRanges)
library(readr)

# Loading Illumina methylation manifest (skipping first 7 rows for metadata)
# if running on windows run this
illumina_df <- read_csv("C:/Users/anjan/OneDrive - Johns Hopkins/Easwaran Lab/Aging_NIH/Annotation Datasets/humanmethylation450_15017482_v1-2.csv", skip = 7)
# if running on macbook, run this
illumina_df <- read_csv("/Users/anjan/Library/CloudStorage/OneDrive-JohnsHopkins/Easwaran Lab/Aging_NIH/Annotation Datasets/humanmethylation450_15017482_v1-2.csv", skip =7)

# Filtering out rows where MAPINFO or CHR is NA (required for GRanges)
illumina_df_filtered <- illumina_df[!is.na(illumina_df$MAPINFO) & !is.na(illumina_df$CHR), ]

# Filtering out non-autosomal chromosomes (we only want autosomes, so exclude X and Y)
illumina_df_filtered <- illumina_df_filtered[!illumina_df_filtered$CHR %in% c("X", "Y"), ]

# Creating GRanges object for Illumina data (adjust Illumina 1-based index to match to UCSC's 0-based)
# Treating the MAPINFO as a range: start = MAPINFO - 1 (convert to 0-based), end = MAPINFO
illumina_gr <- GRanges(
  seqnames = Rle(paste0("chr", illumina_df_filtered$CHR)),  # Prefix with 'chr' to match UCSC format
  ranges = IRanges(start = illumina_df_filtered$MAPINFO - 1, end = illumina_df_filtered$MAPINFO),  # Adjusting start for 0-based compatibility
  Illumina_ID = illumina_df_filtered$IlmnID  # Storing Illumina probe ID as metadata
)

# Checking the resulting GRanges object for Illumina data
illumina_gr


#now the UCSC file
ucsc_df <- merged_data_clean

# View the first few rows to confirm
head(ucsc_df)

# Filtering out non-autosomal chromosomes in UCSC data (we exclude X, Y, and others like 'chrM')
# ucsc_df_filtered <- ucsc_df[!ucsc_df$X.chrom %in% c("chrX", "chrY", "chrM"), ]
ucsc_df_filtered <- ucsc_df[grep("^(chr[0-9]+)$", ucsc_df$X.chrom), ]

# Creating GRanges object for UCSC data (leave as is, as Illumina file has been modified to account for UCSC 0-based indexing)
ucsc_gr <- GRanges(
  seqnames = Rle(ucsc_df_filtered$X.chrom),  # Chromosome names
  ranges = IRanges(start = ucsc_df_filtered$chromStart, end = ucsc_df_filtered$chromEnd),  # UCSC uses 0-based start, 1-based end
  region_type = ucsc_df_filtered$encodeLabel  # Region type (e.g., "promoter-like", "enhancer signature")
)

# Finding overlaps between Illumina and UCSC regions with specific parameters
overlaps <- findOverlaps(
  query = illumina_gr,
  subject = ucsc_gr,
  maxgap = 2L,           # No gaps allowed between ranges
  minoverlap = 0L,       # At least 1 base overlap
  type = "any",          # Overlaps of any type are allowed
  select = "all",        # Return all overlaps
  ignore.strand = TRUE   # Ignore strand information (not important in this case)
)

# Convert GRanges objects to data frames for easier binding
illumina_hits_df <- as.data.frame(illumina_gr[queryHits(overlaps)])  # Illumina probe data
ucsc_hits_df <- as.data.frame(ucsc_gr[subjectHits(overlaps)])  # UCSC enhancer regions data

# Combine metadata from Illumina and UCSC for further analysis
# This step binds the results from the two GRanges objects into a combined dataframe
combined_overlap <- cbind(illumina_hits_df, ucsc_hits_df)

# View the resulting combined overlap data
head(combined_overlap)

# Alternatively, you can also combine the most relevant columns from both datasets if needed:
combined_data <- data.frame(
  Illumina_ID = illumina_hits_df$Illumina_ID,
  Illumina_Chrom = illumina_hits_df$seqnames,  # Illumina chromosome info
  Illumina_Start = illumina_hits_df$start,  # Illumina start position
  Illumina_End = illumina_hits_df$end,  # Illumina end position
  UCSC_Chrom = ucsc_hits_df$seqnames,  # UCSC chromosome info
  UCSC_Start = ucsc_hits_df$start,  # UCSC start position
  UCSC_End = ucsc_hits_df$end,  # UCSC end position
  UCSC_region_type = ucsc_hits_df$region_type  # UCSC region type
)

# View the combined metadata with selected columns
head(combined_data)
length(unique(combined_data$Illumina_ID))

#Saving the files
file_path <- "/Users/anjan/Library/CloudStorage/OneDrive-JohnsHopkins/Easwaran Lab/Aging_NIH/Annotation Datasets"
write.csv(combined_data, paste0(file_path, "/Illumina_Human_hg19_Methylation_UCSC_cCRE_Overlaps.csv"), row.names = FALSE)

# Saving the combined data as an RDS file with a descriptive name
saveRDS(combined_data, paste0(file_path, "/Illumina_Human_hg19_Methylation_UCSC_cCRE_Overlaps.rds"))
