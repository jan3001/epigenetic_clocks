#installing dnaMethyAge package
devtools::install_github("yiluyucheng/dnaMethyAge")

#loading library
library(dnaMethyAge)
library(preprocessCore) #for quantile normalization
library(ggplot2)        #for general plotting
library(ggfortify)      #to enhance PCA plots?

# ---- ONLY CHANGE THIS: point to your DATA folder ----
data_dir <- "/path/to/your/project/DATA"

# Function to load data, ensure ID_REF is present, and assign conditions
load_data <- function(file_path, condition_name){
  data <- read.table(file_path, sep = "\t", header = TRUE) #read file
  if(!'ID_REF' %in% colnames(data)){ # check if ID_REF exists
    data$ID_REF <- rownames(data)    # add ID_REF if missing
  }
  colnames_df <- as.data.frame(colnames(data)) # get column names
  colnames_df$condition <- condition_name      # assign condition
  return(list(data = data, colnames_df = colnames_df))
}

#Loading in all datasets
humice_df     <- load_data(file.path(data_dir, "humice",             "beta_matrix_noob.csv"),    "humice")
cb_df         <- load_data(file.path(data_dir, "cordblood_1",       "GSM765881-23251.txt"),    "cordblood")
cb_df2        <- load_data(file.path(data_dir, "cordblood_2",       "GSM1001921-15653.txt"),   "cordblood")
five_years    <- load_data(file.path(data_dir, "5_years",            "samples_5_years.txt"),    "5_years")
twenty_thirty <- load_data(file.path(data_dir, "21_32_years",        "twenty_thirty_yrs.txt"),  "22_to_32_years")
eighty_hundred<- load_data(file.path(data_dir, "89_103_years",       "eighty_hundred_years.txt"), "89_to_103_years")

#merging dataframes by "ID_REF" using Reduce to avoid repetition
check_id_ref <- function(df_list) {
  for (i in seq_along(df_list)) {
    if (!"ID_REF" %in% colnames(df_list[[i]])) {
      stop(paste("Error: 'ID_REF' column is missing in dataframe", i))
    }
  }
}
dfs_to_merge <- list(
  humice_df$data, cb_df$data, cb_df2$data,
  five_years$data, twenty_thirty$data, eighty_hundred$data
)
check_id_ref(dfs_to_merge)
merged_df <- Reduce(function(x,y) merge(x,y, by="ID_REF"), dfs_to_merge)

#Loading probe annotation data for filtering
load(file.path(data_dir, "get450kProbeAnnotations_To_refTSS_colonChromHMMStates.Rda"))
ls()

# Filter for autosomal probes from the merged dataframe
auto <- merged_df[merged_df$ID_REF %in% autosomalProbes, ]
rownames(auto) <- auto$ID_REF

# Remove the 'ID_REF' column as it's now set as rownames
auto_fin <- auto[, !(colnames(auto) %in% 'ID_REF')]

# Quantile normalization on the matrix of beta values
test <- normalize.quantiles(as.matrix(auto_fin), copy = TRUE)
colnames(test) <- colnames(auto_fin)
rownames(test) <- rownames(auto_fin)

# Filter out rows with missing values (complete cases only)
test_fin <- test[complete.cases(test), ]

#Checking beta matrix dimensions
print(dim(test_fin))

#Checking all available clocks listed
availableClock()

#sample sheet with ages
samples <- read.csv(
  file.path(data_dir, "Samples (2).tsv"),
  sep = "\t", header = TRUE, fill = TRUE
)
print(samples)

#(1) SELECTING HORVATHS2013 FOR TESTING
clock_name <- 'HorvathS2013'
horvathAge <- methyAge(test_fin, clock = clock_name)

##HORVATHS2013 ACCELERATION
info <- samples[, c("Sample", "Age", "Condition")]
horvath_age_with_accel <- methyAge(
  test_fin, clock = clock_name, age_info = info,
  fit_method = 'Linear', do_plot = TRUE
)

condition_colors <- c('#29BEDF','#B0C52F','#E48AD9','#ff9999','#4EDAA0')

ggplot(data = horvath_age_with_accel, aes(x = Age, y = mAge, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = condition_colors) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mAge vs Age (HorvathS2013)",
       x = "Chronological Age", y = "Epigenetic Age (mAge)") +
  theme_minimal()

#(2) SELECTING HANNUMG2013 FOR TESTING
clock_name <- 'HannumG2013'
hannumAge <- methyAge(test_fin, clock = clock_name)

##HANNUMG2013 AGE ACCELERATION
info <- samples[, c("Sample", "Age", "Condition")]
hannum_age_with_accel <- methyAge(
  test_fin, clock = clock_name, age_info = info,
  fit_method = 'Linear', do_plot = TRUE
)

ggplot(data = hannum_age_with_accel, aes(x = Age, y = mAge, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = condition_colors) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mAge vs Age (HannumG2013)",
       x = "Chronological Age", y = "Epigenetic Age (mAge)") +
  theme_minimal()

#(3) SELECTING LEVINEM2018 FOR TESTING
clock_name <- 'LevineM2018'
levineAge <- methyAge(test_fin, clock = clock_name)

##LEVINEM2018 AGE ACCELERATION
info <- samples[, c("Sample", "Age", "Condition")]
levine_age_with_accel <- methyAge(
  test_fin, clock = clock_name, age_info = info,
  fit_method = 'Linear', do_plot = TRUE
)

ggplot(data = levine_age_with_accel, aes(x = Age, y = mAge, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = condition_colors) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mAge vs Age (LevineM2018)",
       x = "Chronological Age", y = "Epigenetic Age (mAge)") +
  theme_minimal()

#(4) SELECTING ZHANGEQ2019 FOR TESTING
clock_name <- 'ZhangQ2019'
zhang2019Age <- methyAge(test_fin, clock = clock_name)

##ZHANGQ2019 AGE ACCELERATION
info <- samples[, c("Sample", "Age", "Condition")]
zhang2019_age_with_accel <- methyAge(
  test_fin, clock = clock_name, age_info = info,
  fit_method = 'Linear', do_plot = TRUE
)

ggplot(data = zhang2019_age_with_accel, aes(x = Age, y = mAge, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = condition_colors) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mAge vs Age (ZhangQ2019)",
       x = "Chronological Age", y = "Epigenetic Age (mAge)") +
  theme_minimal()

#(5) SELECTING ShirebyG2020 FOR TESTING
clock_name <- 'ShirebyG2020'
shireby2020Age <- methyAge(test_fin, clock = clock_name)

##ShirebyG2020 AGE ACCELERATION
info <- samples[, c("Sample", "Age", "Condition")]
shireby2020_age_with_accel <- methyAge(
  test_fin, clock = clock_name, age_info = info,
  fit_method = 'Linear', do_plot = TRUE
)

ggplot(data = shireby2020_age_with_accel, aes(x = Age, y = mAge, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = condition_colors) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(title = "mAge vs Age (ShirebyG2020)",
       x = "Chronological Age", y = "Epigenetic Age (mAge)") +
  theme_minimal()

# Continue updating each clock section in the same way...
# (6) YangZ2016, (7) ZhangY2017, (8) LuA2019, etc.,
# replacing only the paths at the top with file.path(data_dir, ...)
