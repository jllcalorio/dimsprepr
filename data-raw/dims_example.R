# Generate a small simulated DIMS dataset for examples, tests, and vignettes
# Run this script with: source("data-raw/dims_example.R")

set.seed(4821)

n_bio     <- 40L   # 20 Control + 20 Treatment
n_qc      <- 10L
n_samples <- n_bio + n_qc
n_features <- 80L

# Simulate raw intensities (log-normal-ish)
base_intensities <- matrix(
  exp(rnorm(n_samples * n_features, mean = 7, sd = 1.5)),
  nrow = n_samples, ncol = n_features
)
colnames(base_intensities) <- paste0("mz_", sprintf("%04d", seq_len(n_features)))

# Add mild group effect to first 10 features for Treatment
treatment_rows <- 21:40
base_intensities[treatment_rows, 1:10] <-
  base_intensities[treatment_rows, 1:10] * exp(rnorm(20 * 10, mean = 0.5, sd = 0.3))

# Add mild batch effect
batch_vec <- rep(1:2, each = n_samples / 2)
base_intensities[batch_vec == 2, ] <-
  base_intensities[batch_vec == 2, ] * exp(rnorm(n_features, mean = 0.2, sd = 0.1))

# Sprinkle zeros (below-LOD) — ~5% of values
zero_mask <- sample(length(base_intensities), size = round(0.05 * length(base_intensities)))
base_intensities[zero_mask] <- 0

# Make a few features almost entirely missing (to exercise missing-value filter)
base_intensities[, 78:80] <- 0
base_intensities[sample(n_samples, 3), 78] <- exp(rnorm(3, 7, 1))

# Round to mimic instrument output
dims_data <- round(base_intensities, 2)
rownames(dims_data) <- paste0("S", sprintf("%02d", seq_len(n_samples)))

# Metadata
dims_metadata <- data.frame(
  Sample            = rownames(dims_data),
  Group             = c(rep("Control", 20), rep("Treatment", 20), rep("QC", 10)),
  Batch             = as.integer(batch_vec),
  InjectionSequence = seq_len(n_samples),
  SubjectID         = c(paste0("BIO", sprintf("%02d", 1:40)), rep(NA_character_, 10)),
  stringsAsFactors  = FALSE
)
rownames(dims_metadata) <- dims_metadata$Sample

usethis::use_data(dims_data, dims_metadata, overwrite = TRUE)
