# Generate a small simulated DIMS dataset for examples, tests, and vignettes
# Run this script from the package root with: source("data-raw/dims_example.R")
#
# Design: 60 injections in 2 batches of 30. Each batch opens and closes with a
# pooled QC and has a QC every 6th injection (6 QCs per batch, 12 in total),
# so QC-based drift correction (QCRSC) can run on every batch. Biological
# samples alternate Control/Treatment (24 each).

set.seed(4821)

n_batch    <- 2L
per_batch  <- 30L
n_samples  <- n_batch * per_batch
n_features <- 80L
qc_pos     <- c(1L, 7L, 13L, 19L, 25L, 30L)          # within-batch QC positions

batch    <- rep(seq_len(n_batch), each = per_batch)
pos      <- rep(seq_len(per_batch), times = n_batch)
is_qc    <- pos %in% qc_pos
group    <- character(n_samples)
group[is_qc]  <- "QC"
group[!is_qc] <- rep(c("Control", "Treatment"), length.out = sum(!is_qc))

# Feature-level abundance (log-normal) and per-sample dilution
feature_mean <- exp(rnorm(n_features, mean = 7, sd = 1.5))
dilution     <- exp(rnorm(n_samples, 0, 0.15))
dilution[is_qc] <- 1

# Analytical noise on every injection; biological variation on non-QC only
x <- outer(dilution, feature_mean) *
  matrix(exp(rnorm(n_samples * n_features, 0, 0.08)), n_samples, n_features)
x[!is_qc, ] <- x[!is_qc, ] * exp(rnorm(sum(!is_qc) * n_features, 0, 0.25))

# Group effect on first 10 features (Treatment up-regulated)
trt <- group == "Treatment"
x[trt, 1:10] <- x[trt, 1:10] * exp(rnorm(sum(trt) * 10, mean = 0.5, sd = 0.2))

# Batch effect (batch 2 shifted) and within-batch signal drift (decline)
x[batch == 2L, ] <- sweep(x[batch == 2L, ], 2, exp(rnorm(n_features, 0.25, 0.1)), "*")
drift_slope <- runif(n_features, 0.1, 0.4)
x <- x * (1 - outer((pos - 1) / (per_batch - 1), drift_slope))

# ~3% zeros (below detection) in biological samples, a few sparse features
bio_cells <- which(matrix(!is_qc, n_samples, n_features))
x[sample(bio_cells, round(0.03 * length(bio_cells)))] <- 0
x[, 78:80] <- 0
x[sample(n_samples, 3), 78] <- exp(rnorm(3, 7, 1))

dims_data <- round(x, 2)
colnames(dims_data) <- sprintf("mz_%04d", seq_len(n_features))
rownames(dims_data) <- sprintf("S%02d", seq_len(n_samples))

dims_metadata <- data.frame(
  Sample            = rownames(dims_data),
  Group             = group,
  Batch             = as.integer(batch),
  InjectionSequence = seq_len(n_samples),
  SubjectID         = ifelse(is_qc, NA_character_,
                             sprintf("BIO%02d", cumsum(!is_qc))),
  stringsAsFactors  = FALSE
)
rownames(dims_metadata) <- dims_metadata$Sample

save(dims_data,     file = "data/dims_data.rda",     compress = "xz", version = 2)
save(dims_metadata, file = "data/dims_metadata.rda", compress = "xz", version = 2)
