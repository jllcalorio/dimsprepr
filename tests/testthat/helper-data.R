# Shared test fixtures – lightweight, no external deps
# ponytail: one fixture set, every test file reuses it

make_test_data <- function(n_samples = 30, n_features = 20, n_qc = 6, seed = 7744) {
  set.seed(seed)
  n_total <- n_samples + n_qc

  x <- matrix(abs(rnorm(n_total * n_features, mean = 500, sd = 150)),
              nrow = n_total, ncol = n_features)
  colnames(x) <- paste0("F", seq_len(n_features))
  rownames(x) <- paste0("S", seq_len(n_total))

  # Sprinkle some zeros and NAs

  x[sample(length(x), 30)] <- 0
  x[sample(length(x), 15)] <- NA

  meta <- data.frame(
    Sample            = rownames(x),
    Group             = c(rep("Control", n_samples / 2),
                          rep("Treatment", n_samples / 2),
                          rep("QC", n_qc)),
    Batch             = as.integer(rep(1:2, length.out = n_total)),
    InjectionSequence = seq_len(n_total),
    SubjectID         = c(paste0("BIO", seq_len(n_samples)), rep(NA, n_qc)),
    stringsAsFactors  = FALSE
  )
  rownames(meta) <- meta$Sample

  list(x = x, metadata = meta)
}
