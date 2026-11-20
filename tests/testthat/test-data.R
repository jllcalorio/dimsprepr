test_that("bundled datasets are consistent with each other", {
  expect_true(is.matrix(dims_data) && is.numeric(dims_data))
  expect_equal(dim(dims_data), c(60L, 80L))
  expect_identical(rownames(dims_data), dims_metadata$Sample)
  expect_named(dims_metadata, c("Sample", "Group", "Batch", "InjectionSequence", "SubjectID"))
  expect_false(anyNA(dims_data))
  expect_true(all(dims_data >= 0))
})

test_that("every batch starts and ends with a QC and has >= 5 QCs (QCRSC-ready)", {
  for (b in split(dims_metadata, dims_metadata$Batch)) {
    b <- b[order(b$InjectionSequence), ]
    expect_equal(b$Group[c(1, nrow(b))], c("QC", "QC"))
    expect_gte(sum(b$Group == "QC"), 5L)
  }
  expect_true(all(is.na(dims_metadata$SubjectID[dims_metadata$Group == "QC"])))
})
