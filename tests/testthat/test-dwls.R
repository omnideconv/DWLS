sc_object_small <- as.matrix(utils::read.csv("small_test_data/sc_object_small.csv", row.names = 1))
cell_annotations_small <- readr::read_lines("small_test_data/cell_annotations_small.txt")



test_that("DWLS build signature MAST works", {
  signature <- buildSignatureMatrixMAST(sc_object_small, cell_annotations_small, tempdir(), TRUE, NULL)
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single cell matrix",
    object = ncol(signature), expected = length(unique(cell_annotations_small))
  )
  check_signature <- as.matrix(read.csv("test_models/dwls_model_small.csv",
    row.names = 1,
    check.names = FALSE
  ))
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})

test_that("DWLS build signature MAST optimized works", {
  signature <- buildSignatureMatrixMASTOptimized(sc_object_small, cell_annotations_small, tempdir(), TRUE, NULL)
  expect_equal(
    info = "signature matrix has same amount of columns as unique cell types in single cell matrix",
    object = ncol(signature), expected = length(unique(cell_annotations_small))
  )
  check_signature <- as.matrix(read.csv("test_models/dwls_model_small.csv",
                                        row.names = 1,
                                        check.names = FALSE
  ))
  expect_equal(info = "signature matrix is correct", object = signature, expected = check_signature)
})



