# Test helper functions for conditional skipping

skip_if_no_python <- function() {
  if (!reticulate::py_available(initialize = FALSE)) {
    testthat::skip("Python not available")
  }
}

skip_if_no_anndata <- function() {
  # Check basilisk.utils compatibility first (required by anndata)
  # Do NOT load anndata - check basilisk.utils independently
  if (!requireNamespace("basilisk.utils", quietly = TRUE)) {
    testthat::skip("basilisk.utils not installed")
  }
  # Check the function exists before loading anndata
  if (!"activateEnvironment" %in% getNamespaceExports("basilisk.utils")) {
    testthat::skip("basilisk.utils activateEnvironment not available (version mismatch)")
  }
  # Now safe to check if anndata R package is installed
  if (!requireNamespace("anndata", quietly = TRUE)) {
    testthat::skip("anndata R package not installed")
  }
  # Also check Python anndata module
  skip_if_no_python()
  if (!reticulate::py_module_available("anndata")) {
    testthat::skip("anndata Python module not available")
  }
}

skip_if_no_seuratdisk <- function() {
  if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    testthat::skip("SeuratDisk not installed")
  }
}

skip_if_no_seuratwrappers <- function() {
  if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    testthat::skip("SeuratWrappers not installed")
  }
}

skip_if_no_monocle3 <- function() {
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    testthat::skip("monocle3 not installed")
  }
}

skip_if_no_loomr <- function() {
  if (!requireNamespace("loomR", quietly = TRUE)) {
    testthat::skip("loomR not installed")
  }
}
