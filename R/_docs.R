#' @title anndata_to functions
#'
#' @description
#' Functions to convert a specific object type to another specific object type.
#'
#' NOTE: Do NOT try to \code{inheritDotParams sceasy::convertFormat}.
#' The documentation in that function is empty and this will a
#' cryptic error with \link[devtools]{document}.
#' @family anndata_to
#' @param standardise Run \link[EWCE]{standardise_ctd}.
#' @inheritParams converters
#' @inheritParams map_data
#' @inheritParams to_se
#' @inheritParams orthogene::aggregate_rows
#' @inheritParams EWCE::generate_celltype_data
#' @inheritParams EWCE::standardise_ctd
#' @returns Converted object.
#'
#' @name anndata_to
#' @import data.table
#' @import orthogene
#' @importFrom utils getFromNamespace
#' @importFrom EWCE generate_celltype_data standardise_ctd
NULL

#' @title calc_ functions
#' @family calc_
#' @name calc_
#'
#' @description
#' Calculations functions that operate on matrices.
#' @returns Calculation results.
NULL

#' @title to_ functions
#'
#' @description
#' Functions to convert one any object type to another specific object type.
#' @family to_
#' @returns Converted object.
#'
#' @name to_
#' @import data.table
NULL
