#' Offtargets detected within the BreakTag test data included in this package
#'
#' A GRanges object with the coordinates of the enriched regions, with
#' number of breaks detected, the cut site, the enrichment over the
#' non-targeted library and the p-value of the statistical test.
#'
#' It is the result of running `breakinspectoR()` with the example datasets
#' included in this package:
#' `extdata/vegfa.chr6.bed.gz` and 'extdata/nontarget.chr6.bed.gz'
#'
#' @docType data
#' @keywords datasets
#' @name offtargets
#' @usage data(breakinspectoR_examples)
#' @format A GRanges object with the coordinates of the enriched regions,
#' as a result from calling the breakinspectoR analysis with
#' `extdata/vegfa.chr6.bed.gz` and 'extdata/nontarget.chr6.bed.gz'.
NULL

#' Scission profile of the offtargets detected within the BreakTag test data
#' included in this package
#'
#' A GRanges object with the coordinates of the cutsite, the number of
#' breaks on and around the cutsite for the target and non-target libraries,
#' and a p- and q-values.
#'
#' It is the result of running `scission_profile_analysis()` with the example
#' datasets included in this package:
#' `extdata/vegfa.chr6.bed.gz` and 'extdata/nontarget.chr6.bed.gz'
#'
#' @docType data
#' @keywords datasets
#' @name offtargets.scission_profile
#' @usage data(breakinspectoR_examples)
#' @format A GRanges object with the coordinates of the cutsite, the number of
#' breaks on and around the cutsite for the target and non-target libraries,
#' and a p- and q-values, as a result from callilng the breakinspectoR analysis
#' with `extdata/vegfa.chr6.bed.gz` and 'extdata/nontarget.chr6.bed.gz'.
NULL
