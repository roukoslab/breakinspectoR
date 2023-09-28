if(requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {

  offtargets <- breakinspectoR::breakinspectoR(
    target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
    nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
    guide    ="GACCCCCTCCACCCCGCCTC",
    PAM      =c("NGG", "NAG"),
    bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
    cutsiteFromPAM=3,
    verbose  =FALSE)

  test_that("Positive on-target is found in breakinspectoR analysis", {
    expect_is(offtargets, "GRanges")
    expect_true(offtargets[offtargets$mismatches == 0] == GRanges("chr6:43770822-43770841"))
  })

  offtargets_filtered <- breakinspectoR::reduce(offtargets)

  test_that("Positive on-target is found after reducing the list of offtargets", {
    expect_is(offtargets_filtered, "GRanges")
    expect_true(length(offtargets_filtered) <= length(offtargets))
    expect_true(offtargets_filtered[offtargets_filtered$mismatches == 0] == GRanges("chr6:43770822-43770841"))
  })
}

