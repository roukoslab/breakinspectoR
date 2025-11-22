# Offtargets detected within the BreakTag test data included in this package

A GRanges object with the coordinates of the enriched regions, with
number of breaks detected, the cut site, the enrichment over the
non-targeted library and the p-value of the statistical test.

## Usage

``` r
data(breakinspectoR_examples)
```

## Format

A GRanges object with the coordinates of the enriched regions, as a
result from calling the breakinspectoR analysis with
\`extdata/vegfa.chr6.bed.gz\` and 'extdata/nontarget.chr6.bed.gz'.

## Details

It is the result of running \`breakinspectoR()\` with the example
datasets included in this package: \`extdata/vegfa.chr6.bed.gz\` and
'extdata/nontarget.chr6.bed.gz'
