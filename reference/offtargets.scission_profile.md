# Scission profile of the offtargets detected within the BreakTag test data included in this package

A GRanges object with the coordinates of the cutsite, the number of
breaks on and around the cutsite for the target and non-target
libraries, and a p- and q-values.

## Usage

``` r
data(breakinspectoR_examples)
```

## Format

A GRanges object with the coordinates of the cutsite, the number of
breaks on and around the cutsite for the target and non-target
libraries, and a p- and q-values, as a result from callilng the
breakinspectoR analysis with \`extdata/vegfa.chr6.bed.gz\` and
'extdata/nontarget.chr6.bed.gz'.

## Details

It is the result of running \`scission_profile_analysis()\` with the
example datasets included in this package: \`extdata/vegfa.chr6.bed.gz\`
and 'extdata/nontarget.chr6.bed.gz'
