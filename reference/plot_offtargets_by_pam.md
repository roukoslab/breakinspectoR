# plot_offtargets_by_pam Plot the count of offtargets by the different PAMs

plot_offtargets_by_pam Plot the count of offtargets by the different
PAMs

## Usage

``` r
plot_offtargets_by_pam(x, fraction = TRUE)
```

## Arguments

- x:

  A GRanges object containing the results of a breakinspectoR analysis.

- fraction:

  logical indicating whether the plot should show absolute absolute
  numbers or the fraction.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
  offtargets <- breakinspectoR(
    target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
    nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
    guide    ="GACCCCCTCCACCCCGCCTC",
    PAM      =c("NGG", "NAG"),
    bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
    cutsiteFromPAM=3
  )
} # }
data(breakinspectoR_examples)

plot_offtargets_by_pam(offtargets, fraction=TRUE)
```
