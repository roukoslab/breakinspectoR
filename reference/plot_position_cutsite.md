# plot_position_cutsite Plot the frequency of the position of the cutsite relative to the guide.

plot_position_cutsite Plot the frequency of the position of the cutsite
relative to the guide.

## Usage

``` r
plot_position_cutsite(x, guide, pam)
```

## Arguments

- x:

  A GRanges object containing the results of a breakinspectoR analysis.

- guide:

  character with the guide sequence. Only used for drawing.

- pam:

  character with PAM sequence. Only used for drawing.

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

plot_position_cutsite(offtargets, guide="GACCCCCTCCACCCCGCCTC", pam="NGG")
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the breakinspectoR package.
#>   Please report the issue to the authors.
```
