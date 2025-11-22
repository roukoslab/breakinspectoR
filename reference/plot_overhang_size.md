# plot_overhang_size Plots the proportion of blunt and staggered reads by size.

plot_overhang_size Plots the proportion of blunt and staggered reads by
size.

## Usage

``` r
plot_overhang_size(x)
```

## Arguments

- x:

  A named list of GRanges object containing the results of a
  scission_profile_analysis.

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

  offtargets.scission_profile <- scission_profile_analysis(
    x        =offtargets,
    target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
    nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
    bsgenome ="BSgenome.Hsapiens.UCSC.hg38")
} # }
data(breakinspectoR_examples)

# simulate 2 different experiments by picking 25 random offtargets
spcas9 <- sample(offtargets.scission_profile, 25)
lz3    <- sample(offtargets.scission_profile, 25)
plot_overhang_size(list(SpCas9=spcas9, LZ3=lz3))
#> Using L2 as id variables
#> Using L2 as id variables
```
