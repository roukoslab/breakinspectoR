# plot_blunt_rate_density Density Plot for Blunt Rate Distribution.

plot_blunt_rate_density Density Plot for Blunt Rate Distribution.

## Usage

``` r
plot_blunt_rate_density(x)
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
plot_blunt_rate_density(list(SpCas9=spcas9, LZ3=lz3))
#> Picking joint bandwidth of 0.506
```
