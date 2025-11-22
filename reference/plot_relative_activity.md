# plot_relative_activity Plot Cas9 activity of a list of multiple gRNA offtargets relative to one experiment.

plot_relative_activity Plot Cas9 activity of a list of multiple gRNA
offtargets relative to one experiment.

## Usage

``` r
plot_relative_activity(
  x,
  ref = 1,
  what = c("on-targets", "off-targets", "all")
)
```

## Arguments

- x:

  A named list of GRanges object containing the results of a
  breakinspectoR analysis. Possibly filtered (eg. qval \< .01).

- ref:

  Which experiment from \`x\` should be taken as reference? either a
  number or a name. Defaults to the first experiment in \`x\`.

- what:

  Whether to use only signal from on-targets (0 mismatches), off-targets
  (\>0 mm) or all. One of c("on-targets", "off-targets", "all).

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

# simulate 2 different experiments by picking 25 random offtargets
spcas9 <- c(offtargets[offtargets$mismatches == 0], sample(offtargets, 25))
lz3    <- c(offtargets[offtargets$mismatches == 0], sample(offtargets, 25))
plot_relative_activity(list(SpCas9=spcas9, LZ3=lz3), what="all")
```
