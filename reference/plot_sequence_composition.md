# plot_sequence_composition Sequence composition of off-target sites.

plot_sequence_composition Sequence composition of off-target sites.

## Usage

``` r
plot_sequence_composition(x, guide, pam)
```

## Arguments

- x:

  A GRanges object containing the results of a breakinspectoR analysis.

- guide:

  character with the guide sequence. Only used for drawing, the number
  of mismatches is already present in \`x\`.

- pam:

  character with PAM sequence. Only used for drawing, the number of
  mismatches is already present in \`x\`.

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

plot_sequence_composition(offtargets, guide="GACCCCCTCCACCCCGCCTC", pam="NGG")
```
