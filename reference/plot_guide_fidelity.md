# plot_guide_fidelity

Show, for each position of the guide, the fidelity of the offtarget
sequence as the percentage of matching bases.

## Usage

``` r
plot_guide_fidelity(x, guide, pam)
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

Nothing.

## Details

Currently this function has the limitation that does only \*fixed\*
matching, expecting the offtarget sequence and the guide to show the
same character. See ?\`Biostrings::lowlevel-matching\` for more
information.

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

plot_guide_fidelity(offtargets, guide="GACCCCCTCCACCCCGCCTC", pam="NGG")
```
