# reduceOT

Does a very opinionated filtering of the OT detected by breakinspector.

## Usage

``` r
reduceOT(
  x,
  fdr = 0.01,
  qval = 0.01,
  mismatches = 7,
  standard_chromosomes = TRUE,
  cores = getOption("mc.cores", 2L),
  verbose = TRUE
)
```

## Arguments

- x:

  A GRanges object containing the results of a breakinspectoR analysis.

- fdr:

  numerical, keep only hits below this empirical FDR threshold.

- qval:

  numerical, keep only hits below this q-value threshold.

- mismatches:

  numeric, mismatches allowed between the guide and the genomic
  sequence.

- standard_chromosomes:

  logical, constrain the analysis to standard chromosomes only (defaults
  to TRUE).

- cores:

  numerical, number of parallel cores to use.

- verbose:

  logical, keep informing about every step

## Value

a GRanges object with the filtered OT.

## Examples

``` r
if (FALSE) { # \dontrun{
  offtargets <- breakinspectoR(
    target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
    nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
    guide    ="GACCCCCTCCACCCCGCCTC",
    PAM      =c("NGG", "NAG"),
    bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
    eFDR     =FALSE,
    cutsiteFromPAM=3)
} # }
data(breakinspectoR_examples)

offtargets_filtered <- reduceOT(offtargets)
#> Keeping significant hits only (FDR<0.01 & qval<0.01)...FALSE
#>  done
#> Removing remove NNN* protospacers...
#>  done
#> Merging different OT in same protospacer...
#>  done
#> Merging same OT in different protospacers...
#>  done
#> Keeping standard chromosomes only...
#>  done
```
