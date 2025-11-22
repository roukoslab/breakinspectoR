# scission_profile_analysis

Analyses the scission profile, whether the Cas9 protein cut blunt or
staggered around the cut site. A significant p-value means that most
likely Cas9 performed a staggered cut.

## Usage

``` r
scission_profile_analysis(
  x,
  target,
  nontarget = GRanges(),
  region = 3,
  standard_chromosomes = TRUE,
  bsgenome = "BSgenome.Hsapiens.UCSC.hg38",
  qval_cutoff = 0.01,
  cores = getOption("mc.cores", 2L),
  verbose = TRUE
)
```

## Arguments

- x:

  A GRanges object containing the results of a breakinspectoR analysis.

- target:

  GRanges object with coordinates of breaks in the targeted library, or
  a character with the path to a bed file.

- nontarget:

  GRanges object with coordinates of breaks in the non-targeted library,
  or the path to a bed file. Default: empty GRanges() object, which is
  equivalent to not having a control library.

- region:

  region around the cutsite to look for alternative breaks.

- standard_chromosomes:

  logical, constrain the analysis to standard chromosomes only (defaults
  to TRUE).

- bsgenome:

  character, bsgenome to use (eg. BSgenome.Hsapiens.UCSC.hg38).

- qval_cutoff:

  numeric, indicating the q-value cutoff, defaults to 0.01.

- cores:

  Use multiple cores. Defaults to 2 cores, and the function won't
  benefit of adding \>2. Nevertheless, the speedup with 2 cores is
  considerable.

- verbose:

  logical, keep informing about every step

## Value

a GRanges object with the coordinates of the cutsite, the number of
breaks on and around the cutsite for the target and non-target
libraries, and a p- and q-values

## Details

The function looks at the number of breaks detected on and around the
cutsite using \*only\* the breaks detected on the PAM proximal site
(looking at reads on the same strand as the PAM). Then, a binomial test
similar to the one done in the \`breakinspectoR()\` function is
performed to test for a significant enrichment of breaks around the
cutsite in the target compared to the non-target libraries.

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

## this is needed only for the package to install
if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
##
  offtargets.scission_profile <- scission_profile_analysis(
    x        =offtargets,
    target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
    nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
    bsgenome ="BSgenome.Hsapiens.UCSC.hg38")

  mcols(offtargets) <- cbind(mcols(offtargets), mcols(offtargets.scission_profile))
}
#> Loading BSgenome.Hsapiens.UCSC.hg38...
#>  done
#> Importing breaks...
#>  done
#> Counting breaks at each offtarget region...
#>  done
#> Calculating probabilities...
#>  done
```
