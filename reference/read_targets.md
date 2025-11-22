# read_targets

Read in a bed file with positions of breaks

## Usage

``` r
read_targets(x, genome, standard_chromosomes, strandless)
```

## Arguments

- x:

  GRanges object with coordinates of breaks or the path to a bed file.

- genome:

  a BSgenome object.

- standard_chromosomes:

  logical, constrain the analysis to standard chromosomes only (defaults
  to TRUE).

- strandless:

  collapse reads from multiple strands.

## Value

A GRanges object, unstranded and with a \`score\` column

## Examples

``` r
breaks <- read_targets(x=system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
                       genome="BSgenome.Hsapiens.UCSC.hg38",
                       standard_chromosomes=TRUE,
                       strandless=FALSE)
```
