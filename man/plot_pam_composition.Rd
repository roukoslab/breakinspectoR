% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotlib.R
\name{plot_pam_composition}
\alias{plot_pam_composition}
\title{plot_pam_composition
Plots the proportion nucleotide usage per position of the PAM.}
\usage{
plot_pam_composition(x)
}
\arguments{
\item{x}{A GRanges object containing the results of a breakinspectoR
analysis.}
}
\value{
A ggplot object.
}
\description{
plot_pam_composition
Plots the proportion nucleotide usage per position of the PAM.
}
\examples{
\dontrun{
  offtargets <- breakinspectoR(
    target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
    nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
    guide    ="GACCCCCTCCACCCCGCCTC",
    PAM      =c("NGG", "NAG"),
    bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
    cutsiteFromPAM=3
  )
}
data(breakinspectoR_examples)

plot_pam_composition(offtargets)
}
