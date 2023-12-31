% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotlib.R
\name{plot_scission_profile}
\alias{plot_scission_profile}
\title{plot_scission_profile
Plot ratios of blunt vs. staggered signal for each OT detected.}
\usage{
plot_scission_profile(x, type = c("frequency", "absolute"))
}
\arguments{
\item{x}{A GRanges object containing the results of a scission profile
analysis.}

\item{type}{show the ratio as absolute signal or percentage. One of
c("frequency", "absolute").}
}
\value{
A ggplot object.
}
\description{
plot_scission_profile
Plot ratios of blunt vs. staggered signal for each OT detected.
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

  offtargets.scission_profile <- scission_profile_analysis(
    x        =offtargets,
    target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
    nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
    bsgenome ="BSgenome.Hsapiens.UCSC.hg38")

  plot_scission_profile(offtargets.scission_profile)
}
}
