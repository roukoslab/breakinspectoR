#' breakinspectoR
#'
#' @description
#' Looks for enriched DNA regions targeted for cleavage by the CRISPR system
#'
#' @param target GRanges object with coordinates of breaks in the targeted
#' library, or a character with the path to a bed file.
#' @param nontarget GRanges object with coordinates of breaks in the
#' non-targeted library, or the path to a bed file. Default: empty GRanges()
#' object, which is equivalent to not having a control library.
#' @param guide character vector with the guide RNA sequence.
#' @param PAM character vector with the sequence of the protospacer adjacent
#' motif(s).
#' @param mismatches numeric, mismatches allowed between the guide and the
#' genomic sequence.
#' @param cutsiteFromPAM distance of the cutsite from the PAM.
#' @param min_breaks minimum number of breaks in the cutsite to be considered
#' for the analysis.
#' @param standard_chromosomes logical, constrain the analysis to standard
#' chromosomes only (defaults to TRUE).
#' @param bsgenome character, bsgenome to use (eg. BSgenome.Hsapiens.UCSC.hg38).
#' @param eFDR logical, estimate the empirical false discovery rate.
#' @param scale_nontarget logical, scale the nontarget library to the size of
#' the target library. Useful if target and nontarget libraries were sequenced
#' at different depths. Defaults to FALSE.
#' @param seed numeric, seed for random number generator. Disable with NA.
#' @param cores Use multiple cores. Defaults to 2 cores, and the function won't
#' benefit of adding >2. Nevertheless, the speedup with 2 cores is considerable.
#' @param verbose logical, keep informing about every step.
#'
#' @return a GRanges object with the coordinates of the enriched regions, with
#' number of breaks detected, the cut site, the enrichment over the
#' non-targeted library and the p-value of the statistical test.
#'
#' @import BSgenome
#' @import Biostrings
#' @import rtracklayer
#' @importFrom S4Vectors mcols decode
#' @importFrom stats rmultinom
#' @importFrom parallel mclapply
#' @export
#'
#' @examples
#' ## this is needed only for the package to install
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#' ##
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     eFDR     =FALSE,
#'     cutsiteFromPAM=3)
#' }
#'
breakinspectoR <- function(target, nontarget=GRanges(), guide, PAM=c("NAG", "NGG"),
                           mismatches=7, cutsiteFromPAM=3, min_breaks=2,
                           standard_chromosomes=TRUE,
                           bsgenome="BSgenome.Hsapiens.UCSC.hg38",
                           eFDR=TRUE,
                           scale_nontarget=FALSE,
                           seed=666,
                           cores=getOption("mc.cores", 2L),
                           verbose=TRUE) {

  # set seed for random number generation
  if(!is.na(seed)) {
    msg(verbose, "Setting seed for random number generation...", appendLF=FALSE)
    set.seed(seed)
    msg(verbose, " done", appendLF=TRUE)
  }

  # load the reference genome. Will be used almost everywhere...
  msg(verbose, "Loading ", bsgenome, "...", appendLF=FALSE)
  if(!is(bsgenome, "BSgenome")) {
    if (requireNamespace(bsgenome, quietly=TRUE)) {
      genome <- eval(parse(text=paste0(bsgenome, "::", bsgenome)))
    }
  } else {
    genome <- bsgenome
  }
  msg(verbose, " done", appendLF=TRUE)

  # read the bed files
  msg(verbose, "Importing breaks...", appendLF=FALSE)
  breaks <- mclapply(list(target=target, nontarget=nontarget), read_targets,
                     genome, standard_chromosomes, strandless=TRUE, mc.cores=cores)
  msg(verbose, " done", appendLF=TRUE)

  # scale the non-target library, if requested
  # how: draw 10 vectors of multinomially distributed random breaks in target
  #      with probability equal to the non-target breaks, and average them
  if(scale_nontarget) {
    msg(verbose, "Scaling non-target library...", appendLF=FALSE)
    breaks$nontarget$score <- rmultinom(1, sum(breaks$target$score), breaks$nontarget$score)
    breaks$nontarget <- breaks$nontarget[breaks$nontarget$score > 0]
    msg(verbose, " done", appendLF=TRUE)
  }

  # identify all loci in "target" that have at least `min_breaks` breaks
  offtarget_candidates <- breaks$target[breaks$target$score >= min_breaks]
  mcols(offtarget_candidates) <- NULL
  offtarget_candidates <- unique(offtarget_candidates)

  # add the sequence context around the candidate loci, and the PAM
  msg(verbose, "Obtaining sequence context...", appendLF=FALSE)
  s <- get_guide_and_pam(offtarget_candidates,
                         cutsiteFromPAM=cutsiteFromPAM,
                         glen=nchar(guide),
                         pamlen=max(nchar(PAM)),
                         genome=genome)
  msg(verbose, " done", appendLF=TRUE)

  # keep only offtarget candidates with PAM and few mismatches to the guide
  msg(verbose, "Filtering out candidates without PAM or too many mismatches...", appendLF=FALSE)
  guide_sense <- vcountPattern(guide, s$sense_guide, max.mismatch=mismatches, fixed=FALSE) > 0
  guide_as    <- vcountPattern(guide, s$as_guide   , max.mismatch=mismatches, fixed=FALSE) > 0
  pam_sense   <- rowSums(sapply(PAM, function(pam) vcountPattern(pam, s$sense_pam, fixed=FALSE))) > 0
  pam_as      <- rowSums(sapply(PAM, function(pam) vcountPattern(pam, s$as_pam   , fixed=FALSE))) > 0
  match_sense <- guide_sense & pam_sense
  match_as    <- guide_as    & pam_as

  offtarget_candidates$cutsite <- ranges(offtarget_candidates)
  offtarget_candidates$sense   <- ifelse(match_sense, "sense",       ifelse(match_as, "antisense", ""))
  offtarget_candidates$guide   <- ifelse(match_sense, s$sense_guide, ifelse(match_as, s$as_guide, ""))
  offtarget_candidates$pam     <- ifelse(match_sense, s$sense_pam,   ifelse(match_as, s$as_pam, ""))
  ranges(offtarget_candidates)[match_as]    <- s$as_loci[match_as]
  ranges(offtarget_candidates)[match_sense] <- s$sense_loci[match_sense]

  # add mismatches to the guide (take the "max" as it corresponds to the longest match)
  offtarget_candidates$mismatches <- neditAt(DNAString(guide),
                                             DNAStringSet(offtarget_candidates$guide),
                                             fixed=FALSE)

  # this is the final set of offtargets that have the guide+PAM with up to `mismatches`
  offtargets <- offtarget_candidates[(match_sense | match_as) & offtarget_candidates$mismatches <= mismatches]
  msg(verbose, " done", appendLF=TRUE)

  # at this point, return if no on/off targets were detected
  if(length(offtargets) == 0) {
    msg(verbose, "No targets detected", appendLF=FALSE)
    return(offtargets)
  }
  
  # get the number of breaks per position for those regions
  msg(verbose, "Counting breaks at each offtarget region...", appendLF=FALSE)
  x <- do.call(cbind, lapply(breaks, function(x) {

    breaks <- t(sapply(coverage(x, weight=x$score)[offtargets], decode))
    i <- offtargets$sense == "antisense"
    if(any(i)) {
      breaks[i, ] <- sapply(breaks[i, ], rev)
    }

    cutsite <- GRanges(seqnames(offtargets), IRanges(offtargets$cutsite))
    cutsite_breaks <- sapply(coverage(x, weight=x$score)[cutsite], decode)

    data.frame(breaks =rowSums(breaks),
               cutsite_breaks=cutsite_breaks)
  }))

  mcols(offtargets) <- cbind(mcols(offtargets), x)
  msg(verbose, " done", appendLF=TRUE)

  # test each region
  msg(verbose, "Calculating probabilities...", appendLF=FALSE)
  x <- test_offtargets(offtargets$target.breaks   , sum(breaks$target$score),
                       offtargets$nontarget.breaks, sum(breaks$nontarget$score))
  mcols(offtargets) <- cbind(mcols(offtargets), x)
  msg(verbose, " done", appendLF=TRUE)

  # calculate empirical FDR
  if(eFDR) {
    msg(verbose, "Calculating empirical FDR", appendLF=FALSE)
    null <- breaks$target
    null$score <- round(rowMeans(rmultinom(10, size=sum(null$score), prob=sample(null$score))))
    fdr <- breakinspectoR(null, breaks$nontarget, guide, PAM, mismatches,
                          cutsiteFromPAM, min_breaks, standard_chromosomes,
                          genome, eFDR=FALSE, verbose=FALSE)
    f <- Vectorize(function(x, mm) {
      sum(fdr$target.breaks[fdr$mismatches == mm] >= x) / sum(offtargets$target.breaks[offtargets$mismatches == mm] >= x)
    })
    offtargets$fdr <- f(offtargets$target.breaks, offtargets$mismatches)
    offtargets$fdr[offtargets$fdr > 1] <- 1
    msg(verbose, " done", appendLF=TRUE)
  }

  offtargets
}

#' test_offtargets
#'
#' @description
#' Test for enrichment in the offtarget regions vs. non-targeted library
#'
#' @param target_breaks numeric vector with the breaks in the target library.
#' @param target_total numeric with total number of breaks in the whole target
#' library
#' @param nontarget_breaks numeric vector with the breaks in the non-target
#' library.
#' @param nontarget_total numeric with total number of breaks in the whole
#' non-target library
#'
#' @details
#' test each region
#'   * A = number of breaks within the region in the target library
#'   * B = total number of breaks in the whole target library
#'   * C = number of breaks within the region in the nontarget library
#'   * D = total number of breaks in the whole nontarget library
#' the binomial model for calculating the p-Value is:
#'   * number of trials = A+1
#'   * number of successes = B+1
#'   * estimated success probability in each trial =  (C+1) / (D+1)
#' the enrichment is then calculated as: ((A+1) / (B+1)) / ((C+1) / (D+1))
#' note that we add a pseudocount to avoid dividing by 0
#'
#' @return a data.frame with the enrichment, the p-value and the q-value.
#'
#' @importFrom stats pbinom
#' @importFrom qvalue qvalue
#'
test_offtargets <- function(target_breaks, target_total,
                            nontarget_breaks, nontarget_total) {

  # get enrichment and pvalue
  x <- do.call(rbind, Map(target_breaks, nontarget_breaks,
                          f=function(target, nontarget) {
    p <- pbinom(q   =target + 1,
                size=target_total + 1,
                prob=((nontarget + 1) / (nontarget_total + 1)),
                lower.tail=FALSE)
    data.frame(enrichment=((target+1) / (target_total+1)) / ((nontarget+1) / (nontarget_total+1)),
               pval=p)
  }))

  # get qvalue
  x$qval <- tryCatch(qvalue(p=x$pval)$qvalues,
                     # handle error: could not estimate pi0.
                     # See https://support.bioconductor.org/p/105623/
                     error=function(e) tryCatch(qvalue(p=x$pval, lambda=0)$qvalues,
                                                error=function(e) NA))

  x
}

#' get_guide_and_pam
#'
#' @description
#' Get guide and PAM sequences around the cutsite
#'
#' @param x GRanges objects with loci to fetch the sequences
#' @param cutsiteFromPAM expected distance of the cutsite from the PAM
#' @param glen length of the guide sequence
#' @param pamlen length of the pam sequence
#' @param genome a bsgenome object (eg. BSgenome.Hsapiens.UCSC.hg38)
#'
#' @details
#' This function returns the sequence context of the targeted region. This is,
#' the guide+PAM. Since the coordinates of a break sits to the left of the locus
#' where the enzyme cuts regardless of the strand, the distance to the PAM will
#' be different for the two strands. For example:
#'
#' If a guide matches the *sense* strand:
#'   sgRNA=GACCCCCTCCACCCCGCCTC
#'
#'        <------ guide ------>pam
#'                        *|        <- `*`: break coord; `|`: actual cutsite
#' 5' (+) GACCCCCTCCACCCCGC|CTCNGG  <- sense (match)
#' 3' (+) CTGGGGGAGGTGGGGCG|GAGNCC  <- antisense (NO match)
#'                        0|123456  <- distance from cutsite coord
#'
#' If a guide matches the *antisense* strand:
#'             *|                   <- `*`: break  coord; `|`: actual cutsite
#' 5' (+) CCNGAG|GCGGGGTGGAGGGGGTC  <- sense (NO match)
#' 3' (-) GGNCTC|CGCCCCACCTCCCCCAG  <- antisense (match)
#'        543210|                   <- distance from cutsite coord
#'        pam<------ guide ------>
#' 
#'    which reversed would look like:
#'        <------ guide ------>pam
#'                         |*      <- `*`: break coord; `|`: actual cutsite
#' 5' (+) GACCCCCTCCACCCCGC|CTCNGG
#'                         |012345  <- distance from cutsite coord
#'
#'    and therefore the actual cutsite distance from PAM differs in the sense
#'    and antisense strands.
#'
#' @return A list of DNAStrings with the sense and antisense guide and pam, and
#' the loci.
#'
#' @import BSgenome
#' @import GenomicRanges
#' @import IRanges
get_guide_and_pam <- function(x, cutsiteFromPAM, glen, pamlen, genome) {

  # actual cutsite sits always to the left of the actual cutsite, regardless of
  # the strand. We always operate on the sense strand. See details
  cutsiteFromPAM_sense <- cutsiteFromPAM
  cutsiteFromPAM_as    <- cutsiteFromPAM - 1

  # getSeq, catching the error if coordinates were wrong
  getseq_robust <- function(genome, seqnames, start, end, strand, ...) {
      tryCatch({
        getSeq(genome, GRanges(paste0(seqnames, ":", start, "-", end, ":", strand), ...))
      }, error=function(e) NULL)
    }

  list(
    # get sense sequence. PAM is on the right of the cutsite
    sense_guide=getseq_robust(genome,
                              seqnames(x),
                              start =end(x) + cutsiteFromPAM_sense - glen + 1,
                              end   =end(x) + cutsiteFromPAM_sense,
                              strand="+"),
    sense_pam  =getseq_robust(genome,
                              seqnames(x),
                              start =end(x) + cutsiteFromPAM_sense + 1,
                              end   =end(x) + cutsiteFromPAM_sense + pamlen,
                              strand="+"),
    sense_loci =IRanges(start=end(x) + cutsiteFromPAM_sense - glen + 1, # coordinates of the guide
                        end  =end(x) + cutsiteFromPAM_sense),
    # get antisense sequence. PAM is on the left of the cutsite
    as_guide   =reverseComplement(   # reverse complement the sequence from the "+" strand
                  getseq_robust(genome,
                                seqnames(x),
                                start =end(x) - cutsiteFromPAM_as,
                                end   =end(x) - cutsiteFromPAM_as + glen - 1,
                                strand="+")),
    as_pam     =reverseComplement(   # reverse complement the sequence from the "+" strand
                  getseq_robust(genome,
                                seqnames(x),
                                start =end(x) - cutsiteFromPAM_as - pamlen,
                                end   =end(x) - cutsiteFromPAM_as - 1,
                                strand="+")),
    as_loci    =IRanges(start=end(x) - cutsiteFromPAM_as,
                        end  =end(x) - cutsiteFromPAM_as + glen - 1))
}

#' read_targets
#'
#' @description
#' Read in a bed file with positions of breaks
#'
#' @param x GRanges object with coordinates of breaks or the path to a bed file.
#' @param genome a BSgenome object.
#' @param standard_chromosomes logical, constrain the analysis to standard
#' chromosomes only (defaults to TRUE).
#' @param strandless collapse reads from multiple strands.
#'
#' @return A GRanges object, unstranded and with a `score` column
#'
#' @import rtracklayer
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @importFrom BiocGenerics unstrand
#' @importFrom S4Vectors mcols
#' @importFrom methods is
#' @export
#'
#' @examples
#' breaks <- read_targets(x=system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'                        genome="BSgenome.Hsapiens.UCSC.hg38",
#'                        standard_chromosomes=TRUE,
#'                        strandless=FALSE)
read_targets <- function(x, genome, standard_chromosomes, strandless) {
  # load genome if needed
  if(!is(genome, "BSgenome")) {
    if(requireNamespace(genome, quietly=TRUE)) {
      genome <- eval(parse(text=paste0(genome, "::", genome)))
    }
  }

  # read the bed file if it is not a GRanges object
  if(!is(x, "GRanges")) {
    x <- tryCatch({
      rtracklayer::import.bed(x, genome=seqinfo(genome))
    }, error=function(e) {
      if(e$message == "'seqnames' contains sequence names with no entries in 'seqinfo'") {
        message(e)
        message("\nTrying to proceed after dropping unknown chromosomes")
        x <- rtracklayer::import.bed(x)
        keepSeqlevels(x, intersect(seqlevels(genome), seqlevels(x)), pruning.mode="coarse")
      } else {
        stop(e)
      }
    })
  }

  # remove non-standard chromosomes if needed
  if(standard_chromosomes) {
    x <- keepStandardChromosomes(x, pruning.mode="coarse")
  }

  # add a column with the score if this doesn't exits. The score is the number
  # breaks at that locus (assume 1 if not present)
  if(! "score" %in% colnames(mcols(x))) {
    mcols(x)$score <- 1
  }

  # collapse reads from multiple strands
  if(strandless & !all(strand(x) == "*")) {
    xx <- GenomicRanges::reduce(BiocGenerics::unstrand(x), with.revmap=TRUE, min.gapwidth=0L)
    l  <- sapply(xx$revmap, length)
    if(any(l > 1)) {   # this code below is slow, do only if needed (l>1)
      revmap <- data.frame(x.row =unlist(xx$revmap),
                           xx.row=rep(1:length(xx), l))
      revmap$x.score <- x$score[revmap$x.row]
      xx$score <- as.numeric(tapply(revmap$x.score, revmap$xx.row, sum))
      x <- xx[, -1] # remove the revmap column
    }
  }

  x
}

#' reduceOT
#'
#' @description
#' Does a very opinionated filtering of the OT detected by breakinspector.
#'
#' @param x A GRanges object containing the results of a breakinspectoR
#' analysis.
#' @param fdr numerical, keep only hits below this empirical FDR threshold.
#' @param qval numerical, keep only hits below this q-value threshold.
#' @param mismatches numeric, mismatches allowed between the guide and the
#' genomic sequence.
#' @param standard_chromosomes logical, constrain the analysis to standard
#' chromosomes only (defaults to TRUE).
#' @param cores numerical, number of parallel cores to use.
#' @param verbose logical, keep informing about every step
#'
#' @return a GRanges object with the filtered OT.
#'
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @importFrom Biostrings neditAt
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom parallel mclapply
#' @export
#'
#' @examples
#' ## this is needed only for the package to install
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#' ##
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     eFDR     =FALSE,
#'     cutsiteFromPAM=3)
#'
#'   offtargets_filtered <- reduceOT(offtargets)
#' }
#'
reduceOT <- function(x, fdr=.01, qval=.01, mismatches=7, standard_chromosomes=TRUE,
                   cores=getOption("mc.cores", 2L),
                   verbose=TRUE) {

  # significant FDR and q-value
  msg(verbose, paste0("Keeping significant hits only (FDR<", fdr, " & qval<", qval, ")..."), appenLF=FALSE)
  x <- x[x$qval < qval]
  if("fdr" %in% colnames(mcols(x))) {   # only if empirical FDR was calculated
    x <- x[x$fdr  < fdr]
  }
  msg(verbose, " done", appendLF=TRUE)

  # remove NNN* protospacers
  msg(verbose, "Removing remove NNN* protospacers...", appendLF=FALSE)
  calcMM <- Vectorize(function(x, y) Biostrings::neditAt(x, y))
  x <- 
    tryCatch({
      x$mismatches <- calcMM(x$guide, guide)
       x[x$mismatches <= mismatches]
    }, error=function(e) x)
  msg(verbose, " done", appendLF=TRUE)

  # merge different OT in same protospacer
  msg(verbose, "Merging different OT in same protospacer...", appendLF=FALSE)
  x <- 
    tryCatch({
      coord <- paste(seqnames(x), start(x), end(x), x$sense) 
      x <- do.call(rbind, mclapply(split(x, coord), function(x) {
        i <- which.max(x$target.cutsite_breaks) # most often cutsite
        x[i]
      }, mc.cores=cores))
    }, error=function(e) x)
  msg(verbose, " done", appendLF=TRUE)

  # merge same OT in different protospacers: after merging OT in same protospacer,
  # will happen that for 2 different partially overlapping protospacers, the same
  # cutsite is identified as the most common cut site. Merge those protospacers into
  # 1 single entry
  msg(verbose, "Merging same OT in different protospacers...", appendLF=FALSE)
  x <-
    tryCatch({
      coord <- paste(seqnames(x), x$cutsite, x$sense) 
      x <- do.call(rbind, mclapply(split(x, coord), function(x) {
        i <- which.min(x$mismatches) # most fidel sequence to the sgRNA
        x[i]
      }, mc.cores=cores))
    }, error=function(e) x)
  msg(verbose, " done", appendLF=TRUE)

  # remove non-standard chromosomes if needed
  if(standard_chromosomes) {
    msg(verbose, "Keeping standard chromosomes only...", appendLF=FALSE)
    x <- keepStandardChromosomes(x, pruning.mode="coarse")
    msg(verbose, " done", appendLF=TRUE)
  }

  x
}

#' scission_profile_analysis
#'
#' @description
#' Analyses the scission profile, whether the Cas9 protein cut blunt or
#' staggered around the cut site. A significant p-value means that most likely
#' Cas9 performed a staggered cut.
#'
#' @param x A GRanges object containing the results of a breakinspectoR
#' analysis.
#' @param target GRanges object with coordinates of breaks in the targeted
#' library, or a character with the path to a bed file.
#' @param nontarget GRanges object with coordinates of breaks in the
#' non-targeted library, or the path to a bed file. Default: empty GRanges()
#' object, which is equivalent to not having a control library.
#' @param region region around the cutsite to look for alternative breaks.
#' @param standard_chromosomes logical, constrain the analysis to standard
#' chromosomes only (defaults to TRUE).
#' @param bsgenome character, bsgenome to use (eg. BSgenome.Hsapiens.UCSC.hg38).
#' @param qval_cutoff numeric, indicating the q-value cutoff, defaults to 0.01.
#' @param cores Use multiple cores. Defaults to 2 cores, and the function won't
#' benefit of adding >2. Nevertheless, the speedup with 2 cores is considerable.
#' @param verbose logical, keep informing about every step
#'
#' @return a GRanges object with the coordinates of the cutsite, the number of
#' breaks on and around the cutsite for the target and non-target libraries,
#' and a p- and q-values
#'
#' @details
#' The function looks at the number of breaks detected on and around the cutsite
#' using *only* the breaks detected on the PAM proximal site (looking at reads
#' on the same strand as the PAM).
#' Then, a binomial test similar to the one done in the `breakinspectoR()`
#' function is performed to test for a significant enrichment of breaks around
#' the cutsite in the target compared to the non-target libraries.
#' 
#' @import BSgenome
#' @import rtracklayer
#' @importFrom S4Vectors mcols decode
#' @importFrom stats rmultinom
#' @importFrom stats pbinom
#' @importFrom qvalue qvalue
#' @importFrom parallel mclapply
#' @export
#'
#' @examples
#' ## this is needed only for the package to install
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#' ##
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     eFDR     =FALSE,
#'     cutsiteFromPAM=3)
#'
#'   offtargets.scission_profile <- scission_profile_analysis(
#'     x        =offtargets,
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38")
#'
#'   mcols(offtargets) <- cbind(mcols(offtargets), mcols(offtargets.scission_profile))
#' }
#'
scission_profile_analysis <- function(x,
                                      target,
                                      nontarget=GRanges(),
                                      region=3,
                                      standard_chromosomes=TRUE,
                                      bsgenome="BSgenome.Hsapiens.UCSC.hg38",
                                      qval_cutoff=0.01,
                                      cores=getOption("mc.cores", 2L),
                                      verbose=TRUE) {

  # load the reference genome. Will be used almost everywhere...
  msg(verbose, "Loading ", bsgenome, "...", appendLF=FALSE)
  if(!is(bsgenome, "BSgenome")) {
    if (requireNamespace(bsgenome, quietly=TRUE)) {
      genome <- eval(parse(text=paste0(bsgenome, "::", bsgenome)))
    }
  } else {
    genome <- bsgenome
  }
  msg(verbose, " done", appendLF=TRUE)

  # read the bed files
  msg(verbose, "Importing breaks...", appendLF=FALSE)
  breaks <- mclapply(list(target=target, nontarget=nontarget), read_targets,
                     genome, standard_chromosomes, strandless=FALSE, mc.cores=cores)
  breaks$target    <- rep(breaks$target   , breaks$target$score)    # expand "score" into "reads" (individual entries in GRanges)
  breaks$nontarget <- rep(breaks$nontarget, breaks$nontarget$score)
  msg(verbose, " done", appendLF=TRUE)

  # count breaks on and around the target and non-target libraries
  msg(verbose, "Counting breaks at each offtarget region...", appendLF=FALSE)
  if(is(x$cutsite, "IRanges")) {
    x <- GRanges(seqnames(x), x$cutsite, ifelse(x$sense == "sense", "+", "-"))
  } else {
    x <- GRanges(seqnames(x), IRanges(x$cutsite, width=1), ifelse(x$sense == "sense", "+", "-"))
  }
  x.around           <- x
  ranges(x.around)   <- IRanges(start(x) - region, start(x) + region)
  x$on.target        <- countOverlaps(x       , breaks$target)
  x$around.target    <- countOverlaps(x.around, breaks$target)    - x$on.target
  x$on.nontarget     <- countOverlaps(x       , breaks$nontarget)
  x$around.nontarget <- countOverlaps(x.around, breaks$nontarget) - x$on.nontarget
  msg(verbose, " done", appendLF=TRUE)

  # test the enrichment of breaks around the cutsite in the target vs. non-target libs
  msg(verbose, "Calculating probabilities...", appendLF=FALSE)
  test <- function(around.target, around.nontarget) {
    pbinom(q   =around.target + 1,
           size=length(breaks$target) + 1,
           prob=((around.nontarget + 1) / (length(breaks$nontarget) + 1)),
           lower.tail=FALSE)
  }
  vtest <- Vectorize(test)
  x$pval <- vtest(x$around.target, x$around.nontarget)

  # get qvalue
  x$qval <- tryCatch(qvalue(p=x$pval)$qvalues,
                     # handle error: could not estimate pi0.
                     # See https://support.bioconductor.org/p/105623/
                     error=function(e) tryCatch(qvalue(p=x$pval, lambda=0)$qvalues,
                                                error=function(e) NA))
  msg(verbose, " done", appendLF=TRUE)

  # add column informing whether the break was significantly staggered (or not, ==blunt)
  x$scission_profile <- ifelse(x$around.target > x$on.target & x$qval < qval_cutoff, "staggered", "blunt")

  # count breaks at each individual position of the region, and add them as columns
  target_counts_per_position <- do.call(cbind, lapply(-region:region, function(pos) {
    x <- GRanges(seqnames(x), IRanges(start=start(x) + pos, end=start(x) + pos), strand=strand(x))
    countOverlaps(x, breaks$target)
  }))
  target_counts_per_position.rev <- t(apply(target_counts_per_position, 1, rev))  # reversed matrix, for the minus strand
  i <- decode(strand(x)) == "-"
  target_counts_per_position[i, ] <- target_counts_per_position.rev[i, ]
  colnames(target_counts_per_position) <- paste0("cutsite+", as.character(-region:region)) |>
                                          sub("\\+-", "-", x=_) |>
                                          sub("\\+0", "" , x=_)
  x$target_counts_per_position <- target_counts_per_position

  x
}

#' msg
#'
#' @description
#' Write messages to the console
#'
#' @param verbose logical, keep informing about every step
#' @param ... other parameters sent to `message()`
#'
#' @return nothing
msg <- function(verbose, ...) {
  if(verbose) {
    message(...)
  }
}
