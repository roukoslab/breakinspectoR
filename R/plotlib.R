#' plot_guide_fidelity
#'
#' @description
#' Show, for each position of the guide, the fidelity of the offtarget sequence
#' as the percentage of matching bases.
#'
#' @details
#' Currently this function has the limitation that does only *fixed* matching,
#' expecting the offtarget sequence and the guide to show the same character.
#' See ?`Biostrings::lowlevel-matching` for more information.
#'
#' @param x A GRanges object containing the results of a breakinspectoR
#' analysis.
#' @param guide character with the guide sequence. Only used for drawing, the
#' number of mismatches is already present in `x`.
#' @param pam character with PAM sequence. Only used for drawing, the number of
#' mismatches is already present in `x`.
#'
#' @return Nothing.
#'
#' @import pheatmap
#' @importFrom grid grid.newpage grid.text pushViewport viewport
#' @importFrom grDevices palette
#' @importFrom viridis viridis
#' @export
#'
#' @examples
#' \dontrun{
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     cutsiteFromPAM=3
#'   )
#'
#'   plot_guide_fidelity(offtargets, guide="GACCCCCTCCACCCCGCCTC", pam="NGG")
#' }
plot_guide_fidelity <- function(x, guide, pam) {
  guide_pam <- unlist(strsplit(paste0(guide, pam), ""))
  offtarget_seqs <- cbind(as.data.frame(do.call(rbind, strsplit(x$guide, ""))),
                          as.data.frame(do.call(rbind, strsplit(x$pam  , ""))))
  colnames(offtarget_seqs) <- paste("pos", seq_len(ncol(offtarget_seqs)))

  matches <- do.call(rbind, lapply(split(offtarget_seqs, x$mismatches), function(x) {
    rowSums(apply(x, 1, function(x) x == guide_pam)) / nrow(x)  # % of matches per position of the guide
  }))

  # do the plot
  setHook("grid.newpage", action="prepend", function() {
    pushViewport(viewport(x=1, y=1, width=0.98, height=0.98, name="vp", just=c("right","top")))
  })

  annotation_col <- data.frame(PAM      =c(rep("guide", nchar(guide)), rep("PAM", nchar(pam))),
                              nucleotide=guide_pam, row.names=colnames(offtarget_seqs))

  pheatmap(as.matrix(matches),
            color=rev(viridis(100)),
            main="off-target fidelity by number of mismatches",
            labels_col=guide_pam,
            angle_col=0,
            annotation_col=annotation_col,
            annotation_colors=list(PAM=c(guide=NA, PAM=palette()[1])),
            annotation_legend=TRUE,
            annotation_names_col=FALSE,
            cluster_cols=FALSE,
            cluster_rows=FALSE,
            breaks=seq(0, 1, length.out=100))

  setHook("grid.newpage", action="replace", NULL)
  grid.text("guide sequence position", y=-0.01)
  grid.text("number of mismatches", x=-0.01, rot=90)
}

#' plot_sequence_composition
#' Sequence composition of off-target sites.
#'
#' @param x A GRanges object containing the results of a breakinspectoR
#' analysis.
#' @param guide character with the guide sequence. Only used for drawing, the
#' number of mismatches is already present in `x`.
#' @param pam character with PAM sequence. Only used for drawing, the number of
#' mismatches is already present in `x`.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom scales percent
#' @importFrom reshape2 melt
#' @importFrom grDevices palette
#' @export
#'
#' @examples
#' \dontrun{
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     cutsiteFromPAM=3
#'   )
#'
#'   plot_sequence_composition(offtargets, guide="GACCCCCTCCACCCCGCCTC", pam="NGG")
#' }
plot_sequence_composition <- function(x, guide, pam) {

  guide_pam <- unlist(strsplit(paste0(guide, pam), ""))
  offtarget_seqs <- cbind(as.data.frame(do.call(rbind, strsplit(x$guide, ""))),
                          as.data.frame(do.call(rbind, strsplit(x$pam  , ""))))
  colnames(offtarget_seqs) <- paste("pos", seq_len(ncol(offtarget_seqs)))

  x <- reshape2::melt(lapply(offtarget_seqs, table))  # offtarget_seqs is a dataframe, table by column (each column is a position in the guide)
  x$value <- x$value / nrow(offtarget_seqs)
  x$L1 <- factor(x$L1, levels=unique(x$L1)[order(as.numeric(sub("pos ", "", unique(x$L1))))])

  p <- ggplot(x, aes(x=as.numeric(L1), y=value, color=Var1)) +
         geom_point() +
         geom_line() +
         annotate("segment", x=nchar(guide) + 0.5,
                          xend=nchar(guide) + nchar(pam) + 0.5,
                          y=0, yend=0, colour=palette()[1], size=3) +
         scale_color_manual("", values=palette()) +
         scale_y_continuous(labels=scales::percent) +
         scale_x_continuous(breaks=1:length(guide_pam), labels=guide_pam) +
         labs(title="Sequence composition of off-target sites",
                       x="guide sequence", y="") +
         theme_bw() +
         theme(legend.position="bottom")

  # add percentage off mismatches as a bar
  mm <- rowSums(apply(offtarget_seqs, 1, function(x) x != guide_pam)) / nrow(offtarget_seqs)  # % of matches per position of the sgRNA
  mm <- reshape2::melt(mm)
  mm$x <- factor(rownames(mm), levels=rownames(mm))

  p <- p + geom_bar(aes(x=as.numeric(x), y=value, fill="mismatch rate"),
                             data=mm, inherit.aes=FALSE, stat="identity",
                             alpha=1/5, width=1/2) +
           scale_fill_manual("", values="red")

  return(p)
}

#' plot_offtargets_by_pam
#' Plot the count of offtargets by the different PAMs
#'
#' @param x A GRanges object containing the results of a breakinspectoR
#' analysis.
#' @param fraction logical indicating whether the plot should show absolute
#' absolute numbers or the fraction.
#'
#' @import ggplot2
#' @importFrom scales percent
#' @importFrom reshape2 melt
#' @export
#'
#' @return A ggplot object.
#' @examples
#' \dontrun{
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     cutsiteFromPAM=3
#'   )
#'
#'   plot_offtargets_by_pam(offtargets, fraction=TRUE)
#' }
plot_offtargets_by_pam <- function(x, fraction=TRUE) {

  x <- reshape2::melt(table(x$pam))

  if(!fraction) {
    p <- ggplot(x, aes_string(x="Var1", y="value")) +
           geom_bar(stat="identity") +
           labs(title="offtargets detected by PAM", y="", x="") +
           theme_bw() +
           theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  } else {
    x$perc <- x$value / sum(x$value)
    p <- ggplot(x, aes_string(x="Var1", y="perc")) +
           geom_bar(stat="identity") +
           labs(title="offtargets detected by PAM", y="", x="") +
           scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
           theme_bw() +
           theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  }

  return(p)
}


#' plot_mismatch_freq
#' Plot the frequency of mismatches per offtarget.
#'
#' @param x A GRanges object containing the results of a breakinspectoR
#' analysis.
#' @param bins number of bins tu use in the histogram. Defaults to the maximum
#' number of mismatches + 1.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     cutsiteFromPAM=3
#'   )
#'
#'   plot_mismatch_freq(offtargets)
#' }
plot_mismatch_freq <- function(x, bins=1+diff(range(x$mismatches))) {
  x <- as.data.frame(x)

  p <- ggplot(x, aes_string(x="mismatches")) +
         geom_histogram(bins=bins, color="white") +
         labs(title="offtargets detected by number of mismatches",
              y="", x="mismatches") +
         theme_bw()

  return(p)
}

#' plot_position_cutsite
#' Plot the frequency of the position of the cutsite relative to the guide.
#'
#' @param x A GRanges object containing the results of a breakinspectoR
#' analysis.
#' @param guide character with the guide sequence. Only used for drawing.
#' @param pam character with PAM sequence. Only used for drawing.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom scales percent
#' @importFrom reshape2 melt
#' @importFrom grDevices palette
#' @export
#'
#' @examples
#' \dontrun{
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     cutsiteFromPAM=3
#'   )
#'
#'   plot_position_cutsite(offtargets, guide="GACCCCCTCCACCCCGCCTC", pam="NGG")
#' }
plot_position_cutsite <- function(x, guide, pam) {

  guide_pam <- unlist(strsplit(paste0(guide, pam), ""))
  x$target.cutsite <- ifelse(x$sense == "sense", as.numeric(start(x$cutsite)) - as.numeric(start(x) - 1),
                                                 as.numeric(end(x)) - as.numeric(start(x$cutsite)))

  x <- as.data.frame(table(x$target.cutsite))
  x$Freq <- x$Freq / sum(x$Freq)

  p <- ggplot(x, aes(x=as.numeric(as.character(Var1)), y=Freq)) +
         geom_bar(stat="identity") +
         scale_y_continuous(labels=scales::percent, limits=c(0, 1)) +
         scale_x_continuous(breaks=seq_len(length(guide_pam)), labels=guide_pam, limits=c(1, length(guide_pam)) + 0.5) +
         labs(title="Position of the cut site",
              x="guide sequence", y="cut frequency") +
         annotate("segment", x=nchar(guide) + 0.5,
                          xend=nchar(guide) + nchar(pam) + 0.5,
                          y=0, yend=0, colour=palette()[1], size=3) +
         theme_bw()

  return(p)
}

#' plot_genomic_position
#' Plot the genomic position of the breaks, identifying detected
#' on-/off-targets.
#'
#' @param x A GRanges object containing the results of a breakinspectoR
#' analysis. Possibly filtered (eg. qval < .01).
#' @param target GRanges object with coordinates of breaks in the targeted
#' library, or a character with the path to a bed file.
#' @param bsgenome character, bsgenome to use (eg. BSgenome.Hsapiens.UCSC.hg38).
#' @param standard_chromosomes logical, constrain the analysis to standard
#' chromosomes only. Default: TRUE.
#' @param min_breaks minimum number of breaks in the cutsite to be considered
#' for the analysis. Default: 2.
#' @param log_signal display log10 transformed signal. Default: FALSE.
#' @param verbose logical, keep informing about every step.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grDevices palette
#' @export
#'
#' @examples
#' \dontrun{
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     cutsiteFromPAM=3
#'   )
#'
#'   plot_genomic_position(
#'     offtargets,
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38"
#'   )
#' }
plot_genomic_position <- function(x, target, bsgenome,
                                  standard_chromosomes=TRUE,
                                  min_breaks=2, log_signal=FALSE, 
                                  verbose=TRUE) {

  # load the reference genome. Will be used almost everywhere...
  msg(verbose, "Loading ", bsgenome, "...", appendLF=FALSE)
  if (requireNamespace(bsgenome, quietly=TRUE)) {
    genome <- eval(parse(text=paste0(bsgenome, "::", bsgenome)))
  }
  msg(verbose, " done", appendLF=TRUE)

  # read the bed files
  msg(verbose, "Importing breaks...", appendLF=FALSE)
  breaks <- read_targets(target, genome, standard_chromosomes, strandless=TRUE)
  breaks <- breaks[breaks$score >= min_breaks]
  msg(verbose, " done", appendLF=TRUE)

  # scale
  chr_len <- as.numeric(seqlengths(seqinfo(breaks)))
  chr_genomic_starts <- cumsum(c(1, chr_len[-length(chr_len)]))
  names(chr_genomic_starts) <- names(seqlengths(seqinfo(breaks)))
  breaks$genomic_pos <- end(breaks) + chr_genomic_starts[match(as.character(seqnames(breaks)), names(chr_genomic_starts))]

  # annotate
  ranges(x) <- x$cutsite
  breaks$type <- ifelse(breaks %over% x[x$mismatches == 0], "ontarget",
                 ifelse(breaks %over% x, "offtarget", "break"))
  
  # do the plot
  df <- as.data.frame(breaks)

  p <- ggplot(mapping=aes(x=genomic_pos, y=score)) +
         geom_point(data=df, mapping=aes(color="breaks")) +
         scale_color_manual(values=c(breaks="#80808080", offtarget="#FF7F0E", ontarget="#D62728")) +
         labs(x="", y="number of breaks") +
         theme_bw() +
         theme(legend.position="bottom") +
         theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

  if(any(df$type == "offtarget")) {
    p <- p + geom_point(data=subset(df, type=="offtarget"), mapping=aes(color="offtarget"), shape=21, size=1.5)
  }
  if(any(df$type == "ontarget")) {
    p <- p + geom_point(data=subset(df, type=="ontarget"), mapping=aes(color="ontarget"), size=4)
  }

  # calculate the axis breaks and labels (as in select `selectMethod("autoplot", "GRanges")`)
  chr_genomic_mids <- chr_genomic_starts + chr_len / 2
  p <- p + scale_x_continuous(breaks=chr_genomic_mids, labels=names(seqlengths(seqinfo(breaks))))

  # log transform y-axis signal
  if(log_signal) {
    p <- p + scale_y_log10()
  }

  return(p)
}

#' manhattan_plot
#' Manhattan plot showing on/off-targets detected organized by chromosomal
#' position with bar height representing read count (number of breaks).
#'
#' @param x A GRanges object containing the results of a breakinspectoR
#' analysis. Possibly filtered (eg. qval < .01).
#' @param bsgenome character, bsgenome to use (eg. BSgenome.Hsapiens.UCSC.hg38).
#' @param standard_chromosomes logical, constrain the analysis to standard.
#' chromosomes only. Default: TRUE.
#' @param cutsite_breaks_only display breaks of cutsite only, instead of using
#' breaks accumulated in the complete protospacer. Default: TRUE
#' @param subtract_nontarget subtract breaks of non-target library.
#' Default: FALSE
#' @param log_signal display log2 transformed signal. Default: FALSE.
#' @param verbose logical, keep informing about every step.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grDevices palette
#' @export
#'
#' @examples
#' \dontrun{
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     cutsiteFromPAM=3
#'   )
#'
#'   manhattan_plot(
#'     offtargets,
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38"
#'   )
#' }
manhattan_plot <- function(x,
                           bsgenome,
                           standard_chromosomes=TRUE,
                           cutsite_breaks_only=TRUE,
                           subtract_nontarget=FALSE,
                           log_signal=FALSE, 
                           verbose=TRUE) {

  # load the reference genome. Will be used almost everywhere...
  msg(verbose, "Loading ", bsgenome, "...", appendLF=FALSE)
  if (requireNamespace(bsgenome, quietly=TRUE)) {
    genome <- eval(parse(text=paste0(bsgenome, "::", bsgenome)))
  }
  msg(verbose, " done", appendLF=TRUE)

  # drop non-standard chromosomes
  if(standard_chromosomes) {
    genome <- keepStandardChromosomes(seqinfo(genome))
  } else {
    genome <- seqinfo(genome)
  }


  # scale
  chr_len <- as.numeric(seqlengths(genome))
  chr_genomic_starts <- cumsum(c(1, chr_len[-length(chr_len)]))
  names(chr_genomic_starts) <- names(seqlengths(genome))
  i <- match(as.character(seqnames(x)), names(chr_genomic_starts))
  x$continuous_chr_cutsite <- start(x$cutsite) + chr_genomic_starts[i]

  # calculate the signal to plot
  if(cutsite_breaks_only) {
    x$score <- x$target.cutsite_breaks
    if(subtract_nontarget) {
      x$score <- x$score - x$nontarget.cutsite_breaks
    }
  } else {
    x$score <- x$target.breaks
    if(subtract_nontarget) {
      x$score <- x$score - x$nontarget.breaks
    }
  }


  # do the plot
  df <- as.data.frame(x)
  df$chr <- as.character(seqnames(x))

  p <- ggplot(df, aes(x=continuous_chr_cutsite, y=score, color=chr, fill=chr)) +
         geom_bar(stat="identity") +
         labs(x="", y="number of breaks") +
         theme_minimal() +
         theme(legend.position="none") +
         theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

  # add an arrow to the on-target
  if(any(df$mismatches == 0)) {
    ot    <- df$mismatches == 0
    x     <- df$continuous_chr_cutsite[ot]
    yend  <- df$score[ot]
    if(log_signal) {
      y <- 10^(log10(yend) + diff(log10(range(df$score))) / 10)
    } else {
      y <- yend + diff(range(df$score)) / 10
    }
    p  <- p +
            geom_point(data=df[ot, ], color="blue", size=3) +
            annotate("segment", x=x, xend=x, y=y, yend=yend, arrow=arrow(), color="blue")
  }

  # calculate the axis breaks and labels (as in select `selectMethod("autoplot", "GRanges")`)
  chr_genomic_mids <- chr_genomic_starts + chr_len / 2
  p <- p + scale_x_continuous(limits=c(1, max(chr_genomic_starts + chr_len)),
                              breaks=chr_genomic_mids,
                              labels=names(seqlengths(genome)))

  # log transform y-axis signal
  if(log_signal) {
    p <- p + scale_y_log10()
  }

  return(p)
}

#' plot_scission_profile
#' Plot ratios of blunt vs. staggered signal for each OT detected.
#' 
#'
#' @param x A GRanges object containing the results of a scission profile
#' analysis.
#' @param type show the ratio as absolute signal or percentage. One of
#' c("frequency", "absolute").
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom grDevices palette
#' @export
#'
#' @examples
#' \dontrun{
#'   offtargets <- breakinspectoR(
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     guide    ="GACCCCCTCCACCCCGCCTC",
#'     PAM      =c("NGG", "NAG"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38",
#'     cutsiteFromPAM=3
#'   )
#'
#'   offtargets.scission_profile <- scission_profile_analysis(
#'     x        =offtargets,
#'     target   =system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR"),
#'     nontarget=system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR"),
#'     bsgenome ="BSgenome.Hsapiens.UCSC.hg38")
#'
#'   plot_scission_profile(offtargets.scission_profile)
#' }
plot_scission_profile <- function(x, type=c("frequency", "absolute")) {
  type <- match.arg(type)

  if(is(x, "GRanges")) {
    x <- as.data.frame(x)
  }

  x$name <- paste(x$seqnames, x$start, x$end)
  x      <- reshape2::melt(x[, c("name", "on.target", "around.target", "scission_profile")],
                            id.vars=c("name", "scission_profile"))
  signal <- tapply(x$value, x$name, sum)
  ratio  <- tapply(x$value, x$name, function(x) (x[1])/(x[2]))
  x$name <- factor(x$name, levels=names(ratio)[order(ratio, signal, na.last=FALSE)])

  p <- ggplot(x, aes(x=name, y=value, fill=variable)) +
         geom_bar(stat="identity", position=ifelse(type == "absolute", "stack", "fill"), color=NA, size=0) +
         geom_rug(aes(y=0, color=scission_profile), sides="b") +
         labs(title="DSB unique UMIs", x="", y="") +
         scale_fill_manual("signal", values=palette()) +
         scale_color_manual("scission profile", values=c(NA, "red")) +
         theme_bw() +
         theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
               panel.grid.major=element_blank(), panel.grid.minor=element_blank())

  p
}
