# breakinspectoR
A companion R package to the BreakTag protocol for the identification of CRISPR offtargets.
breakinspectoR is an R package which performs a guided search toward putative on-/off-targets.

## Overall description

![breakinspectoR workflow](assets/breakinspectoR_wf.png)

### Initial preprocessing

Initial preprocessing is typically done in a linux cluster using the [Breaktag pipeline](https://github.com/roukoslab/breaktag). It includes the following steps:

1. scanning for reads (single- or paired-end) containing the expected 8-nt UMI followed by the 8-nt sample barcode in the 5' end of read 1.
2. alignment of reads to reference genome with BWA, with a seed length of 19 and default scoring/penalty values for mismatches, gaps and read clipping.
3. reads mapped with a minimum quality score Q (defaults to Q=60) are retained.
4. close spatial consecutive reads within a window of 30 nucleotides and UMI differing with up to 2 mismatches are considered PCR duplicates and only one is kept.
    
The resulting reads are aggregated per position and reported as a BED file.
    
### breakinspectoR analysis

The analysis path consists of the following steps:

1. from the previously generated BED file, identify stacks of read ends as the candidate loci of being CRISPR-edited.
2. obtain the sequence context of the candidate loci, and keep only those which are at most `N` nucleotides upstream from a PAM and contain up to M mismatches with regard to the gDNA guide sequence. Defaults are the canonical N=3, M=7, PAM="NGG".
3. count number of reads (== signal or DSBreaks) in the targeted and the non-targeted control library.
4. test if the enrichment of reads we see in the targeted library is significant compared to the non-targeted control library. Here breakinspectoR will perform a binomial test with the following criteria:
    1. preconditions:
        1. A = number of breaks within the region in the target library
        2. B = total number of breaks in the whole target library
        3. C = number of breaks within the region in the nontarget library
        4. D = total number of breaks in the whole nontarget library
    2. the binomial model for calculating the p-Value is:
        1. number of trials = A+1
        2. number of successes = B+1
        3. estimated success probability in each trial =  (C+1) / (D+1)
    3. the enrichment is then calculated as: ((A+1) / (B+1)) / ((C+1) / (D+1))
    4. note that we add a pseudocount to avoid dividing by 0
    5. q-values and local False Discovery Rate values are estimated for FDR control.
5.additionally, breakinspectoR implements a complimentary [and elaborated] method to estimate the false discovery of targets. To summarize, breakinspectoR reshuffles the signal in the target library using several multinomially distributed random number vectors sampled with equal probabilities to the signal in the originally detected offtargets. Then, breakinspectoR analysis is done in the reshuffled target library vs. the non-target library, and an FDR is estimated comparing the signal of each offtarget called in the original target library to the targets detected in the reshuffled target library (where no targets were expected to be called).
6. breakinspectoR includes several handy visualizations to further analyze and summarize the on-/off-targets detected. Some of these functions include the analysis of fidelity of the gDNA, sequence composition of target regions, frequency of mismatches per position of the protospacer, or the genomic distribution of the targeted regions. 
                  
## Example usage

This is a simple example using the demo data for human chr6 included with the package.
The experiment identifies offtargets generated by the VEGFA site 2 sgRNA with CRISPR/Cas9.

Call the breakinspectoR analysis to find offtargets enriched in the targeted library
compared to the non-targeted. We'll stick to the default 7 mismatches allowed to the
guide, with the expected cut site 3 bp away from the PAM (Cas9).

```R
target_file     <- system.file("extdata/vegfa.chr6.bed.gz", package="breakinspectoR")
non_target_file <- system.file("extdata/nontarget.chr6.bed.gz", package="breakinspectoR")
guide           <- "GACCCCCTCCACCCCGCCTC"
PAM             <- c(canonical="NGG", "NAG")
bsgenome        <- "BSgenome.Hsapiens.UCSC.hg38"

offtargets <- breakinspectoR(
    target        =target_file,
    nontarget     =non_target_file,
    guide         =guide,
    PAM           =PAM,
    bsgenome      =bsgenome,
    cutsiteFromPAM=3,
    verbose       =FALSE
)
```

The analysis will take few seconds to run. Afterwards we have a comprehensive table with
few hundred enriched offtarget loci, which we can summarize using the accompanying
plotting functions:

```R
plot_position_cutsite(offtargets, guide=guide, pam=PAM["canonical"])

plot_mismatch_freq(offtargets)

plot_offtargets_by_pam(offtargets)

plot_sequence_composition(offtargets, guide=guide, pam=PAM["canonical"])

plot_guide_fidelity(offtargets, guide=guide, pam=PAM["canonical"])

manhattan_plot(offtargets, bsgenome=bsgenome)
```

## Input files

Input BED files describing coordinates and number of DSB are expected for the "target" and "non-target" libraries.
It's possible to run breakinspectoR without the non-target library, nevertheless it is advised to include such experiment to calculate a p-value and control the false discovery rate.

```
chr6  148074  148075  .  2  +
chr6  148093  148094  .  1  -
chr6  148240  148241  .  1  -
chr6  148503  148504  .  1  +
chr6  148636  148637  .  1  -
chr6  148697  148698  .  1  -
chr6  149009  149010  .  1  -
chr6  149363  149364  .  1  +
chr6  150252  150253  .  2  +
chr6  150263  150264  .  1  +
```

These files are typically created with the [Breaktag pipeline](https://github.com/roukoslab/breaktag), although any conforming BED file is accepted in breakinspectoR.
