# get_guide_and_pam

Get guide and PAM sequences around the cutsite

## Usage

``` r
get_guide_and_pam(x, cutsiteFromPAM, glen, pamlen, genome)
```

## Arguments

- x:

  GRanges objects with loci to fetch the sequences

- cutsiteFromPAM:

  expected distance of the cutsite from the PAM

- glen:

  length of the guide sequence

- pamlen:

  length of the pam sequence

- genome:

  a bsgenome object (eg. BSgenome.Hsapiens.UCSC.hg38)

## Value

A list of DNAStrings with the sense and antisense guide and pam, and the
loci.

## Details

This function returns the sequence context of the targeted region. This
is, the guide+PAM. Since the coordinates of a break sits to the left of
the locus where the enzyme cuts regardless of the strand, the distance
to the PAM will be different for the two strands. For example:

If a guide matches the \*sense\* strand: sgRNA=GACCCCCTCCACCCCGCCTC

\<—— guide ——\>pam \*\| \<- \`\*\`: break coord; \`\|\`: actual cutsite
5' (+) GACCCCCTCCACCCCGC\|CTCNGG \<- sense (match) 3' (+)
CTGGGGGAGGTGGGGCG\|GAGNCC \<- antisense (NO match) 0\|123456 \<-
distance from cutsite coord

If a guide matches the \*antisense\* strand: \*\| \<- \`\*\`: break
coord; \`\|\`: actual cutsite 5' (+) CCNGAG\|GCGGGGTGGAGGGGGTC \<- sense
(NO match) 3' (-) GGNCTC\|CGCCCCACCTCCCCCAG \<- antisense (match)
543210\| \<- distance from cutsite coord pam\<—— guide ——\>

which reversed would look like: \<—— guide ——\>pam \|\* \<- \`\*\`:
break coord; \`\|\`: actual cutsite 5' (+) GACCCCCTCCACCCCGC\|CTCNGG
\|012345 \<- distance from cutsite coord

and therefore the actual cutsite distance from PAM differs in the sense
and antisense strands.
