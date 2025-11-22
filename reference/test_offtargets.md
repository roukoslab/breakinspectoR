# test_offtargets

Test for enrichment in the offtarget regions vs. non-targeted library

## Usage

``` r
test_offtargets(target_breaks, target_total, nontarget_breaks, nontarget_total)
```

## Arguments

- target_breaks:

  numeric vector with the breaks in the target library.

- target_total:

  numeric with total number of breaks in the whole target library

- nontarget_breaks:

  numeric vector with the breaks in the non-target library.

- nontarget_total:

  numeric with total number of breaks in the whole non-target library

## Value

a data.frame with the enrichment, the p-value and the q-value.

## Details

test each region \* A = number of breaks within the region in the target
library \* B = total number of breaks in the whole target library \* C =
number of breaks within the region in the nontarget library \* D = total
number of breaks in the whole nontarget library the binomial model for
calculating the p-Value is: \* number of trials = A+1 \* number of
successes = B+1 \* estimated success probability in each trial = (C+1) /
(D+1) the enrichment is then calculated as: ((A+1) / (B+1)) / ((C+1) /
(D+1)) note that we add a pseudocount to avoid dividing by 0
