# 02120 qPCR primers and covid-19
#qPCR
# qPCR primer specificity across SARS-CoV-2 variants

02-120 Programming for Scientists, Fall 2026.

## Scientific question

Do qPCR primer sets designed early in the pandemic still match the genomes of
Alpha, Delta, and Omicron? Which targets keep matching, which acquire
mismatches, and where in the primer do those mismatches fall?

## What the program will do

Given a primer pair and a genome, predict every PCR product the pair would
amplify. A k-mer index locates candidate binding sites on both strands,
mismatches are allowed and scored by distance from the 3' end, and sites
facing each other at a plausible distance are paired into products. A harness
runs published primer sets against each variant genome and reports which
targets survive.

## Checks

- CDC N1, N2, and N3 must each give one product on the reference genome, at
  the amplicon sizes published in Lu et al. (2020).
- Primer sites planted by hand in short synthetic sequences, including
  reverse-strand and single-mismatch cases.
- Nearest-neighbor melting temperatures compared against Primer3 at matched
  salt and oligo concentrations.

## Data
see data folder
## Running it
first model does not account for annealing temperature, but includes pseudoprimers, solve-by-hand cases, and checks CDC primers
second model does attempt to account for annealing temperature but still needs some work
**first website applies our primer specificity model in reverse: it finds possible primers for a given genome and target gene while accounting for specificity, mismatch allowance, GC content, 3' distance of mismatches, mismatch severity, etc.
Limitation: it does not yet include annealing temperature, but the website's first version was tested against reverse-engineered test cases, CDC primers, and we generated some possible cases to see how the primers our website produces compares with NCBI BLAST **
