# Project Diary: Do early SARS-CoV-2 qPCR primers still match later variants?

02-120 Programming for Scientists · Isabel Guo

One log per version of the model, see also `modification documentation.md`.
The numbers in brackets, like **(#2)**, point to its matching numbered items.

---

## Log 1 · 2026-10-29 · Version 1: First model, factor: specificity

**Commit:** `8e3f1c3`

### What the model does, in plain words

A PCR test works by two short pieces of DNA, called **primers**, sticking to
the virus's genome and copying the stretch between them. If the virus
mutates where a primer sticks, the primer may not stick anymore, and the test
can miss the virus.

This first model asks one simple question: **"Where would these two primers
stick, and would they copy anything?"**

1. **Make a lookup table of the genome.** Like the index at the back of a
   book, it records where every short 5-letter word appears, so the program
   doesn't have to read the whole genome for every primer.
2. **Find where each primer sticks.** It looks on both strands of the DNA
   (DNA is double-stranded, and a primer can stick to either one) and lets
   you choose how many wrong letters (mismatches) are still OK. That setting
   is marked with a star banner so it's easy to find **(#1)**.
3. **Say where the mistakes are.** Each mismatch is counted from the primer's
   **3' end**, the end where copying starts. A mistake there matters most.
4. **Pair the primers into products.** A copy (product) is only made when
   one primer points right, the other points left, and they face each other
   at a sensible distance.

It was tested on **made-up primers and a made-up DNA sequence**, so every
answer could be worked out by hand first (7 test cases).

**Not included yet:** melting temperature, primer dimers, real genomes.

### Core functions
| Function | What it does |
|---|---|
| `reverse_complement` | Gives the other DNA strand, so primers can be found on both strands |
| `build_kmer_index` | The lookup table: where every 5-letter word appears |
| `find_binding_sites` | Finds every spot a primer sticks within the mismatch allowance |
| `distances_from_3prime` | Counts each mismatch's distance from the 3' end |
| `pair_into_products` | Pairs facing sites into PCR products |
| `predict_products` | Runs the whole model in one call |

---

## Log 2 · 2026-10-29 · Version 2: Handling bad sequencing data

**Commits:** `81d909e`, `504d494`

### What I did
- Stopped counting unknown bases (N) as mutations. Before, missing data
  could make a primer look lost. Each primer now gets one of three outcomes:
  binds / mismatched / not assessable **(#2)**.
- Made the model find primer sites that run off either end of a sequence,
  and mark them "not assessable" **(#3)**.
- Added a safeguard so a long run of Ns can't make every primer
  "not assessable" **(#4)**.
- Wrote `test_cases_specificity.txt`, which explains every test case, and
  added 4 new tests (Cases 7-10) **(#5)**.
- Renamed the bad-data note in the code to "HANDLING BAD SEQUENCING DATA"
  **(#6)**.

### New core functions
| Function | What it does |
|---|---|
| `compare_site` | Splits differences into real mismatches vs unknown bases |
| `candidate_starts` | Also checks spots near Ns and at the sequence ends |
| `primer_outcome` / `assay_outcome` | Give the three outcomes for each primer and the pair |

---

## Log 3 · 2026-10-30 · Version 3: Explaining failures

**Commit:** `fcf677c`

### What I did
- Fixed a bug: if one primer can't bind anywhere, the answer is now
  "no product", not "not assessable" **(#7)**.
- Added a warning for mismatches near the 3' end **(#8)**.
- Added a probe check. A qPCR signal needs the probe to bind inside the
  product; if it doesn't, the outcome is "no signal" **(#9)**.
- Allowed degenerate bases in primers, e.g. Y = C or T, used by the CDC N3
  probe **(#10)**.
- Showed each mismatch as a base change, e.g. "primer C, genome G" **(#11)**.
- Added a near-miss report: where a failed primer almost bound, including
  one insertion or deletion **(#12)**.
- Printed positions counting from 1, like the CDC table **(#13)**.
- Added 10 new tests (Cases 11-19) **(#14)**.
- While testing, found that the near-miss search can't always find a
  near-miss **(#15)**, and that on random DNA it found fake matches. I added
  stricter rules so it no longer does **(#16, #17)**.

### New core functions
| Function | What it does |
|---|---|
| `run_model` | Replaces `predict_products`; returns sites, products, outcomes, near-misses, warnings |
| `expand_degenerate` / `base_matches` | Handle degenerate primer letters |
| `describe_changes` | Writes each mismatch as a base change |
| `find_near_miss` | Finds the closest site for a failed primer |
| `gapped_pairs` / `placements` | Line a primer up across an insertion or deletion |
| `probe_inside` / `product_status` | Check the probe and label each product |

---

## Log 4 · 2026-10-30 · Version 4: Smarter logic and the first real data

**Commit:** `8af20dc`

### What I did
- Made the model say "binds" when missing data can't change the answer. If a
  site still fits even when every N counts as a mismatch, it binds **(#18)**.
- Made the model read ambiguity codes in the genome, e.g. genome R (A or G)
  can never match a primer C **(#19)**.
- Graded mismatches as severe or mild, using published PCR studies. The 3'
  warning is now graded high / low / possible **(#20)**.
- Let a 1-base gap near the primer's 5' end still bind **(#21)**.
- Ran the model on real data for the first time. CDC N1, N2 and N3 on the
  reference genome each give exactly one product in the right place, and
  every probe binds. The human RNase P control gives none **(#22)**.
- Found an open question: N1 comes out at 72 bp, but the primer table says
  73 bp, even though its own coordinates span 72 **(#23)**.
- Proved the fast search is complete: it matched a slow check of every
  position in 600 searches. I also saved the random-DNA check as a test
  **(#24, #25)**.
- Split the model and the tests into separate files **(#26)**.
- Fixed small wording in the report ("1 base", not "1 bases") **(#27)**.

### New core functions
| Function | What it does |
|---|---|
| `genome_vs_primer` | Decides match, mismatch or unknown for any genome letter |
| `mismatch_severity` | Grades a mismatch severe or mild |
| `three_prime_risk` | Grades the 3' end: high, low or possible |
| `find_indel_sites` | Finds sites that bind through a small 5' gap |
| `read_fasta` | Reads a real genome file |
| `check_viral_assay` (in `cdc_reference_check.py`) | Checks one CDC assay on the reference genome |
| `brute_force_sites` (in the test file) | The slow every-position search used to prove the fast one |

---

## Log 5 · 2026-10-30 · Documentation (no model change)

**Commits:** `bfaff5d`, then the restore below

### What I did
- Started this diary.
- Turned the modification documentation from a `.txt` into a `.md` file, so
  it displays with headings on GitHub **(#28)**.
- Lined this diary up with the modification documentation so the two can be
  read side by side **(#29)**.
- Added a "Things to work on for next time" list, and set the rule that both
  documents get updated with every commit **(#30)**.
- My editor saved an older copy of this diary over the new one; I restored
  the full version and kept my shorter intro **(#31)**.

### Things to work on for next time
- Settle the N1 size question (72 vs 73 bp) by checking Lu et al. (2020)
  **(#23)**.
- Run the CDC primers on the Alpha, Delta and Omicron genomes, which is the
  real question of the project.
- Add the positive control from the proposal: a primer across Alpha's spike
  69-70 deletion must fail on the Alpha genome.
- Start the next factor: melting temperature with the nearest-neighbor model,
  checked against the Primer3 values in `data/primers/`.
- Make the mismatch allowance something you can set when running the program,
  not only in the code.

---

## Where things stand
Still open after version 4 (details at the end of
`modification documentation.md`): the N1 size question, a few search limits,
and running the model on the Alpha, Delta and Omicron genomes.
