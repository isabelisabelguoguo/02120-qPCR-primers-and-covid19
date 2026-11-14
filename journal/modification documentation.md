# Modification Documentation: First Model, Factor: Specificity

(From Log 9 on, this file also covers the second model, factor: annealing temperature.)

Every problem found in the specificity model and how it was fixed, in the
order it happened, with the technical details.

**Reading alongside `diary.md`:** both files use the same logs in the same
order. The diary says what was done in plain words; this file gives the
problem, the fix, and the functions and settings involved. Items are
numbered, and the diary points to them as **(#n)**.

**Files** (all in `journal/`):

| File | Contents |
|---|---|
| `first_model_factor_specificity.py` | the model |
| `test_first_model_factor_specificity.py` | hand-solvable tests + search checks |
| `cdc_reference_check.py` | CDC panel on the real reference genome |
| `test_cases_specificity.txt` | every test case and its outcome |
| `second_model_annealing_temp.py` | second model: Tm and binding at the annealing temperature (Log 9) |
| `test_second_model_annealing_temp.py` | hand-solvable Tm tests + Primer3 check |
| `variant_check_temperature.py` | CDC N1/N2/N3 Tm on every variant |
| `primer_design.py` | designs primers + probe for a target gene, checked with the first model (Log 11) |
| `test_primer_design.py` | fake-genome tests for the designer |
| `test_cases_primer_design.txt` | every designer test case and its outcome, plus the CDC check |
| `troubleshooting.txt` | one line per problem overcome, with the functions each fix changed |
| `test_cases_web_specificity.txt` | 5 website test cases (GC levels and off-target copies) with sequences |
| `check_website_port.py` | checks the website's JavaScript designer gives the same results as the Python |
| `../website/` | the website: `index.html`, `primer_design.js`, `example_data.js` (Log 13) |
| `diary.md` | plain-language log of each version |

---

## Log 1 · 2026-10-29 · Version 1: First model, factor: specificity

**Commit:** `8e3f1c3`

**Technical overview.** Specificity only: no melting temperature, no real data.
- `build_kmer_index` maps every k-mer (k = 5) of the target to its positions.
- `find_binding_sites` seeds by the pigeonhole principle: a primer split into
  (max_mismatches + 1) pieces has at least one exact piece at any true site.
  Each seed hit is then checked over the full primer, on both strands (the
  "-" strand by searching the primer's reverse complement).
- `distances_from_3prime` reports mismatches as distance from the 3' end
  (1 = 3'-terminal base); the 3' end is the right end of a "+" site and the
  left end of a "-" site.
- `pair_into_products` pairs each "+" site with each downstream "-" site
  within `min_length`–`max_length` (default 20–1000 bp).
- Tests: 7 hand-solvable cases on fake 16-base primers and a poly-A-padded
  fake target.

**#1. Problem:** stating the seeding limit:
(primer length) // (max_mismatches + 1) must be ≥ k.

---

## Log 2 · 2026-10-29 · Version 2: Handling bad sequencing data

**Commits:** `81d909e`, `504d494`

**Technical overview.** Read alignment (Salmon, STAR, TopHat) is not needed:
the genomes are finished NCBI consensus sequences. Poor sequencing still
leaves N bases and missing ends in them. GC content was rejected as a quality
signal (FASTA has no quality scores; SARS-CoV-2 has real AT-rich regions).

**#2. Problem:** an N in the genome was counted as a mismatch, so missing data
looked like a mutation and a primer could be wrongly reported as lost.
**Fix:** `compare_site` replaces `mismatch_offsets` and returns mismatches and
unknown positions separately; unknowns never count toward the allowance. New
`primer_outcome` / `assay_outcome` give three outcomes: binds, mismatched,
not assessable.

**#3. Problem:** primer sites that ran off the start or end of a sequence were
never found.
**Fix:** `candidate_starts` adds every placement hanging off either end
(positions outside the sequence count as unknown) and every placement
overlapping an N, since k-mer seeds can't land on Ns. Such sites are "not
assessable" with the reason "runs off the start/end of the sequence".

**#4. Problem** (found while fixing #2): a long run of Ns would make every
primer "not assessable".
**Fix:** a not-assessable site needs at least `MIN_KNOWN_BASES` = 10 known
bases within the allowance.

**#5. Problem:** no plain-language record of what each test checks.
**Fix:** `test_cases_specificity.txt`, with every case's exact target,
outcome and meaning. Tests added: Cases 7-10 (Ns in a site, cut end, cut
start, real mismatches next to missing data).

**#6. Problem:** the bad-data note in the code was titled "NOTE TO READERS".
**Fix:** renamed "HANDLING BAD SEQUENCING DATA".

---

## Log 3 · 2026-10-30 · Version 3: Explaining failures

**Commit:** `fcf677c`

**Technical overview.** `run_model` replaces `predict_products` and returns
one result: `sites`, `products` (each with a `status`), `outcomes`,
`near_misses`, `assay`, `warnings`.

**#7. Problem:** logic bug. In Case 10 the reverse primer could not bind
anywhere, yet the assay was "not assessable" only because of an N in the
forward site.
**Fix:** `pair_into_products` now pairs "binds" and "not assessable" sites
into possible products, and `product_status` labels each one. The assay is
"not assessable" only when a possible product depends on missing data;
otherwise "no product".

**#8. Problem:** a 3' mismatch counted the same as a 5' one, although 3'
mismatches are the ones that block extension.
**Fix:** sites with a mismatch within `THREE_PRIME_WINDOW` = 3 bases of the
3' end are flagged, and their products produce a warning (new
`3' END RISK WINDOW` banner).

**#9. Problem:** the qPCR probe was ignored, although a product gives no
signal unless the probe binds inside it.
**Fix:** optional `probe` argument; `probe_inside` checks for probe sites
fully inside each product. Primers binding without the probe give the new
assay outcome "no signal". Assay priority: product > not assessable >
no signal > no product.

**#10. Problem:** degenerate primer bases (e.g. the Y in the CDC N3 probe)
were counted as mismatches and could not be used as seeds.
**Fix:** an `IUPAC` table; `base_matches` accepts any base a letter stands
for, and `expand_degenerate` / `seed_hits` look up every plain k-mer a
degenerate seed stands for. `clean_oligo` rejects non-DNA letters.

**#11. Problem:** mismatches gave only a position, not the base change.
**Fix:** `describe_changes` reports "primer C, genome G", both bases in the
primer's 5'→3' direction (complemented for "-" sites).

**#12. Problem:** a failed primer gave no clue where it almost bound, and
insertions/deletions (e.g. Alpha's spike Δ69–70) looked like many mismatches.
**Fix:** `find_near_miss` runs for every "mismatched" oligo. It anchors on
every exact k-mer the primer shares with the genome, then `placements` /
`gapped_pairs` try no gap or one gap of up to `MAX_INDEL` = 6 bases.
Score = substitutions + 1 per gap (one event). The best placement is reported
with its changes and gap, or "none close".

**#13. Problem:** positions were 0-based, unlike the CDC table and genome
browsers.
**Fix:** stored 0-based, printed 1-based by `print_report`.

**#14. Problem:** untested edge cases.
**Fix:** Cases 11-19 (self-pairing, overlapping sites, allowance too high →
`ValueError`, palindromic primer, degenerate base, probe binds / fails,
insertion, deletion).

### Found while building version 3

**#15. Problem:** near-miss anchoring needs one exact 5-base stretch. In
Case 10 three mismatches 4 bases apart leave none, so nothing was found.
**Fix:** documented as a limit; Case 10 expects "none close". A 16-base primer
is always found with ≤ 2 events, a 20-base primer with ≤ 3.

**#16. Problem:** on a random 30,000-base sequence, a CDC primer got a fake
"6-base deletion + 1 mismatch" near-miss: a long deletion leaves few bases to
compare, so chance matches are easy.
**Fix:** `MIN_MATCHED_BASES` = 14 primer bases must match exactly.

**#17. Problem:** chance matches still appeared (a deletion + 2–3 mismatches),
because a gap can sit in many places.
**Fix:** `NEAR_MISS_MAX_EVENTS_WITH_INDEL` = 2 (the gap + ≤ 1 substitution);
`NEAR_MISS_MAX_EVENTS` = 4 without a gap. Result: 0 false near-misses on 10
random genomes; a planted 6-base deletion in a 20-base primer is still found.
A 16-base primer can therefore show deletions of ≤ 2 bases: Case 18 changed
from a 4-base to a 2-base deletion, and Case 18b (6 bases, refused) was added.

---

## Log 4 · 2026-10-30 · Version 4: Smarter logic and the first real data

**Commit:** `8af20dc`

**Technical overview.** Four logic fixes, the first real-data check, tests
that verify the search itself, and the model and tests split into two files.

**#18. Problem:** "not assessable" was too cautious. With 2 mismatches
allowed, a site with 0 mismatches + 1 N binds whatever the N is, but was
reported "not assessable".
**Fix:** in `find_binding_sites`, if mismatches + unknowns ≤ max_mismatches,
the outcome is "binds" with the reason "binds even if every unknown base is a
mismatch". Case 7b added; in Case 10 the forward primer now binds.

**#19. Problem:** genome ambiguity codes were all treated like N. A genome R
(A or G) can never match a primer C.
**Fix:** `genome_vs_primer` returns match (every option accepted), mismatch
(no option accepted) or unknown (some accepted); `compare_site` uses it.
Case 20 added.

**#20. Problem:** mismatches were all-or-nothing; a mild G–T counted the same
as a severe C–C.
**Fix:** `mismatch_severity` grades primer–template pairs A–A, G–A, A–G and
C–C as severe and all others as mild (Stadhouders et al. 2010; Kwok et al.
1990); degenerate letters are severe if any option is. `three_prime_risk`
grades each site high (severe in window) / low (mild only) / possible
(unknown in window); `worst_risk` combines a product's two sites. The
allowance still counts each mismatch as 1, keeping the seeding guarantee.
Case 3b updated, Case 21 added.

**#21. Problem:** any insertion or deletion meant failure, though a 1-base gap
near the 5' end usually still extends.
**Fix:** `find_indel_sites` adds "binds" sites with one gap of
≤ `MAX_BINDING_INDEL` = 1 base, more than half the primer length from the 3'
end. The gap counts as 1 mismatch, and placements overlapping an existing
ungapped site are skipped. Shared helpers `anchors_for` and `describe_indel`
were factored out of `find_near_miss`. Case 22 added.

**#22. Problem:** the model had never been run on real data.
**Fix:** `cdc_reference_check.py` (`read_assays`, `check_viral_assay`,
`check_negative_control`) runs N1, N2, N3 with probes and the RNase P control
on NC_045512.2 (`read_fasta`). Allowance 2. Results: one product per N assay
at the verified positions, 0 primer mismatches, every probe binds (including
N3's degenerate probe), and no RP product.

**#23. Problem** (found by #22): N1 comes out at 72 bp, but the table lists
73 bp. The table's own published (28303–28374) and verified (28287–28358)
coordinates both span 72 bases.
**Status: RESOLVED in #36** (the table's 73 was wrong; 72 is correct).

**#24. Problem:** nothing proved the fast k-mer search finds every site.
**Fix:** Case 23 compares plain sites with `brute_force_sites` (every start
position, both strands, same classification rules) on 200 random 300-base
targets × allowances 0–2, with mutated primer copies, N/Y bases and cut ends.
Identical on all 600 searches (1,361 sites).

**#25. Problem:** the random-genome check from #17 was only run by hand.
**Fix:** saved as Case 24 (CDC N1 primers + probe on a random 30,000-base
sequence: no sites, no near-misses).

**#26. Problem:** the model and tests shared one ~1,000-line file.
**Fix:** tests moved to `test_first_model_factor_specificity.py`, which
imports the model. `read_fasta` added to the model for real genomes.

**#27. Problem:** the report read "deletion of 1 bases (primer bases 15-15 …)".
**Fix:** `format_indel` now writes "deletion of 1 base (primer base 15 from 3' …)".

---

## Log 5 · 2026-10-30 · Documentation (no model change)

**Commits:** `bfaff5d`, `526c2a5` (the restore in #31), `1a78e66` (#32), `289fb03` (the session log)

**#28. Problem:** this file was plain text (`.txt`), so its headings and labels
showed as raw characters on GitHub.
**Fix:** converted to Markdown (`.md`), with references in the model header
and `diary.md` updated. `test_cases_specificity.txt` stays `.txt` because its
`^` markers must line up under DNA letters.

**#29. Problem:** this file and `diary.md` had different structures, so they
were hard to read side by side.
**Fix:** both use the same logs, dates, titles and order; diary bullets point
to items here as **(#n)**.

**#30. Problem:** no record of next steps, and the documents could fall behind
what is on GitHub.
**Fix:** the latest diary log ends with "Things to work on for next time"
(added at the end of each session). From now on, every commit and push also
updates `diary.md` and this file. A log's commit ID is filled in at the next
update.

**#31. Problem:** the editor saved an older copy of `diary.md` over the
committed one, deleting Log 5, the next-time list and every **(#n)** link,
and pointing back to the removed `.txt` file.
**Fix:** restored the committed version with `git checkout`, kept the user's
shorter intro (pointing to the `.md` file), and filled in Log 5's commit ID.

**#32. Problem:** the README title ("02120 project1") did not say what the
project is about, and #1 here repeated the `USER SETTING` banner details.
**Fix:** renamed the README title to "02120 qPCR primers and covid-19" and
shortened #1 to the seeding limit.

---

## Log 6 · 2026-10-30 · First variant results (quick first pass)

**Commits:** `f0b44dd`, `761fe07` (#35), `5834b32` (#36), `7892816` (#37 and the session logs)

**#33. Problem:** the model had only been run on the reference genome, not on
the variants the project asks about.
**Fix:** new `variant_check.py` reuses `read_assays`, `PRIMER_TABLE`,
`REFERENCE` and `MAX_MISMATCHES` from `cdc_reference_check.py` and calls
`run_model` (allowance 2, probe included) for N1/N2/N3 on the reference and
every `data/variants/*.fasta`. For each oligo, `fewest_mismatches` reports
the lowest mismatch count among its binding sites; the table also shows the
assay outcome, product size and 3' risk. Result: every assay gives a product
on Alpha, Delta and Omicron (N1 72, N2 67, N3 72 bp, as on the reference);
the only mismatch is 1 in the N1 probe on Omicron; no 3' risk.
Limit: one isolate per variant.
**⚠ Flag (#42):** the genome labelled "Alpha" here is really Beta (B.1.351),
so the "Alpha" row of these results is a Beta result.

**#34. Problem:** the genome files had mangled names
(`data:reference:NC_045512.2.fasta.`, `OmicronBA.1_OL672836.1.`) and the
variants had no `.fasta` extension.
**Fix:** renamed with `git mv` to `NC_045512.2.fasta`,
`Alpha_B.1.1.7_MZ068155.fasta`, `Delta_B.1.617.2_MZ555945.fasta` and
`Omicron_BA.1_OL672836.1.fasta`. `cdc_reference_check.py` finds the reference
with a glob, so it needed no change (same results).
**⚠ Flag (#42):** `Alpha_B.1.1.7_MZ068155.fasta` was mislabelled; it is a
Beta genome and is now `Beta_B.1.351_MZ068155.fasta`.

**#35. Problem:** `.DS_Store` files (3 committed) kept showing as modified,
even though `.gitignore` lists `.DS_Store` (it only stops new files), and
`journal/__pycache__/` appeared after running the scripts.
**Fix:** added `__pycache__/` to `.gitignore`; the committed `.DS_Store` files were removed from the repo with
`git rm --cached` (kept on disk).

**#36. Problem:** the N1 size question from #23 (72 bp computed vs 73 bp in
the table) was still open, so `cdc_reference_check.py` reported 1 FAIL.
**Fix:** the table's own coordinates (published 28303–28374, verified
28287–28358) both give 72 bases, the reference genome gives 72, and the JRC
SARS-CoV-2 diagnostics database lists the CDC N1 product as 72 bp. The 73 was
a data-entry error: `pub_amplicon_bp` for N1 in
`data/primers/CDC primer sequences.tsv` changed from 73 to 72. All CDC
reference checks now pass.

**#37. Problem:** #33 found 1 mismatch in the N1 probe on Omicron but did
not say where it is or whether it is a real lineage mutation.
**Finding:** a one-off comparison of the reference and
`Omicron_BA.1_OL672836.1.fasta` at the N1 probe site (28309–28332) puts it at
**28311, C→T**, base 3 from the probe's 5' end. This is the known Omicron
mutation **N:P13L** (nucleocapsid codon 13, proline → leucine), so it is a
lineage change, not an isolate quirk. No code changed.
Note (from reading `find_binding_sites` and `product_status`): mismatch
severity and `three_prime_risk` are labels only. The binds/product decision
uses the mismatch count alone, so a severe 3' mismatch still gives
"product" with a `[3' RISK: HIGH]` warning.

---

## Log 7 · 2026-11-09 · Realistic test primers and k = 6

**Commit:** `0e0e449`

**Note:** most primers used in real (biological) PCR and qPCR are 18-30
bases long: long enough to be unique in a genome and to anneal at about
55-65 °C, short enough to bind quickly and cheaply. The CDC SARS-CoV-2
primers and probes are 18-24 bases. The tests and k = 6 are built around
this range.

**#38. Problem:** every hand-solvable test used 16-base fake primers, shorter
than real qPCR primers (18-30 bases; the CDC ones are 18-24). The 16 was never
justified; it was simply the longest length that fit k = 5 with 2 mismatches
(16 // 3 = 5).
**Fix:** `test_first_model_factor_specificity.py` rebuilt on two fake sets:
set A (`FORWARD_A` 20, `REVERSE_A` 24, `PROBE_A` 26 bases; 70 bp amplicon,
used by most cases) and set B (`FORWARD_B` 18, `REVERSE_B` 30 bases; 73 bp
amplicon, for length cases), plus an 18-base `PALINDROME` and a 20-base
`OVERLAP_REVERSE`. Every expected answer was redone by hand. Changes beyond
the numbers:
- Case 1b (new): exact match with set B, 73 bp.
- Cases 3, 3c, 21a: the planted 3' and 5' mismatches are now A-C (mild).
- Case 10: 4 mismatches 5 apart are needed to hide a 24-base site from the
  near-miss search (3 always leave a long enough stretch).
- Case 13: the error now comes from the 18-base primer (18 // 4 = 4).
- Case 18b now uses the 18-base primer (6-base deletion leaves 12 < 14
  matched bases, "none close"); new Case 18c finds the same 6-base deletion in
  the 30-base primer (24 bases still match). Whether a deletion is recognised
  depends on primer length.
- Case 22b: the 3'-half deletion moved to 4 from the 3' end; at 2 or 3 the
  neighbouring Cs made it look like 0 or 1 mismatch.
- Case 23 now checks primers of 18, 20, 24 and 30 bases (600 searches agree
  with the brute-force search, 1129 sites).
`test_cases_specificity.txt` rewritten with the new targets and markers;
the model's comments no longer use 16-base examples.

**#39. Problem:** k = 5 was a leftover of the 16-base test primers. With real
primers of 18 bases or more, a longer seed is allowed and much faster: a
random 6-mer hits about 7 places in a 30,000-base genome, a 5-mer about 29.
**Fix:** `run_model` default `k=6`, with the reason written above
`build_kmer_index`. The seeding rule (length // (max_mismatches + 1) >= k)
now needs at least 18 bases with the default 2 mismatches; every CDC oligo
qualifies (shortest: N2 reverse, 18). The tests take ~10 s instead of
~18 s. `cdc_reference_check.py` and `variant_check.py` give identical results
with k = 5 and k = 6.

**#40. Problem:** the near-miss and small-gap searches also start from exact
k-base stretches, so k = 6 weakens them: a near-miss is now guaranteed with
up to 2 events for 18- and 20-base primers (3 for 24, 4 for 30), down from 3
for a 20-base primer with k = 5.
**Fix:** documented as a limit (Section 5 comment) and shown by new Case
10b: 3 mismatches cut the 20-base forward site into exact stretches of
5, 5, 5 and 2 bases; with k = 6 the near-miss is "none close", with `k=5` it
is found (3 events, 3, 9 and 15 from the 3' end).

---

## Log 8 · 2026-11-09 · 3' risk changes the result; a mislabelled genome

**Commit:** `e6793bf`

**#41. Problem:** a severe mismatch at a primer's 3' end only printed a
warning; the product still counted as a plain "product" (#37), so a pair
likely to amplify poorly looked as good as a perfect one.
**Fix:** new product status **"at risk"**: both primers bind (and the probe,
if given), but a primer site has `three_prime_risk` "high" (a severe mismatch
within `THREE_PRIME_WINDOW` = 3 bases of the 3' end). `product_status`
returns it after the "no signal" and "not assessable" checks. The new
`STATUS_ORDER` = product > at risk > not assessable > no signal > no product
is used by `assay_outcome` and is meant to rank primer pairs on the future
website (at-risk pairs flagged and suggested last). LOW and POSSIBLE risk do
not change the status. The probe's 3' risk is ignored: a TaqMan probe is
never extended. Warnings now cover "at risk" products too;
`cdc_reference_check.py` and `variant_check.py` count "at risk" products as
products. Case 3b now expects "at risk"; new Case 3d (two copies, one with
the severe 3' mismatch) gives two "at risk" products and one "product", and
the assay is "product". Cases 21a/21b check that LOW and POSSIBLE stay
"product". The CDC and variant results do not change (no HIGH 3' risk).

**#42. Problem:** `Alpha_B.1.1.7_MZ068155.fasta` is not Alpha. At the spike
signature positions it carries all of Beta's (B.1.351) changes (D80A
A21801C, D215G, K417N, E484K, N501Y, Δ242-244) and none of Alpha's (ΔH69/V70,
ΔY144, A570D, P681H, T716I). The NCBI record (MZ068155.1,
SARS-CoV-2/human/PHL/COVID53909/2021) names no lineage; the label came from
the file name. Found while testing a spike Δ69-70 assay, which should fail on
Alpha but gave a product on this genome.
**Fix:** renamed with `git mv` to `Beta_B.1.351_MZ068155.fasta`
(`variant_check.py` now labels it "Beta"; its header comment flags the
change). Every "Alpha" result in Log 6 (#33) is a Beta result; flagged there.
Still to do: confirm the lineage with Nextclade, and add a real Alpha genome
for the Δ69-70 positive control.

**Note (no code change):** the session log for this commit
(`af9073fb-….jsonl`) includes a 15-question quiz on the model plus review
rounds. These rules were confirmed against the code:
- 3' risk uses `from_3prime <= THREE_PRIME_WINDOW`, so 1–3 bases from the
  3' end are inside the window.
- "At risk" needs both a severe pair and a position inside the window.
- Status is decided in the order mismatch count, then window, then severity.
- The 1-base indel rule (`from_3prime[0] > length // 2`) is separate from
  3' risk.
- "No signal" means the product is made, but the probe is mismatched.

---

## Log 9 · 2026-11-10 · Second model: annealing temperature

**Commit:** `b585a97`

**#43. Problem:** binding was yes/no. A mismatch that the first model
allows can still weaken the duplex enough that the oligo falls off at the
annealing temperature, and the first model cannot see that (Still open, #37).
**Fix:** new `second_model_annealing_temp.py`. It runs `run_model`, then for
every site builds the duplex (oligo 5'->3' over the template base it really
faces) and sums nearest-neighbor dH/dS (`duplex_energy`):
- matched stacks from `santalucia_1998_nn.tsv`, plus initiation at each end
  pair and the symmetry term for palindromes
- salt correction as in Primer3:
  `Na_eq = mono + 120*sqrt(Mg - dNTP)` (von Ahsen 2001), then
  `dS += 0.368*(N-1)*ln(Na_eq)`
- `Tm = 1000*dH / (dS + R*ln(C/4)) - 273.15`

Checked against `primer3_reference_tm.tsv`: all 11 CDC oligos match to
0.01 C (Case T5). Conditions and `ANNEALING_TEMP_C` = 55 C (CDC: 95 C 3 s,
55 C 30 s, 45 cycles) are a USER SETTINGS box.

**#44. Problem:** the parameter table had only the 10 matched stacks, so a
variant's real site (with its mismatch) could not be scored.
**Fix:** added two tables to `data/thermodynamics/`, with citations in their
headers:
- `allawi_santalucia_internal_mismatch_nn.tsv`: 51 stacks; Allawi &
  SantaLucia 1997-98, Peyret 1999
- `santalucia_peyret_terminal_mismatch_nn.tsv`: 48 stacks; SantaLucia &
  Peyret 2001

Both were copied from Biopython 1.88's `DNA_IMM1`/`DNA_TMM1`.

A terminal mismatch takes its stack from the terminal table and moves the
initiation term to the next pair in. Two mismatches side by side, indel
sites, unknown genome bases and mismatches at a degenerate oligo letter are
"not computed" with a reason, and the first model's status stands.

Checked on 3000 random 18-30 base duplexes:
- internal mismatches match Biopython `Tm_NN` to 0.004 C
- terminal ones differ by up to 0.4 C, because Biopython puts the initiation
  term on the mismatched base

Primer3's `ntthal` agrees on perfect duplexes (0.06 C) and usually on single
mismatches (SD 1.5 C), but gives 14-16 C drops for some single A-A/C-C
mismatches, far beyond measured values. It is not used as a reference.

**#45. Problem:** Tm alone does not say whether an oligo is still bound at
the annealing temperature, and in PCR the oligo is in excess over the
template.
**Fix:** `fraction_bound` = K*C / (1 + K*C) with K = exp(-dG/RT) at the
annealing temperature. It is 50% at 1000*dH / (dS + R*ln C), about 2 C above
the Primer3 Tm. Each site reports Tm, dTm (vs a perfect match), margin
(Tm - Ta) and fraction bound.

The new `status_at_annealing` starts from the first model's status and can
only make it worse:
- **"weak binding"** if a primer site is bound less than
  `MIN_FRACTION_BOUND` = 0.5
- **"no signal"** if the primers hold but no probe site inside the product
  does

`STATUS_ORDER_AT_TA` = product > at risk > weak binding > not assessable >
no signal > no product.

Tests (`test_second_model_annealing_temp.py`):
- T1-T4: hand sums; T3 is an internal A-G at -2.59 C
- T4: a severe 3'-terminal G-A drops Tm by only 0.47 C. It is 3' risk HIGH
  in the first model, so the two models answer different questions.
- T4b: tandem mismatches are not computed
- T4c: all 240 internal single mismatches lower Tm
- T4d: fraction bound
- T6: 68 C gives "weak binding"
- T7: a probe with 2 mild mismatches falls from 69.1 to 60.1 C and gives
  "no signal" at 64 C
- T8: an indel site is not computed

**Result** (`variant_check_temperature.py`): on Beta, Delta and Omicron every
CDC oligo keeps its reference Tm except the Omicron N1 probe (C28311T, #37):
71.2 -> 64.3 C (dTm -6.8 C). It is still ~9 C above 55 C and ~100% bound, so
all three assays give "product" at 55 C.

---

## Log 10 · 2026-11-11 · Check-in practice (no model change)

**Commits:** `6da5ebe`, then `ade42ee` (#46, the primer image and the
session logs)

The model itself did not change. This log was practice for explaining the
first model at the check-in, then one change to the demo (#46).

**#46. Problem:** the demo printed only the final report, so I couldn't see
which functions run, in what order, to get from the genome to the result.
**Fix:** Section 11 (the `if __name__ == "__main__":` demo) now runs
`run_model` once on the 120-base sample genome with a mutation 2 bases from
the forward primer's 3' end. Every function in the file is wrapped (`traced`)
so the run prints `running <function>`, indented under its caller. Repeat
calls are grouped into one line with a count, e.g. `genome_vs_primer (x2577)`,
in order of first call; without grouping the trace is about 12,000 lines.
The last line prints `return output = ` and the full dictionary `run_model`
returns. The demo runs only when the file is run directly, not when the
tests or `variant_check.py` import it. All tests still pass.

Also added to the diary: the NHGRI primer diagram
(`images/primer_binding_nhgri.jpg`) and what `strand +` / `strand −` mean.

**To do next time** (limits and next steps reviewed in the practice):
- **One genome per variant:** run each variant on several isolates (#33).
- **No real Alpha (B.1.1.7) genome:** the file saved as Alpha was Beta (#42).
- **Seed limit:** `candidate_starts` and the near-miss search need an exact
  6-base stretch (k = 6), so a site whose differences cut the primer into
  pieces all shorter than 6 bases is not found (#15, #40).
- **Yale spike Δ69-70 assay:** save its primers and probe in
  `data/primers/` and run them through both models.
- **Tm at indel sites:** add bulge-loop parameters so
  `second_model_annealing_temp.py` can score insertion/deletion sites (#44).
- **Run button:** VS Code's Run button quotes the path twice because the
  folder name `programming project` has a space, so the file never runs.
  Rename the folder (e.g. `programming_project`) and recreate `.venv`.

---

## Log 11 · 2026-11-13 · Primer design from the first model

**Commit:** `19b27c0`

**#47. Problem:** the models only check existing primers. For the planned
website, a user should be able to give a genome, a target gene and a
mismatch allowance and get primers back.
**Fix:** new `primer_design.py`, `design_primers(genome, target,
max_mismatches)`:
1. `clean_sequence` accepts plain or FASTA text (drops `>` lines,
   whitespace, lowercase); `find_target` finds the target on either strand.
2. `possible_designs`: every `PRIMER_LENGTH` = 20 window of the target with
   known bases only (`KNOWN_BASES`) is a forward primer; its reverse
   complement is a reverse primer. Pairs giving `PRODUCT_LENGTH` = 60-150 bp
   get a `PROBE_LENGTH` = 24 probe centred in the gap. Sorted by product
   length closest to `IDEAL_PRODUCT` = 100, then start.
3. `specific`: each oligo must have exactly one site from
   `find_binding_sites` on the whole genome (one shared k-mer index, results
   cached per oligo), then `run_model` must give exactly one product, status
   "product", at the designed coordinates.

`MAX_ALLOWED_MISMATCHES` = 2: with k = 6, `find_binding_sites` needs
length // (allowance + 1) >= 6, so 20-base primers allow at most 2 (#15,
#40). Output: 5'-forward-3', 3'-reverse-5' (written backwards to line up
under the genome) and the order sequence, binding positions (1-based),
lengths, product, probe position.

The first draft also filtered on Tm (58-62 C, probe 6-12 C above, via the
second model's `perfect_tm`), GC 40-60% and base runs. All of it was removed
so the designer uses first-model information only.

**#48. Problem:** the three designs returned were one pair shifted by one
base, so all three products overlapped.
**Fix:** a design is skipped if its product overlaps an already chosen one.

**#49. Problem:** with first-model information only, nearly every window
passes, including GC-rich ones that bind too strongly and mis-prime.
**Fix:** Section 2b, `gc_fraction` (known bases only, so N does not lower
it) and `gc_level`: ok <= 60% (`GC_MEDIUM_ABOVE`), medium <= 70%
(`GC_TOO_HIGH_ABOVE`), too high above that.
- `check_target_gc` runs first in `design_primers`, after the input checks
  and before `find_target`: too high raises `ValueError` with the GC value
  and the reason; medium adds a warning.
- `usable` excludes oligos above 70%. Medium-GC oligos are kept, listed in
  each design's `medium_gc`, and get one warning each.
- Sort key is now (number of medium-GC oligos, product distance from 100,
  start), so ok-GC designs come first.
- The command line prints `ERROR: ...` and exits 1 on a `ValueError`.

**#50. Problem:** `gc_percent` rounded to whole numbers, so a 69.6% target
printed as "70 %" while labelled medium.
**Fix:** one decimal place.

**#51. Tests** (`test_primer_design.py`, cases in
`test_cases_primer_design.txt`): a 4,400-base seeded random genome with a
400-base target at 2001-2400. Every returned design is checked against all
the rules and by an independent `run_model` call. 12 cases: basic, first
design by hand (F 2001, R 2081-2100, probe 2039, 100 bp), reverse-complement
target, FASTA input, N block, GC-rich patch, two copies (no design), a near
copy with every 10th base changed (no design at allowance 2, 3 at
allowance 1), 80% GC target (error), 80% GC target missing from the genome
(GC error first), 65.5% GC target (warning) and bad input. All pass.
The first run failed in Case 11: the test counted warnings containing
"design", and the target warning contains "designs". Fixed in the test
(`startswith("design ")`); the designer was unchanged.

**#52. Result: CDC N1/N2/N3** on NC_045512.2, target = published product
+/- 50 bases, allowance 2:
- The published oligos pass every designer rule: 1 site each, one product,
  probe inside, GC 40.0-58.3%. The N3 probe has a Y, which the designer
  never proposes (known bases only).
- The designer does not reproduce them. Its first choices are 100 bp
  products starting at the target's first base (e.g. N1: F 28237-28256,
  R 28317-28336) vs the published 67-72 bp. Causes: fixed 20/24-base
  lengths (CDC: 18-24 / 23-24) and `IDEAL_PRODUCT` = 100.

---

## Log 12 · 2026-11-13 · Primer design: lengths, the 40-60% GC rule, ranking

**Commits:** `3f629ea`, then `24f6899` (`troubleshooting.txt`, #58)

**#53. Problem:** every primer was 20 bases and every probe 24, so most CDC
oligos (18-24 / 23-24) could not even be among the options. A length range
alone would give 7 x 13 length combinations per position and millions of
designs, with nothing to choose between them.
**Fix:** `PRIMER_LENGTH` = (18, 24), `PROBE_LENGTH` = (18, 30). 18 is the
minimum for allowance 2 with k = 6 (18 // 3 = 6; #15, #40). `choose_length`
gives each anchor one window: the usable length with GC closest to
`IDEAL_GC` = 0.5, the shorter on a tie (strict `<` while lengths grow).
Forward primers and probes are anchored at their 5' end and grow right;
reverse primers are anchored at their 5' end, which is the window's right
end, and grow left (`all_windows`). `centred_probe` picks the probe whose
middle is closest to the middle of the gap. N gene: about 100,000 designs,
as before; the demo takes about 2 s.

**#54. Problem:** oligos were only limited from above (<= 70% GC), and
60-70% oligos were used with a warning. Stable binding needs 40-60%.
**Fix:** `OLIGO_GC` = (0.40, 0.60); `usable` now requires it, so the
per-oligo medium flag (`medium_gc`) was removed. Target bands in
`gc_level`/`check_target_gc`: > 70% error, > 60% "medium" warning, < 40%
(`TARGET_GC_LOW_BELOW`) "low" warning only, since SARS-CoV-2 is ~38% GC.
If `possible_designs` returns nothing, a warning says no 40-60% GC
stretches fit together into a product.

**#55. Problem:** with more than one design, the order was only product
length then position.
**Fix:** Section 3b, `match_penalty` = sum over forward, reverse and probe of
|GC% - 50| + |product - ideal| / 10; designs sorted by (penalty, start).
Each design prints its penalty. There is no "flagged" second key yet: with
every oligo at 40-60%, no single design carries a flag; target warnings
apply to all designs equally.

**#56. Problem:** the product length was fixed at 60-150 bp around 100.
**Fix:** `design_primers(..., ideal_product=100)` and a 4th command-line
argument. `product_range` = (max(`MIN_PRODUCT` = 54, ideal - 50), ideal +
50), so the default range is 54-150. Errors: below 54 bp, or longer than the
target. Warning outside `QPCR_PRODUCT` = 60-150 bp. `specific` passes
`max_length=max(1000, product)` to `run_model`, so products over 1000 bp are
still seen. Tried on the N gene: 70, 300 and 1000 bp give designs of
exactly that length in 1-3 s.

**#57. Tests and CDC result.** `test_primer_design.py` now has 15 cases (all
pass). Changed: length checks use the ranges; Case 2 is the ranking (hand
example 5 + 0 + 5 + 2 = 12, and penalties in order); Case 8's near copy
changes every 12th base, since an 18-base primer could cover only one change
in every 10th; Case 6 checks 40-60%. New: Case 12 (28% GC target: low
warning, 1 design), Case 13 (`ATTA` x 50: no design, GC named as the
reason), Case 14 (product 70 and 300 bp). `test_cases_primer_design.txt`
rewritten as version 2 with every case's real output.
CDC N1/N2/N3: all oligos pass the rules (40.0-58.3% GC). At 100 bp the N2
reverse primer chosen is `ATGCGCGACATTCCGAAGAA`, the CDC reverse plus 2 bases
at its 5' end. At 70 bp the products are 70 bp but at other positions:
many designs score 0, and the earliest wins the tie.

**#58. Problem:** the problems and fixes were only in this file, in date
order and in full detail, so there was no quick overview.
**Fix:** new `troubleshooting.txt`: one short bullet per problem overcome
("problem -> fix"), grouped by topic, each with a "Functions:" line naming
the functions it changed or is directly about (all names checked against
the code). Part 1 covers the first model (#1-#46, leaving out the
journal/repo items #5, #6, #28-#32, #34, #35); Part 2 covers
`primer_design.py` (#47-#57), with #52 marked not yet overcome. A first
version was a runnable `troubleshooting.py` that re-checked 15 fixes live;
it was replaced by the plain-text summary before being committed.

---

## Log 13 · 2026-11-13 · The website, and product length as a range

**Commit:** filled in at the next update

**#59. Problem:** the planned website has to run in a browser, which cannot
run Python.
**Fix:** `website/primer_design.js` translates `primer_design.py` and the
first-model functions it uses (`build_kmer_index`, `find_binding_sites`,
`candidate_starts`, `compare_site`, `find_indel_sites`, `pair_into_products`,
`probe_inside`, `product_status`, ...) function by function; each JS function
names the Python one it copies. Near-miss search and printing are left out.
One addition: `checkLetters` gives an error for non-DNA letters (the Python
raises a `KeyError`). N/ambiguous positions are found once per genome
instead of on every call (same answer, faster). `website/index.html`: genome
box (paste or open a FASTA file), target box, allowance (0-2), product
range; ranked design cards with an SVG map of the product, the reverse
primer shown 3'->5' with its order sequence, and copy buttons.
`website/example_data.js`: NC_045512.2 and its N gene for the example.
Published as a private artifact (https://claude.ai/artifact/XKd9SSGhpFHH1WXYemP3Bh).
`journal/check_website_port.py` runs both versions on the same inputs (via
macOS `osascript -l JavaScript`, so no Node) and compares designs,
positions, penalties, counts and warnings/errors: 35 of 35 identical.

**#60. Problem:** copying long sequences from the terminal didn't work.
**Fix:** buttons on the page load each test case and run it.
`test_cases_web_specificity.txt`: 5 cases (genome = 100 random + 200-base
target + 100 random, cases 2-4 + 50 random + a second copy):
1. good GC (46.0 %), no copy -> 2 designs
2. low GC (30.5 %), copy with every 6th base changed (3-4 mismatches per
   primer) -> warning, 1 design; the copy never counts
3. high-but-ok GC (66.0 %), every 12th base changed (1-2 mismatches) ->
   warning; a design at allowance 0 and 1, none at 2
4. good GC, exact copy -> no design
5. 89.0 % GC -> error

The same sequences are in `example_data.js` (`TEST_EXAMPLES`), and
`check_website_port.py` runs each at allowance 0, 1 and 2.

**#61. Problem:** the product length was one number with +/- 50 bp around
it; testing the website, a range was more natural.
**Fix:** `design_primers(..., min_product=60, max_product=150)` (CLI: 4th
and 5th arguments). `PRODUCT_RANGE` = (60, 150) replaces `IDEAL_PRODUCT` and
`PRODUCT_SPREAD`; `product_range` was removed. `match_penalty` now adds
(product - min_product) / 10, so shorter products rank higher. Errors:
shortest < 54 bp, longest < shortest, target shorter than the shortest.
Warning when the range goes outside 60-150. Same change in
`primer_design.js` and the page (two boxes). Tests: Case 2 (60 bp -> 0,
150 bp -> 9; the hand example is still 12), Case 14 (70-90 and 250-350),
Case 15 (53-150, 150-60, 59-base target); 15 cases pass.
`test_cases_primer_design.txt` rewritten as version 3. CDC check with the
default range: the designer's first choices are 60 bp products; N1's
(28275-28334) is next to the CDC N1 product (28287-28358).

**Process note (not a model problem):** `save_journal.py` reads stdin, so
inside a combined command it waited forever and swallowed a heredoc commit
message (no commit was made). It is now run as
`python3 .claude/hooks/save_journal.py < /dev/null`.

---

## Still open (known limits after version 4)
- Near-miss and small-gap searches need an exact 6-base stretch (#15, #40).
- Near-miss limits (#16, #17) were tuned on random sequence, not real genomes
  with repeats.
- Severity is a two-level rule (#20); real effects also depend on
  neighbouring bases and positions beyond the 3' window.
- Sites running off the start print negative positions (e.g. -3-48).
- Input is not cleaned of pasted spaces or line breaks.
- The k-mer index is rebuilt on every run; the test suite takes ~10 s.
- The mismatch allowance can only be changed in the code.
- Variants checked with one isolate each (#33).
- No real Alpha (B.1.1.7) genome: the file saved as Alpha was Beta (#42).
- Severity still does not change "binds"; only a HIGH 3' risk changes the
  product status, to "at risk" (#41). Mild mismatches near the 3' end and
  severe ones elsewhere are warnings only.
- Binding at the annealing temperature is an equilibrium estimate (#45), not
  a predicted Ct or efficiency.
- Tm is not computed for indel sites (no bulge-loop parameters) or two
  mismatches side by side (no tandem-mismatch parameters) (#44).
- Tm uses 250 nM for every oligo (the Primer3 reference); CDC uses 500 nM
  primers and 125 nM probe (about 1 C difference).
- Hairpins and primer dimers are not modelled.
- Primer design (#47-#57): ranks only by GC and product length, so many
  designs tie at a perfect 0 and the earliest wins; it does not reproduce
  the CDC designs; never proposes degenerate bases (CDC N3 probe has a Y);
  no hairpin or primer-dimer check; RNase P not checked (needs the human
  RPP30 segment from NCBI).

## Planned (not built yet)
- **Primer design, next:** compare with NCBI BLAST (our primers against
  NCBI's databases) and Primer-BLAST (its designs for the N gene); report
  each published design's rank; move the CDC check into `journal/`; add the
  second model's Tm to the website.
Factor 2's basic plan (nearest-neighbor Tm, Primer3's salt correction,
mismatch parameters, comparison with the 55 C annealing step, fraction
bound) was built in Log 9 (#43-#45). Still to do for it:
- **Indels and tandem mismatches:** add bulge-loop (SantaLucia & Hicks 2004)
  and tandem-mismatch parameters so those sites get a Tm too.
- **CDC concentrations:** 500 nM primers, 125 nM probe, and the master mix's
  real Mg2+ (check the CDC instructions for use).
- **Probe-specific check (TaqMan).** The probe is cleaved only if it is bound
  when the polymerase arrives. CDC probes sit about 10 C above the primers. A
  probe mismatch (e.g. Omicron N1, -6.8 C) can lower the signal and give a
  later Ct without losing the product. Turn fraction bound into an expected
  Ct shift.
- **Primer suggestions on the website:** rank pairs by `STATUS_ORDER_AT_TA`,
  with "at risk" and "weak binding" flagged and suggested last (#41, #45).
- **Later factors:** hairpins and primer-dimers (Primer3 ΔG columns already in
  `primer3_reference_tm.tsv`), the RT step (reverse primer starts the cDNA;
  RNA secondary structure at the site), and the CDC run conditions (primer and
  probe concentrations, cycling).
- **More genomes per variant**, and confirm C28311T across Omicron BA.1
  sequences (#33, #37).
- **Spike-targeting comparison assay.** Hypothesis: oligos in regions under
  strong selection (spike, especially the receptor-binding domain and
  antibody sites) gain mismatches across variants faster than oligos in
  functionally constrained regions (CDC N1/N2/N3 in the N gene). The change
  is driven by selection, not GC content, and shows up over months as
  variants replace one another, not within a few replication cycles.
  Plan: add a spike primer/probe set to `data/primers/` (the Yale spike
  Δ69-70 assay, Vogels et al. 2021, was found in Log 8), run it with both
  models alongside N1/N2/N3, and compare mismatch counts, `three_prime_risk`
  and dTm per variant. Its probe across Alpha's ΔH69/V70 doubles as the
  positive control (should fail on Alpha and Omicron BA.1; needs a real
  Alpha genome, #42).
