# Modification Documentation: First Model, Factor: Specificity

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
**Status: OPEN.** Left as a FAIL; check the 73 bp value against Lu et al.
(2020) before changing it.

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

**Commits:** `bfaff5d`, `526c2a5` (the restore in #31), then #32 (ID filled in at the next update)

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

## Still open (known limits after version 4)
- N1 amplicon size: 72 bp computed vs 73 bp in the table (#23).
- Near-miss and small-gap searches need an exact 5-base stretch (#15).
- Near-miss limits (#16, #17) were tuned on random sequence, not real genomes
  with repeats.
- Severity is a two-level rule (#20); real effects also depend on
  neighbouring bases and positions beyond the 3' window.
- Sites running off the start print negative positions (e.g. -3-48).
- Input is not cleaned of pasted spaces or line breaks.
- The k-mer index is rebuilt on every run; the test suite takes ~14 s.
- The mismatch allowance can only be changed in the code.
- Not yet run on the Alpha, Delta and Omicron genomes.
