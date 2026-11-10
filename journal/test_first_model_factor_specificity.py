# =============================================================================
# TESTS FOR: FIRST MODEL: FACTOR: SPECIFICITY   (version 4)
# =============================================================================
# The model itself lives in first_model_factor_specificity.py; this file only
# tests it. Every case uses FAKE primers and a FAKE target built by hand, so
# each expected answer is worked out on paper first ("Cases solvable by hand"
# in the proposal). journal/test_cases_specificity.txt lists every case.
#
# Run it with:  python3 journal/test_first_model_factor_specificity.py
# =============================================================================

import random

from first_model_factor_specificity import (
    compare_site, print_report, reverse_complement, run_model, MIN_KNOWN_BASES,
)


# =============================================================================
# SECTION 1: FAKE PRIMERS AND FAKE TARGET BUILDING BLOCKS
# -----------------------------------------------------------------------------
# Every fake primer and probe is 18-30 bases long, the range of real qPCR
# oligos (the CDC primers and probes are 18-24 bases).
#
#   Set A (used by most cases): 20-base forward, 24-base reverse, 26-base probe
#   Set B (length-limit cases): 18-base forward, 30-base reverse
#   plus an 18-base palindrome (Case 14) and a 20-base overlapping reverse
#   primer (Case 12)
#
# The fake target is assembled from pieces of known length, so every expected
# position and product size below can be added up by hand.
#
# Filler is poly-A: no primer contains AAAAA or TTTTT, so the filler cannot
# match any 6-base seed and never creates accidental binding sites.
# =============================================================================

# -- Set A: 20 + 26 + 24 = 70 bp amplicon
FORWARD_A = "ACGGTCTAGCATTGCAGTCC"                  # 20 bases
REVERSE_A = "TCCGATAGTCAGCTTGGACAGTCA"              # 24 bases
REVERSE_SITE_A = reverse_complement(REVERSE_A)      # "TGACTGTCCAAGCTGACTATCGGA"
INSERT_A = "CGTACGTTGCATGCCTAGGTGATCAG"             # 26 bases between the sites
PROBE_A = INSERT_A                                  # the probe binds the whole insert

# -- Set B: 18 + 25 + 30 = 73 bp amplicon
FORWARD_B = "CAGTGCTTACGGATCGAG"                    # 18 bases
REVERSE_B = "GTCATCGACTGTAGCTTCGAGACGTTAGCC"        # 30 bases
REVERSE_SITE_B = reverse_complement(REVERSE_B)      # "GGCTAACGTCTCGAAGCTACAGTCGATGAC"
INSERT_B = "TGACCGTATCAGGTCATGCAGTCCA"              # 25 bases between the sites


def filler(n):
    """n bases of background sequence that no primer can bind."""
    return "A" * n


def amplicon():
    """Set A: forward site + insert + reverse site = 20 + 26 + 24 = 70 bp product.

    Amplicon offsets (0-based): forward site 0-19, insert 20-45, reverse site 46-69.
    """
    return FORWARD_A + INSERT_A + REVERSE_SITE_A


def amplicon_b():
    """Set B: forward site + insert + reverse site = 18 + 25 + 30 = 73 bp product.

    Amplicon offsets (0-based): forward site 0-17, insert 18-42, reverse site 43-72.
    """
    return FORWARD_B + INSERT_B + REVERSE_SITE_B


def mutate(seq, offset, new_base):
    """Return seq with the base at offset replaced (used to plant mismatches or Ns)."""
    return seq[:offset] + new_base + seq[offset + 1:]


def change_base(seq, offset):
    """Plant a guaranteed mismatch: A becomes C, anything else becomes A."""
    return mutate(seq, offset, "C" if seq[offset] == "A" else "A")


# =============================================================================
# SECTION 2: CASES SOLVABLE BY HAND: CORE (version 1)
# -----------------------------------------------------------------------------
# Each case states its expected answer (worked out on paper), runs the model,
# and checks the model agrees. A failed check stops the program.
# Positions in the checks are 0-based (as stored); the report prints 1-based.
# =============================================================================

def check(condition, message):
    if not condition:
        raise AssertionError(message)
    print(f"  PASS: {message}")


def product_lengths(result):
    """Lengths of the products whose status is 'product', smallest first."""
    return sorted(p["length"] for p in result["products"] if p["status"] == "product")


def site_of(result, primer):
    """The first site found for one primer."""
    return [s for s in result["sites"] if s["primer"] == primer][0]


# -- Case 1: exact match, one product ----------------------------------------
# a) Set A. target = 10 filler + 70 bp amplicon + 10 filler
#    forward binds "+" at 10..30, reverse binds "-" at 56..80, product 10..80 = 70 bp
# b) Set B, primers of 18 and 30 bases. Same layout:
#    forward "+" at 10..28, reverse "-" at 53..83, product 10..83 = 73 bp
def case_exact_match():
    target = filler(10) + amplicon() + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 1a: exact match, set A", result)
    check(product_lengths(result) == [70], "exactly one product, 70 bp (20 + 26 + 24)")
    product = result["products"][0]
    check((product["start"], product["end"]) == (10, 80), "product spans 10-80 (printed 11-80)")

    target = filler(10) + amplicon_b() + filler(10)
    result = run_model(FORWARD_B, REVERSE_B, target, max_mismatches=0)
    print_report("Case 1b: exact match, set B (18- and 30-base primers)", result)
    check(product_lengths(result) == [73], "exactly one product, 73 bp (18 + 25 + 30)")
    product = result["products"][0]
    check((product["start"], product["end"]) == (10, 83), "product spans 10-83 (printed 11-83)")


# -- Case 2: reverse-strand matches ------------------------------------------
# Reverse-complement the whole Case 1a target. The product is the same 70 bp,
# but now the REVERSE primer is the "+" site and the FORWARD primer is "-".
def case_reverse_strand():
    target = reverse_complement(filler(10) + amplicon() + filler(10))
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 2: whole target on the other strand", result)
    check(product_lengths(result) == [70], "exactly one product, still 70 bp")
    check(result["products"][0]["left"]["primer"] == "reverse", "reverse primer is now the '+' site")
    check(result["products"][0]["right"]["primer"] == "forward", "forward primer is now the '-' site")


# -- Case 3: introduced mismatches at the 3' end -----------------------------
# Forward site: last base (amplicon offset 19) C -> G, distance 1 from 3' end.
# Reverse site: first plus-strand base (amplicon offset 46) T -> C. The reverse
# primer's 3' end is the LEFT end of a "-" site, so this is also distance 1.
# In the reverse primer's own direction: its 3' base A now faces genome G
# (the complement of plus-strand C).
# Severity: the forward primer's C sits opposite template C (complement of
# genome G): C-C is SEVERE. The reverse primer's A sits opposite template C
# (complement of genome G): A-C is MILD.
def case_mismatches():
    seq = mutate(mutate(amplicon(), 19, "G"), 46, "C")
    target = filler(10) + seq + filler(10)

    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 3a: mismatches, tolerance 0", result)
    check(product_lengths(result) == [], "no product when no mismatches are allowed")
    check(result["outcomes"]["forward"] == "mismatched", "forward is 'mismatched'")
    check(result["assay"] == "no product", "assay is 'no product'")
    check(result["near_misses"]["forward"]["changes"] ==
          [{"from_3prime": 1, "primer": "C", "genome": "G", "severity": "severe"}],
          "forward near-miss: 1 from 3' end, primer C, genome G, severe")
    check(result["near_misses"]["reverse"]["changes"] ==
          [{"from_3prime": 1, "primer": "A", "genome": "G", "severity": "mild"}],
          "reverse near-miss: 1 from 3' end, primer A, genome G, mild")

    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 3b: mismatches, tolerance 1", result)
    check([(p["length"], p["status"]) for p in result["products"]] == [(70, "at risk")],
          "one 70 bp product when 1 mismatch is allowed, status 'at risk' (severe 3' mismatch)")
    check(result["assay"] == "at risk", "assay is 'at risk', not 'product'")
    check(result["products"][0]["left"]["mismatch_from_3prime"] == [1],
          "forward mismatch reported at the 3'-terminal base")
    check(result["products"][0]["right"]["mismatch_from_3prime"] == [1],
          "reverse mismatch reported at the 3'-terminal base")
    check(result["products"][0]["left"]["three_prime_risk"] == "high", "forward site: 3' risk HIGH (severe C-C)")
    check(result["products"][0]["right"]["three_prime_risk"] == "low", "reverse site: 3' risk LOW (mild A-C)")
    check(result["products"][0]["three_prime_risk"] == "high", "product takes the worse grade: HIGH")
    check(len(result["warnings"]) == 1, "one 3' risk warning is printed")


# -- Case 3d: an "at risk" product ranks below a clean one -------------------
# Case 4 layout (two amplicon copies, 30 filler between), but copy 1's forward
# site has the severe 3'-terminal mismatch of Case 3 (amplicon offset 19 C -> G).
# Tolerance 1. Copy 1 at 10..80, copy 2 at 110..180:
#   copy1 F x copy1 R:  70 bp, forward 3' risk HIGH -> at risk
#   copy1 F x copy2 R: 170 bp, same forward site    -> at risk
#   copy2 F x copy2 R:  70 bp, both sites perfect   -> product
# The assay takes the best status: product > at risk > ..., so "product".
def case_at_risk_ranking():
    target = filler(10) + mutate(amplicon(), 19, "G") + filler(30) + amplicon() + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 3d: one at-risk copy, one clean copy", result)
    check([(p["start"], p["length"], p["status"]) for p in result["products"]] ==
          [(10, 70, "at risk"), (10, 170, "at risk"), (110, 70, "product")],
          "copy 1 products are 'at risk', copy 2 product is clean")
    check(result["assay"] == "product", "assay is 'product': a clean product outranks 'at risk'")


# -- Case 3c: mismatch far from the 3' end -----------------------------------
# Change the forward site's FIRST base (offset 0) A -> G: distance 20 from the
# 3' end. Primer A opposite template C (complement of genome G): A-C, mild.
def case_mismatch_5prime():
    target = filler(10) + mutate(amplicon(), 0, "G") + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 3c: mismatch at the 5' end", result)
    check(result["products"][0]["left"]["changes"] ==
          [{"from_3prime": 20, "primer": "A", "genome": "G", "severity": "mild"}],
          "forward mismatch 20 from 3' end: primer A, genome G, mild (A-C)")
    check(result["products"][0]["three_prime_risk"] is None, "a 5' mismatch gives no 3' risk")


# -- Case 4: multiple products -----------------------------------------------
# Two copies of the amplicon: 10 filler + amplicon + 30 filler + amplicon + 10.
# Copy 1 at 10..80, copy 2 at 110..180. Products with max_length=1000:
#   copy1 F x copy1 R = 70 bp, copy2 F x copy2 R = 70 bp,
#   copy1 F x copy2 R = 180 - 10 = 170 bp  (the long "outer" product)
# copy2 F x copy1 R faces away -> no product.
def case_multiple_products():
    target = filler(10) + amplicon() + filler(30) + amplicon() + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 4: two copies of the target", result)
    check(len(result["sites"]) == 4, "4 binding sites (2 forward, 2 reverse)")
    check(product_lengths(result) == [70, 70, 170], "products of 70, 70 and 170 bp")

    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0, max_length=100)
    print_report("Case 4b: same, max product length 100", result)
    check(product_lengths(result) == [70, 70], "the 170 bp product is dropped when the limit is 100")


# -- Case 5: primers facing away from each other -----------------------------
# Put the reverse site FIRST and the forward site after it:
# they point away from each other, so nothing is amplified.
def case_facing_away():
    target = filler(10) + REVERSE_SITE_A + INSERT_A + FORWARD_A + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 5: primers facing away", result)
    check(len(result["sites"]) == 2, "both primers still bind")
    check(result["assay"] == "no product", "but no product forms")


# -- Case 6: no binding site at all ------------------------------------------
def case_no_site():
    result = run_model(FORWARD_A, REVERSE_A, filler(200), max_mismatches=2)
    print_report("Case 6: target without the primer sites", result)
    check(len(result["sites"]) == 0, "no sites")
    check(result["assay"] == "no product", "assay is 'no product'")
    check(result["near_misses"]["forward"] is None, "no near-miss either: nothing is close")


# =============================================================================
# SECTION 3: CASES SOLVABLE BY HAND: BAD SEQUENCING DATA (version 2)
# -----------------------------------------------------------------------------
# Each case imitates a trace that poor sequencing leaves in a finished genome
# and checks the model says "not assessable" instead of "mismatched".
# =============================================================================

# -- Case 7: Ns inside a primer site -----------------------------------------
# Replace amplicon offsets 5 and 6 (inside the forward site) with N.
# On the "+" strand, offset 5 is 20-5 = 15 and offset 6 is 20-6 = 14 bases
# from the 3' end. With tolerance 0 the Ns must NOT count as mismatches.
def case_n_in_site():
    seq = mutate(mutate(amplicon(), 5, "N"), 6, "N")
    target = filler(10) + seq + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 7: Ns inside the forward site", result)
    forward_site = site_of(result, "forward")
    check(forward_site["mismatches"] == 0, "Ns are not counted as mismatches")
    check(forward_site["unknown_from_3prime"] == [14, 15], "Ns reported 14 and 15 bases from the 3' end")
    check(result["outcomes"]["forward"] == "not assessable", "forward is 'not assessable'")
    check(result["outcomes"]["reverse"] == "binds", "reverse still binds")
    check(result["assay"] == "not assessable", "assay is 'not assessable'")


# -- Case 7b: same Ns, but 2 mismatches allowed (version 4) ------------------
# 0 real mismatches + 2 unknown bases = 2, which still fits an allowance of 2:
# even if both Ns were mismatches the primer would bind, so the missing data
# cannot change the answer. The site is "binds", not "not assessable".
def case_n_cannot_change_answer():
    seq = mutate(mutate(amplicon(), 5, "N"), 6, "N")
    target = filler(10) + seq + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=2)
    print_report("Case 7b: Ns inside the forward site, tolerance 2", result)
    check(site_of(result, "forward")["outcome"] == "binds", "forward binds: the Ns cannot change the answer")
    check(result["assay"] == "product", "assay is 'product'")


# -- Case 8: sequence ends partway through a primer site ---------------------
# target = 10 filler + forward site + insert + first 12 bases of the reverse
# site. The reverse site would span 56..80 but the sequence stops at 68, so
# 12 of its 24 bases are present (>= 10 required) and all of them match.
def case_truncated_end():
    target = filler(10) + FORWARD_A + INSERT_A + REVERSE_SITE_A[:12]
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 8: sequence stops inside the reverse site", result)
    check(site_of(result, "reverse")["start"] == 56, "reverse site would span 56-80")
    check(site_of(result, "reverse")["reason"] == "runs off the end of the sequence",
          "reason is 'runs off the end'")
    check(result["outcomes"]["reverse"] == "not assessable", "reverse is 'not assessable'")
    check(result["assay"] == "not assessable", "assay is 'not assessable'")


# -- Case 9: sequence starts partway through a primer site -------------------
# Drop the first 4 bases of the amplicon: 16 of the forward site's 20 bases
# remain at the very start of the sequence, so it "runs off the start".
def case_truncated_start():
    target = amplicon()[4:] + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 9: sequence starts inside the forward site", result)
    forward_site = site_of(result, "forward")
    check(forward_site["start"] == -4, "forward site would have started 4 bases before the sequence")
    check(forward_site["reason"] == "runs off the start of the sequence", "reason is 'runs off the start'")
    check(result["assay"] == "not assessable", "assay is 'not assessable'")


# -- Case 10: real mismatches next to missing data ---------------------------
# Four real mismatches in the reverse site (amplicon offsets 50, 55, 60, 65 =
# site offsets 4, 9, 14, 19) plus an N in the forward site (offset 3).
# Tolerance 1: the reverse primer has more real differences than allowed, so
# it is "mismatched". Because the reverse primer cannot bind anywhere, no
# product is possible whatever the N hides: the assay is "no product", NOT
# "not assessable" (fixed in version 3).
# The 4 mismatches, 5 bases apart, cut the 24-base reverse site into five
# stretches of 4 bases, so the near-miss search (which needs an exact 6-base
# stretch, k = 6) reports "none close". This shows that limit honestly.
# Version 4: the forward site has 0 mismatches + 1 N, which fits an allowance
# of 1 even if the N is a mismatch, so the forward primer "binds".
def case_mismatched_beats_unknown():
    seq = amplicon()
    for offset in (50, 55, 60, 65):
        seq = change_base(seq, offset)
    seq = mutate(seq, 3, "N")
    target = filler(10) + seq + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 10: reverse truly mismatched, N elsewhere", result)
    check(result["outcomes"]["reverse"] == "mismatched", "reverse is 'mismatched'")
    check(result["outcomes"]["forward"] == "binds", "forward binds: 1 N fits an allowance of 1")
    check(result["assay"] == "no product", "assay is 'no product': the reverse primer cannot bind")
    check(result["near_misses"]["reverse"] is None,
          "near-miss not found: no exact 6-base stretch left (known limit)")


# -- Case 10b: a near-miss that k = 6 misses but k = 5 would find -----------
# Three real mismatches in the 20-base forward site at offsets 5, 11 and 17
# cut it into exact stretches of 5, 5, 5 and 2 bases (offsets 0-4, 6-10,
# 12-16, 18-19). Tolerance 1, so the forward primer is "mismatched".
#   k = 6: no stretch is 6 long -> near-miss "none close"
#   k = 5: each 5-base stretch is a seed -> near-miss found, 3 events,
#          20 - 5 = 15, 20 - 11 = 9 and 20 - 17 = 3 bases from the 3' end
# This is the cost of the faster k = 6 index.
def case_near_miss_needs_k():
    seq = amplicon()
    for offset in (5, 11, 17):
        seq = change_base(seq, offset)
    target = filler(10) + seq + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 10b: 3 mismatches 6 apart, k = 6", result)
    check(result["outcomes"]["forward"] == "mismatched", "forward is 'mismatched'")
    check(result["near_misses"]["forward"] is None, "k = 6: no exact 6-base stretch, near-miss not found")

    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1, k=5)
    print_report("Case 10b: same target, k = 5", result)
    miss = result["near_misses"]["forward"]
    check(miss is not None and miss["events"] == 3, "k = 5: near-miss found, 3 events")
    check([c["from_3prime"] for c in miss["changes"]] == [3, 9, 15], "mismatches 3, 9 and 15 from the 3' end")


# =============================================================================
# SECTION 4: CASES SOLVABLE BY HAND: EDGE CASES (version 3)
# -----------------------------------------------------------------------------
# Situations the earlier cases did not cover.
# =============================================================================

# -- Case 11: a primer pairs with itself -------------------------------------
# Put the forward primer's reverse complement where the reverse site would be.
# The forward primer now binds "+" at 10..30 AND "-" at 56..76, so it makes a
# 20 + 26 + 20 = 66 bp product on its own; the reverse primer binds nowhere.
def case_self_pairing():
    target = filler(10) + FORWARD_A + INSERT_A + reverse_complement(FORWARD_A) + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 11: forward primer pairs with itself", result)
    product = result["products"][0]
    check((product["left"]["primer"], product["right"]["primer"]) == ("forward", "forward"),
          "the product is forward x forward")
    check(product_lengths(result) == [66], "it is 66 bp")
    check(result["outcomes"]["reverse"] == "mismatched", "the reverse primer does not bind")


# -- Case 12: overlapping primer sites ---------------------------------------
# A 20-base reverse primer whose site starts INSIDE the forward site:
# site = last 5 bases of the forward site + first 15 bases of the insert.
# Forward site 10..30, overlapping reverse site 25..45 -> product 10..45 = 35 bp.
OVERLAP_REVERSE = reverse_complement(FORWARD_A[-5:] + INSERT_A[:15])


def case_overlapping_sites():
    target = filler(10) + FORWARD_A + INSERT_A[:15] + filler(10)
    result = run_model(FORWARD_A, OVERLAP_REVERSE, target, max_mismatches=0)
    print_report("Case 12: primer sites overlap by 5 bases", result)
    check(site_of(result, "reverse")["start"] == 25, "reverse site starts at 25, inside the forward site")
    check(product_lengths(result) == [35], "one 35 bp product")


# -- Case 13: mismatch allowance too high for the primer length --------------
# 3 mismatches -> 4 pieces. The 30-base reverse primer gives 30 // 4 = 7-base
# pieces (fine), but the 18-base forward primer gives 18 // 4 = 4-base pieces,
# shorter than k = 6, so the search cannot guarantee it finds every site.
# It must stop with an error that names the forward primer.
def case_allowance_too_high():
    print("\n--- Case 13: 3 mismatches allowed on an 18-base primer ---")
    try:
        run_model(FORWARD_B, REVERSE_B, amplicon_b(), max_mismatches=3)
        check(False, "an error is raised")
    except ValueError as error:
        print(f"  error message: {error}")
        check(str(error).startswith("forward: 18 bases"),
              "an error names the 18-base forward primer instead of a silently incomplete search")


# -- Case 14: palindromic primer binds both strands at the same spot ---------
# GACTAGCTTAAGCTAGTC (18 bases) is its own reverse complement, so ONE place in
# the target is a "+" site and a "-" site at once (both 10..28). The two sites
# face each other, but the "product" would be only 18 bp (< 20 minimum).
PALINDROME = "GACTAGCTTAAGCTAGTC"


def case_palindrome():
    target = filler(10) + PALINDROME + filler(10)
    result = run_model(PALINDROME, REVERSE_A, target, max_mismatches=0)
    print_report("Case 14: palindromic primer", result)
    strands = sorted((s["strand"], s["start"]) for s in result["sites"] if s["primer"] == "forward")
    check(strands == [("+", 10), ("-", 10)], "same spot found on both strands")
    check(result["assay"] == "no product", "no product (18 bp is below the 20 bp minimum)")


# -- Case 15: degenerate base in the primer ----------------------------------
# Forward primer with Y (= C or T) at offset 9: ACGGTCTAGYATTGCAGTCC.
# Genome C or T at that spot -> 0 mismatches. Genome A -> 1 mismatch,
# 20 - 9 = 11 bases from the 3' end, reported as "primer Y, genome A".
# Severity: template T faces primer C or T (C-T, T-T): neither is severe -> mild.
DEGENERATE_FORWARD = mutate(FORWARD_A, 9, "Y")


def case_degenerate_primer():
    for genome_base in ("C", "T"):
        target = filler(10) + mutate(amplicon(), 9, genome_base) + filler(10)
        result = run_model(DEGENERATE_FORWARD, REVERSE_A, target, max_mismatches=0)
        print_report(f"Case 15: primer Y, genome {genome_base}", result)
        check(site_of(result, "forward")["mismatches"] == 0, f"Y matches genome {genome_base}")

    target = filler(10) + mutate(amplicon(), 9, "A") + filler(10)
    result = run_model(DEGENERATE_FORWARD, REVERSE_A, target, max_mismatches=1)
    print_report("Case 15: primer Y, genome A", result)
    check(site_of(result, "forward")["changes"] ==
          [{"from_3prime": 11, "primer": "Y", "genome": "A", "severity": "mild"}],
          "Y does not match genome A: 11 from 3' end, primer Y, genome A, mild")


# -- Case 16: probe binds inside the product ---------------------------------
# Probe = the 26-base insert, which sits between the primer sites.
def case_probe_binds():
    target = filler(10) + amplicon() + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, probe=PROBE_A, max_mismatches=0)
    print_report("Case 16: probe binds inside the product", result)
    check(result["products"][0]["probe"] == "binds", "probe binds inside the 70 bp product")
    check(result["assay"] == "product", "assay is 'product'")


# -- Case 17: primers bind but the probe does not ----------------------------
# Three real mismatches in the insert (amplicon offsets 24, 32, 40), tolerance 1.
# Both primers still bind and copy the DNA, but the probe cannot bind, so a
# qPCR machine would see no fluorescence: "no signal".
def case_probe_mismatched():
    seq = amplicon()
    for offset in (24, 32, 40):
        seq = change_base(seq, offset)
    target = filler(10) + seq + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, probe=PROBE_A, max_mismatches=1)
    print_report("Case 17: probe has 3 mismatches", result)
    check(result["outcomes"]["forward"] == "binds" and result["outcomes"]["reverse"] == "binds",
          "both primers bind")
    check(result["products"][0]["probe"] == "mismatched", "the probe does not bind")
    check(result["assay"] == "no signal", "assay is 'no signal'")


# -- Case 18: deletion inside a primer site ----------------------------------
# Delete 2 bases ("TA", amplicon offsets 6-7) from the 20-base forward site.
# Without gaps this looks like many mismatches; the near-miss search should
# explain it as ONE event: a 2-base deletion of primer bases 14 and 13 from
# the 3' end. (Gaps longer than 1 base are never allowed to bind.)
def case_deletion():
    seq = amplicon()
    target = filler(10) + seq[:6] + seq[8:] + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 18: 2-base deletion in the forward site", result)
    miss = result["near_misses"]["forward"]
    check(result["outcomes"]["forward"] == "mismatched", "forward is 'mismatched'")
    check(miss["events"] == 1, "near-miss is 1 event, not many mismatches")
    check(miss["indel"] == {"kind": "deletion", "length": 2, "from_3prime": (13, 14)},
          "explained as a 2-base deletion, 13-14 from the 3' end")


# -- Case 18b: a 6-base deletion is too long for an 18-base primer ----------
# Set B. Delete 6 bases (offsets 6-11) from the 18-base forward site: only
# 18 - 6 = 12 primer bases would still match, fewer than MIN_MATCHED_BASES
# (14), so the model refuses to guess.
def case_deletion_too_long():
    seq = amplicon_b()
    target = filler(10) + seq[:6] + seq[12:] + filler(10)
    result = run_model(FORWARD_B, REVERSE_B, target, max_mismatches=1)
    print_report("Case 18b: 6-base deletion in an 18-base primer site", result)
    check(result["outcomes"]["forward"] == "mismatched", "forward is 'mismatched'")
    check(result["near_misses"]["forward"] is None, "no near-miss: too few matching bases to be sure")


# -- Case 18c: the same 6-base deletion IS found in a 30-base primer --------
# Set B. Delete 6 bases from the middle of the 30-base reverse site (amplicon
# offsets 55-60 = site offsets 12-17): 30 - 6 = 24 bases still match (>= 14),
# so it is explained as ONE 6-base deletion, like Alpha's spike 69-70 deletion.
def case_deletion_long_primer():
    seq = amplicon_b()
    target = filler(10) + seq[:55] + seq[61:] + filler(10)
    result = run_model(FORWARD_B, REVERSE_B, target, max_mismatches=1)
    print_report("Case 18c: 6-base deletion in a 30-base primer site", result)
    miss = result["near_misses"]["reverse"]
    check(result["outcomes"]["reverse"] == "mismatched", "reverse is 'mismatched'")
    check(miss["events"] == 1, "near-miss is 1 event")
    check((miss["indel"]["kind"], miss["indel"]["length"]) == ("deletion", 6), "explained as a 6-base deletion")


# -- Case 19: insertion inside a primer site ---------------------------------
# Insert TTT into the middle of the 24-base reverse site (after amplicon offset 57).
def case_insertion():
    seq = amplicon()
    target = filler(10) + seq[:58] + "TTT" + seq[58:] + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 19: 3-base insertion in the reverse site", result)
    miss = result["near_misses"]["reverse"]
    check(result["outcomes"]["reverse"] == "mismatched", "reverse is 'mismatched'")
    check(miss["events"] == 1, "near-miss is 1 event")
    check((miss["indel"]["kind"], miss["indel"]["length"]) == ("insertion", 3), "explained as a 3-base insertion")


# =============================================================================
# SECTION 5: CASES SOLVABLE BY HAND: VERSION 4 LOGIC FIXES
# -----------------------------------------------------------------------------
# Genome ambiguity codes, graded 3' risk, and binding through a small
# insertion/deletion in the primer's 5' half. (Case 7b above covers the
# "missing data cannot change the answer" fix.)
# =============================================================================

# -- Case 20: ambiguity codes in the genome ----------------------------------
# The forward primer has C at offset 9 (11 bases from the 3' end).
#   genome Y (C or T): might match, might not      -> unknown
#   genome R (A or G): cannot match a C             -> certain mismatch
#                      (template Y = C or T; C-C is severe -> severe)
#   primer Y vs genome Y: every option is accepted  -> certain match
def case_genome_ambiguity():
    target = filler(10) + mutate(amplicon(), 9, "Y") + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 20a: primer C, genome Y", result)
    check(site_of(result, "forward")["unknown_from_3prime"] == [11], "genome Y is unknown, not a mismatch")
    check(result["outcomes"]["forward"] == "not assessable", "forward is 'not assessable'")

    target = filler(10) + mutate(amplicon(), 9, "R") + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 20b: primer C, genome R, tolerance 0", result)
    check(result["outcomes"]["forward"] == "mismatched", "genome R is a certain mismatch: 'mismatched'")
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 20b: primer C, genome R, tolerance 1", result)
    check(site_of(result, "forward")["changes"] ==
          [{"from_3prime": 11, "primer": "C", "genome": "R", "severity": "severe"}],
          "reported as primer C, genome R, severe")

    target = filler(10) + mutate(amplicon(), 9, "Y") + filler(10)
    result = run_model(DEGENERATE_FORWARD, REVERSE_A, target, max_mismatches=0)
    print_report("Case 20c: primer Y, genome Y", result)
    forward_site = site_of(result, "forward")
    check(forward_site["outcome"] == "binds" and forward_site["unknown_from_3prime"] == [],
          "primer Y fully accepts genome Y: a certain match")


# -- Case 21: graded 3' risk -------------------------------------------------
# a) Only the reverse site's 3'-terminal mismatch (amplicon offset 46 T -> C):
#    A-C is mild, so the product's 3' risk is LOW.
# b) The forward site's 3'-terminal base is N, tolerance 1: the site binds
#    even if the N is a mismatch, but the 3' end is uncertain: POSSIBLE.
def case_graded_risk():
    target = filler(10) + mutate(amplicon(), 46, "C") + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 21a: one mild mismatch at the 3' end", result)
    check(result["products"][0]["three_prime_risk"] == "low", "product 3' risk is LOW")
    check(result["assay"] == "product", "a LOW 3' risk does not make it 'at risk'")

    target = filler(10) + mutate(amplicon(), 19, "N") + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 21b: N at the forward primer's 3' end", result)
    check(site_of(result, "forward")["outcome"] == "binds", "forward binds even if the N is a mismatch")
    check(result["products"][0]["three_prime_risk"] == "possible", "product 3' risk is POSSIBLE")
    check(result["assay"] == "product", "a POSSIBLE 3' risk does not make it 'at risk'")


# -- Case 22: 1-base deletion in the primer's 5' half ------------------------
# Delete amplicon offset 2 (a G, 18 bases from the 3' end; the G beside it
# makes 17 an equally good answer). Both are more than 20 // 2 = 10 from the
# 3' end, so the 3' half is intact and with 1 mismatch allowed the forward
# primer still binds, counting the gap as 1 mismatch.
# The product loses 1 base: 69 bp. With 0 mismatches allowed it does not bind.
def case_small_deletion_5prime():
    seq = amplicon()
    target = filler(10) + seq[:2] + seq[3:] + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 22: 1-base deletion in the 5' half, tolerance 1", result)
    forward_site = site_of(result, "forward")
    check((forward_site["indel"]["kind"], forward_site["indel"]["length"]) == ("deletion", 1),
          "forward binds through a 1-base deletion")
    check(forward_site["indel"]["from_3prime"][0] in (17, 18), "the deletion is 17 or 18 from the 3' end (5' half)")
    check(product_lengths(result) == [69], "product is 69 bp (one base shorter)")

    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=0)
    print_report("Case 22: same, tolerance 0", result)
    check(result["outcomes"]["forward"] == "mismatched", "with 0 allowed the gap is too much")


# -- Case 22b: 1-base deletion in the primer's 3' half -----------------------
# Delete amplicon offset 16 (a G, 4 bases from the 3' end; its neighbours T and
# A differ, so the spot is unambiguous). Read without a gap, the site's last
# 4 bases GTCC face TCCC: 2 mismatches, more than the allowance of 1. A gap
# that close to the 3' end is not allowed to bind; the near-miss search
# explains it as a 1-base deletion 4 from the 3' end.
def case_small_deletion_3prime():
    seq = amplicon()
    target = filler(10) + seq[:16] + seq[17:] + filler(10)
    result = run_model(FORWARD_A, REVERSE_A, target, max_mismatches=1)
    print_report("Case 22b: 1-base deletion in the 3' half", result)
    check(result["outcomes"]["forward"] == "mismatched", "forward is 'mismatched'")
    check(result["near_misses"]["forward"]["indel"] ==
          {"kind": "deletion", "length": 1, "from_3prime": (4, 4)},
          "near-miss: 1-base deletion 4 bases from the 3' end")


# =============================================================================
# SECTION 6: CHECKING THE SEARCH ITSELF
# -----------------------------------------------------------------------------
# The k-mer search is fast because it skips most positions. These checks make
# sure it never skips a real site, and never invents near-misses.
# =============================================================================

def brute_force_sites(primer, target, max_mismatches):
    """The slow way: compare the primer at EVERY position on both strands."""
    found = set()
    for strand, pattern in (("+", primer), ("-", reverse_complement(primer))):
        for start in range(-len(pattern) + 1, len(target)):
            mismatches, unknowns = compare_site(target, start, pattern)
            if len(mismatches) <= max_mismatches and len(pattern) - len(unknowns) >= MIN_KNOWN_BASES:
                found.add((strand, start))
    return found


def random_test_target(rng, primer):
    """A 300-base random target with mutated primer copies, Ns/Ys and maybe a cut end."""
    target = [rng.choice("ACGT") for _ in range(300)]
    for _ in range(4):
        copy = list(primer if rng.random() < 0.5 else reverse_complement(primer))
        for _ in range(rng.randint(0, 3)):
            copy[rng.randrange(len(copy))] = rng.choice("ACGT")
        position = rng.randrange(-6, 300 - len(copy) + 6)
        for i, base in enumerate(copy):
            if 0 <= position + i < len(target):
                target[position + i] = base
    for i in range(len(target)):
        if rng.random() < 0.02:
            target[i] = rng.choice("NY")
    return "".join(target)


# -- Case 23: fast search finds exactly what the slow search finds -----------
# Primers of 18, 20, 24 and 30 bases x 50 random targets each x mismatch
# allowances 0, 1 and 2. Plain (ungapped) sites from the model must equal the
# slow every-position search exactly.
def case_search_complete():
    print("\n--- Case 23: fast search vs slow every-position search ---")
    rng = random.Random(23)
    searches, compared = 0, 0
    for primer in (FORWARD_B, FORWARD_A, REVERSE_A, REVERSE_B):
        for _ in range(50):
            target = random_test_target(rng, primer)
            for max_mismatches in (0, 1, 2):
                result = run_model(primer, FORWARD_A, target, max_mismatches=max_mismatches)
                fast = {(s["strand"], s["start"]) for s in result["sites"]
                        if s["primer"] == "forward" and s["indel"] is None}
                slow = brute_force_sites(primer, target, max_mismatches)
                if fast != slow:
                    raise AssertionError(f"{primer}: fast {sorted(fast)} != slow {sorted(slow)} for {target}")
                searches += 1
                compared += len(slow)
    check(True, f"{searches} searches agree with the slow search ({compared} sites compared)")


# -- Case 24: no invented near-misses on a random genome ---------------------
# A random 30,000-base sequence (the size of SARS-CoV-2) should not contain
# the CDC N1 primers or probe, and the near-miss search must not invent a
# "close" site there. (Before the MIN_MATCHED_BASES rule, it did.)
def case_random_genome():
    print("\n--- Case 24: CDC N1 primers on a random 30,000-base sequence ---")
    rng = random.Random(24)
    genome = "".join(rng.choice("ACGT") for _ in range(30000))
    result = run_model("GACCCCAAAATCAGCGAAAT", "TCTGGTTACTGCCAGTTGAATCTG", genome,
                       probe="ACCCCGCATTACGTTTGGTGGACC")
    check(result["sites"] == [], "no binding sites")
    check(all(miss is None for miss in result["near_misses"].values()), "no near-misses invented")


# =============================================================================
# SECTION 7: RUN EVERY CASE
# =============================================================================

if __name__ == "__main__":
    case_exact_match()
    case_reverse_strand()
    case_mismatches()
    case_mismatch_5prime()
    case_at_risk_ranking()
    case_multiple_products()
    case_facing_away()
    case_no_site()
    case_n_in_site()
    case_n_cannot_change_answer()
    case_truncated_end()
    case_truncated_start()
    case_mismatched_beats_unknown()
    case_near_miss_needs_k()
    case_self_pairing()
    case_overlapping_sites()
    case_allowance_too_high()
    case_palindrome()
    case_degenerate_primer()
    case_probe_binds()
    case_probe_mismatched()
    case_deletion()
    case_deletion_too_long()
    case_deletion_long_primer()
    case_insertion()
    case_genome_ambiguity()
    case_graded_risk()
    case_small_deletion_5prime()
    case_small_deletion_3prime()
    case_search_complete()
    case_random_genome()
    print("\nAll cases passed.")
