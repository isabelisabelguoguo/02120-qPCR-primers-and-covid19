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
# The fake target is assembled from pieces of known length, so every expected
# position and product size below can be added up by hand.
#
# Filler is poly-A: it cannot match any 5-base seed of these primers, so it
# never creates accidental binding sites.
# =============================================================================

FAKE_FORWARD = "GTTCAGAGCTACTGAC"   # 16 bases
FAKE_REVERSE = "TCCGATAGTCAGCTTG"   # 16 bases

# Where the reverse primer sits on the plus strand (as its reverse complement)
REVERSE_SITE = reverse_complement(FAKE_REVERSE)   # "CAAGCTGACTATCGGA"

INSERT = "CGTACGTTGCATGCCTAGGT"     # 20 bases between the two primer sites
FAKE_PROBE = INSERT                 # the probe binds the whole insert


def filler(n):
    """n bases of background sequence that no primer can bind."""
    return "A" * n


def amplicon():
    """forward site + insert + reverse site = 16 + 20 + 16 = 52 bp product."""
    return FAKE_FORWARD + INSERT + REVERSE_SITE


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
# target = 10 filler + 52 bp amplicon + 10 filler
# forward binds "+" at 10..26, reverse binds "-" at 46..62, product 10..62 = 52 bp
def case_exact_match():
    target = filler(10) + amplicon() + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 1: exact match", result)
    check(product_lengths(result) == [52], "exactly one product, 52 bp (16 + 20 + 16)")
    product = result["products"][0]
    check((product["start"], product["end"]) == (10, 62), "product spans 10-62 (printed 11-62)")


# -- Case 2: reverse-strand matches ------------------------------------------
# Reverse-complement the whole Case 1 target. The product is the same 52 bp,
# but now the REVERSE primer is the "+" site and the FORWARD primer is "-".
def case_reverse_strand():
    target = reverse_complement(filler(10) + amplicon() + filler(10))
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 2: whole target on the other strand", result)
    check(product_lengths(result) == [52], "exactly one product, still 52 bp")
    check(result["products"][0]["left"]["primer"] == "reverse", "reverse primer is now the '+' site")
    check(result["products"][0]["right"]["primer"] == "forward", "forward primer is now the '-' site")


# -- Case 3: introduced mismatches at the 3' end -----------------------------
# Forward site: last base (amplicon offset 15) C -> G, distance 1 from 3' end.
# Reverse site: first plus-strand base (amplicon offset 36) C -> T. The reverse
# primer's 3' end is the LEFT end of a "-" site, so this is also distance 1.
# In the reverse primer's own direction: its 3' base G now faces genome A.
# Severity (version 4): the forward primer's C sits opposite template C
# (complement of genome G): C-C is SEVERE. The reverse primer's G sits
# opposite template T (complement of genome A): G-T is MILD.
def case_mismatches():
    seq = mutate(mutate(amplicon(), 15, "G"), 36, "T")
    target = filler(10) + seq + filler(10)

    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 3a: mismatches, tolerance 0", result)
    check(product_lengths(result) == [], "no product when no mismatches are allowed")
    check(result["outcomes"]["forward"] == "mismatched", "forward is 'mismatched'")
    check(result["assay"] == "no product", "assay is 'no product'")
    check(result["near_misses"]["forward"]["changes"] ==
          [{"from_3prime": 1, "primer": "C", "genome": "G", "severity": "severe"}],
          "forward near-miss: 1 from 3' end, primer C, genome G, severe")
    check(result["near_misses"]["reverse"]["changes"] ==
          [{"from_3prime": 1, "primer": "G", "genome": "A", "severity": "mild"}],
          "reverse near-miss: 1 from 3' end, primer G, genome A, mild")

    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 3b: mismatches, tolerance 1", result)
    check(product_lengths(result) == [52], "one product when 1 mismatch is allowed")
    check(result["products"][0]["left"]["mismatch_from_3prime"] == [1],
          "forward mismatch reported at the 3'-terminal base")
    check(result["products"][0]["right"]["mismatch_from_3prime"] == [1],
          "reverse mismatch reported at the 3'-terminal base")
    check(result["products"][0]["left"]["three_prime_risk"] == "high", "forward site: 3' risk HIGH (severe C-C)")
    check(result["products"][0]["right"]["three_prime_risk"] == "low", "reverse site: 3' risk LOW (mild G-T)")
    check(result["products"][0]["three_prime_risk"] == "high", "product takes the worse grade: HIGH")
    check(len(result["warnings"]) == 1, "one 3' risk warning is printed")


# -- Case 3c: mismatch far from the 3' end -----------------------------------
# Change the forward site's FIRST base (offset 0): distance 16 from the 3' end.
def case_mismatch_5prime():
    target = filler(10) + mutate(amplicon(), 0, "C") + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 3c: mismatch at the 5' end", result)
    check(result["products"][0]["left"]["changes"] ==
          [{"from_3prime": 16, "primer": "G", "genome": "C", "severity": "mild"}],
          "forward mismatch 16 from 3' end: primer G, genome C, mild (G-G)")
    check(result["products"][0]["three_prime_risk"] is None, "a 5' mismatch gives no 3' risk")


# -- Case 4: multiple products -----------------------------------------------
# Two copies of the amplicon: 10 filler + amplicon + 30 filler + amplicon + 10.
# Copy 1 at 10..62, copy 2 at 92..144. Products with max_length=1000:
#   copy1 F x copy1 R = 52 bp, copy2 F x copy2 R = 52 bp,
#   copy1 F x copy2 R = 144 - 10 = 134 bp  (the long "outer" product)
# copy2 F x copy1 R faces away -> no product.
def case_multiple_products():
    target = filler(10) + amplicon() + filler(30) + amplicon() + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 4: two copies of the target", result)
    check(len(result["sites"]) == 4, "4 binding sites (2 forward, 2 reverse)")
    check(product_lengths(result) == [52, 52, 134], "products of 52, 52 and 134 bp")

    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0, max_length=100)
    print_report("Case 4b: same, max product length 100", result)
    check(product_lengths(result) == [52, 52], "the 134 bp product is dropped when the limit is 100")


# -- Case 5: primers facing away from each other -----------------------------
# Put the reverse site FIRST and the forward site after it:
# they point away from each other, so nothing is amplified.
def case_facing_away():
    target = filler(10) + REVERSE_SITE + INSERT + FAKE_FORWARD + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 5: primers facing away", result)
    check(len(result["sites"]) == 2, "both primers still bind")
    check(result["assay"] == "no product", "but no product forms")


# -- Case 6: no binding site at all ------------------------------------------
def case_no_site():
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, filler(200), max_mismatches=2)
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
# On the "+" strand, offset 5 is 16-5 = 11 and offset 6 is 16-6 = 10 bases
# from the 3' end. With tolerance 0 the Ns must NOT count as mismatches.
def case_n_in_site():
    seq = mutate(mutate(amplicon(), 5, "N"), 6, "N")
    target = filler(10) + seq + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 7: Ns inside the forward site", result)
    forward_site = site_of(result, "forward")
    check(forward_site["mismatches"] == 0, "Ns are not counted as mismatches")
    check(forward_site["unknown_from_3prime"] == [10, 11], "Ns reported 10 and 11 bases from the 3' end")
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
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=2)
    print_report("Case 7b: Ns inside the forward site, tolerance 2", result)
    check(site_of(result, "forward")["outcome"] == "binds", "forward binds: the Ns cannot change the answer")
    check(result["assay"] == "product", "assay is 'product'")


# -- Case 8: sequence ends partway through a primer site ---------------------
# target = 10 filler + forward site + insert + first 12 bases of the reverse
# site. The reverse site would span 46..62 but the sequence stops at 58, so
# 12 of its 16 bases are present (>= 10 required) and all of them match.
def case_truncated_end():
    target = filler(10) + FAKE_FORWARD + INSERT + REVERSE_SITE[:12]
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 8: sequence stops inside the reverse site", result)
    check(site_of(result, "reverse")["reason"] == "runs off the end of the sequence",
          "reason is 'runs off the end'")
    check(result["outcomes"]["reverse"] == "not assessable", "reverse is 'not assessable'")
    check(result["assay"] == "not assessable", "assay is 'not assessable'")


# -- Case 9: sequence starts partway through a primer site -------------------
# Drop the first 4 bases of the amplicon: 12 of the forward site's 16 bases
# remain at the very start of the sequence, so it "runs off the start".
def case_truncated_start():
    target = amplicon()[4:] + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 9: sequence starts inside the forward site", result)
    forward_site = site_of(result, "forward")
    check(forward_site["start"] == -4, "forward site would have started 4 bases before the sequence")
    check(forward_site["reason"] == "runs off the start of the sequence", "reason is 'runs off the start'")
    check(result["assay"] == "not assessable", "assay is 'not assessable'")


# -- Case 10: real mismatches next to missing data ---------------------------
# Three real mismatches in the reverse site (amplicon offsets 40, 44, 48)
# plus an N in the forward site. Tolerance 1: the reverse primer has more
# real differences than allowed, so it is "mismatched". Because the reverse
# primer cannot bind anywhere, no product is possible whatever the N hides:
# the assay is "no product", NOT "not assessable" (fixed in version 3).
# The 3 mismatches are spaced 4 bases apart, cutting the reverse site into
# stretches of at most 4 bases, so the near-miss search (which needs an exact
# 5-base stretch) reports "none close". This shows that limit honestly.
# Version 4: the forward site has 0 mismatches + 1 N, which fits an allowance
# of 1 even if the N is a mismatch, so the forward primer now "binds".
def case_mismatched_beats_unknown():
    seq = amplicon()
    for offset in (40, 44, 48):
        seq = change_base(seq, offset)
    seq = mutate(seq, 3, "N")
    target = filler(10) + seq + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 10: reverse truly mismatched, N elsewhere", result)
    check(result["outcomes"]["reverse"] == "mismatched", "reverse is 'mismatched'")
    check(result["outcomes"]["forward"] == "binds", "forward binds: 1 N fits an allowance of 1")
    check(result["assay"] == "no product", "assay is 'no product': the reverse primer cannot bind")
    check(result["near_misses"]["reverse"] is None,
          "near-miss not found: no exact 5-base stretch left (known limit)")


# =============================================================================
# SECTION 4: CASES SOLVABLE BY HAND: EDGE CASES (version 3)
# -----------------------------------------------------------------------------
# Situations the earlier cases did not cover.
# =============================================================================

# -- Case 11: a primer pairs with itself -------------------------------------
# Put the forward primer's reverse complement where the reverse site would be.
# The forward primer now binds "+" at 10..26 AND "-" at 46..62, so it makes a
# 52 bp product on its own; the reverse primer binds nowhere.
def case_self_pairing():
    target = filler(10) + FAKE_FORWARD + INSERT + reverse_complement(FAKE_FORWARD) + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 11: forward primer pairs with itself", result)
    product = result["products"][0]
    check((product["left"]["primer"], product["right"]["primer"]) == ("forward", "forward"),
          "the product is forward x forward")
    check(product_lengths(result) == [52], "it is 52 bp")
    check(result["outcomes"]["reverse"] == "mismatched", "the reverse primer does not bind")


# -- Case 12: overlapping primer sites ---------------------------------------
# A second reverse primer whose site starts INSIDE the forward site:
# site = last 5 bases of the forward site + first 11 bases of the insert.
# Forward site 10..26, overlapping reverse site 21..37 -> product 10..37 = 27 bp.
OVERLAP_REVERSE = reverse_complement(FAKE_FORWARD[-5:] + INSERT[:11])


def case_overlapping_sites():
    target = filler(10) + FAKE_FORWARD + INSERT[:11] + filler(10)
    result = run_model(FAKE_FORWARD, OVERLAP_REVERSE, target, max_mismatches=0)
    print_report("Case 12: primer sites overlap by 5 bases", result)
    check(site_of(result, "reverse")["start"] == 21, "reverse site starts at 21, inside the forward site")
    check(product_lengths(result) == [27], "one 27 bp product")


# -- Case 13: mismatch allowance too high for the primer length --------------
# 16 bases // (3 + 1) pieces = 4 bases per piece, shorter than k = 5,
# so the search cannot guarantee it finds every site. It must stop with an error.
def case_allowance_too_high():
    print("\n--- Case 13: 3 mismatches allowed on a 16-base primer ---")
    try:
        run_model(FAKE_FORWARD, FAKE_REVERSE, amplicon(), max_mismatches=3)
        check(False, "an error is raised")
    except ValueError as error:
        print(f"  error message: {error}")
        check(True, "an error is raised instead of a silently incomplete search")


# -- Case 14: palindromic primer binds both strands at the same spot ---------
# GCATGCATATGCATGC is its own reverse complement, so ONE place in the target
# is a "+" site and a "-" site at once (both 10..26). The two sites face each
# other, but the "product" would be only 16 bp (< 20 minimum): no product.
PALINDROME = "GCATGCATATGCATGC"


def case_palindrome():
    target = filler(10) + PALINDROME + filler(10)
    result = run_model(PALINDROME, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 14: palindromic primer", result)
    strands = sorted((s["strand"], s["start"]) for s in result["sites"] if s["primer"] == "forward")
    check(strands == [("+", 10), ("-", 10)], "same spot found on both strands")
    check(result["assay"] == "no product", "no product (16 bp is below the 20 bp minimum)")


# -- Case 15: degenerate base in the primer ----------------------------------
# Forward primer with Y (= C or T) at offset 8: GTTCAGAGYTACTGAC.
# Genome C or T at that spot -> 0 mismatches. Genome A -> 1 mismatch,
# 16 - 8 = 8 bases from the 3' end, reported as "primer Y, genome A".
DEGENERATE_FORWARD = mutate(FAKE_FORWARD, 8, "Y")


def case_degenerate_primer():
    for genome_base in ("C", "T"):
        target = filler(10) + mutate(amplicon(), 8, genome_base) + filler(10)
        result = run_model(DEGENERATE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
        print_report(f"Case 15: primer Y, genome {genome_base}", result)
        check(site_of(result, "forward")["mismatches"] == 0, f"Y matches genome {genome_base}")

    target = filler(10) + mutate(amplicon(), 8, "A") + filler(10)
    result = run_model(DEGENERATE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 15: primer Y, genome A", result)
    check(site_of(result, "forward")["changes"] ==
          [{"from_3prime": 8, "primer": "Y", "genome": "A", "severity": "mild"}],
          "Y does not match genome A: 8 from 3' end, primer Y, genome A, mild")


# -- Case 16: probe binds inside the product ---------------------------------
# Probe = the 20-base insert, which sits between the primer sites.
def case_probe_binds():
    target = filler(10) + amplicon() + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, probe=FAKE_PROBE, max_mismatches=0)
    print_report("Case 16: probe binds inside the product", result)
    check(result["products"][0]["probe"] == "binds", "probe binds inside the 52 bp product")
    check(result["assay"] == "product", "assay is 'product'")


# -- Case 17: primers bind but the probe does not ----------------------------
# Three real mismatches in the insert (amplicon offsets 18, 25, 31), tolerance 1.
# Both primers still bind and copy the DNA, but the probe cannot bind, so a
# qPCR machine would see no fluorescence: "no signal".
def case_probe_mismatched():
    seq = amplicon()
    for offset in (18, 25, 31):
        seq = change_base(seq, offset)
    target = filler(10) + seq + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, probe=FAKE_PROBE, max_mismatches=1)
    print_report("Case 17: probe has 3 mismatches", result)
    check(result["outcomes"]["forward"] == "binds" and result["outcomes"]["reverse"] == "binds",
          "both primers bind")
    check(result["products"][0]["probe"] == "mismatched", "the probe does not bind")
    check(result["assay"] == "no signal", "assay is 'no signal'")


# -- Case 18: deletion inside a primer site ----------------------------------
# Delete 2 bases from the middle of the forward site (amplicon offsets 6-7).
# Without gaps this looks like many mismatches; the near-miss search should
# explain it as ONE event: a 2-base deletion. (2 is the longest deletion a
# 16-base primer can show: 16 - 2 = 14 = MIN_MATCHED_BASES.)
def case_deletion():
    seq = amplicon()
    target = filler(10) + seq[:6] + seq[8:] + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 18: 2-base deletion in the forward site", result)
    miss = result["near_misses"]["forward"]
    check(result["outcomes"]["forward"] == "mismatched", "forward is 'mismatched'")
    check(miss["events"] == 1, "near-miss is 1 event, not many mismatches")
    check((miss["indel"]["kind"], miss["indel"]["length"]) == ("deletion", 2), "explained as a 2-base deletion")


# -- Case 18b: a 6-base deletion is too long for a 16-base primer -----------
# Delete 6 bases (offsets 5-10): only 10 primer bases would still match,
# fewer than MIN_MATCHED_BASES (14), so the model refuses to guess.
def case_deletion_too_long():
    seq = amplicon()
    target = filler(10) + seq[:5] + seq[11:] + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 18b: 6-base deletion in a 16-base primer site", result)
    check(result["near_misses"]["forward"] is None, "no near-miss: too few matching bases to be sure")


# -- Case 19: insertion inside a primer site ---------------------------------
# Insert TTT into the middle of the reverse site (after amplicon offset 43).
def case_insertion():
    seq = amplicon()
    target = filler(10) + seq[:44] + "TTT" + seq[44:] + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
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
# The forward primer has C at offset 8 (8 bases from the 3' end).
#   genome Y (C or T): might match, might not      -> unknown
#   genome R (A or G): cannot match a C             -> certain mismatch
#   primer Y vs genome Y: every option is accepted  -> certain match
def case_genome_ambiguity():
    target = filler(10) + mutate(amplicon(), 8, "Y") + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 20a: primer C, genome Y", result)
    check(site_of(result, "forward")["unknown_from_3prime"] == [8], "genome Y is unknown, not a mismatch")
    check(result["outcomes"]["forward"] == "not assessable", "forward is 'not assessable'")

    target = filler(10) + mutate(amplicon(), 8, "R") + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 20b: primer C, genome R, tolerance 0", result)
    check(result["outcomes"]["forward"] == "mismatched", "genome R is a certain mismatch: 'mismatched'")
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 20b: primer C, genome R, tolerance 1", result)
    check(site_of(result, "forward")["changes"] ==
          [{"from_3prime": 8, "primer": "C", "genome": "R", "severity": "severe"}],
          "reported as primer C, genome R, severe")

    target = filler(10) + mutate(amplicon(), 8, "Y") + filler(10)
    result = run_model(DEGENERATE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 20c: primer Y, genome Y", result)
    forward_site = site_of(result, "forward")
    check(forward_site["outcome"] == "binds" and forward_site["unknown_from_3prime"] == [],
          "primer Y fully accepts genome Y: a certain match")


# -- Case 21: graded 3' risk -------------------------------------------------
# a) Only the reverse site's 3'-terminal mismatch (amplicon offset 36 C -> T):
#    G-T is mild, so the product's 3' risk is LOW.
# b) The forward site's 3'-terminal base is N, tolerance 1: the site binds
#    even if the N is a mismatch, but the 3' end is uncertain: POSSIBLE.
def case_graded_risk():
    target = filler(10) + mutate(amplicon(), 36, "T") + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 21a: one mild mismatch at the 3' end", result)
    check(result["products"][0]["three_prime_risk"] == "low", "product 3' risk is LOW")

    target = filler(10) + mutate(amplicon(), 15, "N") + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 21b: N at the forward primer's 3' end", result)
    check(site_of(result, "forward")["outcome"] == "binds", "forward binds even if the N is a mismatch")
    check(result["products"][0]["three_prime_risk"] == "possible", "product 3' risk is POSSIBLE")


# -- Case 22: 1-base deletion in the primer's 5' half ------------------------
# Delete amplicon offset 2 (a T, 14 bases from the 3' end; the T beside it
# makes 15 an equally good answer). The 3' half is intact, so with 1 mismatch
# allowed the forward primer still binds, counting the gap as 1 mismatch.
# The product loses 1 base: 51 bp. With 0 mismatches allowed it does not bind.
def case_small_deletion_5prime():
    seq = amplicon()
    target = filler(10) + seq[:2] + seq[3:] + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 22: 1-base deletion in the 5' half, tolerance 1", result)
    forward_site = site_of(result, "forward")
    check((forward_site["indel"]["kind"], forward_site["indel"]["length"]) == ("deletion", 1),
          "forward binds through a 1-base deletion")
    check(forward_site["indel"]["from_3prime"][0] > 8, "the deletion is in the 5' half")
    check(product_lengths(result) == [51], "product is 51 bp (one base shorter)")

    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 22: same, tolerance 0", result)
    check(result["outcomes"]["forward"] == "mismatched", "with 0 allowed the gap is too much")


# -- Case 22b: 1-base deletion in the primer's 3' half -----------------------
# Delete amplicon offset 13 (a G, 3 bases from the 3' end). A gap that close
# to the 3' end is not allowed to bind; the near-miss search explains it.
def case_small_deletion_3prime():
    seq = amplicon()
    target = filler(10) + seq[:13] + seq[14:] + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 22b: 1-base deletion in the 3' half", result)
    check(result["outcomes"]["forward"] == "mismatched", "forward is 'mismatched'")
    check(result["near_misses"]["forward"]["indel"] ==
          {"kind": "deletion", "length": 1, "from_3prime": (3, 3)},
          "near-miss: 1-base deletion 3 bases from the 3' end")


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


def random_test_target(rng):
    """A 300-base random target with mutated primer copies, Ns/Ys and maybe a cut end."""
    target = [rng.choice("ACGT") for _ in range(300)]
    for _ in range(4):
        copy = list(FAKE_FORWARD if rng.random() < 0.5 else reverse_complement(FAKE_FORWARD))
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
# 200 random targets x mismatch allowances 0, 1 and 2. Plain (ungapped) sites
# from the model must equal the slow every-position search exactly.
def case_search_complete():
    print("\n--- Case 23: fast search vs slow every-position search ---")
    rng = random.Random(23)
    compared = 0
    for _ in range(200):
        target = random_test_target(rng)
        for max_mismatches in (0, 1, 2):
            result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=max_mismatches)
            fast = {(s["strand"], s["start"]) for s in result["sites"]
                    if s["primer"] == "forward" and s["indel"] is None}
            slow = brute_force_sites(FAKE_FORWARD, target, max_mismatches)
            if fast != slow:
                raise AssertionError(f"fast {sorted(fast)} != slow {sorted(slow)} for {target}")
            compared += len(slow)
    check(True, f"600 searches agree with the slow search ({compared} sites compared)")


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
    case_multiple_products()
    case_facing_away()
    case_no_site()
    case_n_in_site()
    case_n_cannot_change_answer()
    case_truncated_end()
    case_truncated_start()
    case_mismatched_beats_unknown()
    case_self_pairing()
    case_overlapping_sites()
    case_allowance_too_high()
    case_palindrome()
    case_degenerate_primer()
    case_probe_binds()
    case_probe_mismatched()
    case_deletion()
    case_deletion_too_long()
    case_insertion()
    case_genome_ambiguity()
    case_graded_risk()
    case_small_deletion_5prime()
    case_small_deletion_3prime()
    case_search_complete()
    case_random_genome()
    print("\nAll cases passed.")
