# =============================================================================
# TESTS FOR: PRIMER DESIGN   (version 6)
# =============================================================================
# The designer lives in primer_design.py; this file only tests it. Every case
# uses a FAKE genome and a FAKE target gene, built from seeded random bases
# (the same every run) or by hand, so each expected answer can be worked out
# on paper first.
#
# The basic fake genome:
#   2,000 random bases  +  400-base target  +  2,000 random bases
#   genome positions (0-based): target = 2000-2399
#
# Every design is also checked a second, independent way: run_model from the
# first model is run on the whole genome with the designed oligos, and must
# find exactly one product, at the designed place.
#
# Cases 1-15 test the GC rules and the specificity check WITHOUT temperature
# (annealing_c=None, as in version 3), so their expected answers are
# unchanged. Cases 16-20 test the temperature rules added in version 4;
# cases 21-22 the primer-dimer check and user conditions of version 5;
# case 23 the hairpin check of version 6; case 24 the structure_checks
# switch (Version 2 vs Version 3 of the website).
#
# Run it with:  python3 journal/test_primer_design.py
# =============================================================================

import random

from first_model_factor_specificity import reverse_complement, run_model
import glob
import os

from first_model_factor_specificity import read_fasta
from primer_design import (
    OLIGO_GC, PRIMER_LENGTH, PROBE_LENGTH, gc_fraction, match_penalty, temperature_penalty,
)
from primer_design import design_primers as design_with_temperature
import second_model_annealing_temp
from second_model_annealing_temp import (
    PROJECT, QPCR_CONDITIONS, assay_temperature_check, conditions_for, perfect_tm,
)


def design_primers(*args, **options):
    """Cases 1-15: the designer without temperature (version 3 behaviour)."""
    return design_with_temperature(*args, annealing_c=None, **options)


# =============================================================================
# SECTION 1: FAKE GENOMES AND TARGETS
# =============================================================================

BACKGROUND = 2000
TARGET_LENGTH = 400


def random_bases(n, seed, gc=0.5):
    """n seeded random bases with about this share of G/C."""
    rng = random.Random(seed)
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(n))


def basic_genome(target):
    """2,000 random bases + target + 2,000 random bases."""
    return random_bases(BACKGROUND, seed=1) + target + random_bases(BACKGROUND, seed=2)


TARGET = random_bases(TARGET_LENGTH, seed=3, gc=0.45)
GENOME = basic_genome(TARGET)

SHIFT = {"A": "C", "C": "G", "G": "T", "T": "A"}   # a base always changed to a different one


def check(condition, message):
    if not condition:
        raise AssertionError(message)
    print(f"  PASS: {message}")


def check_errors(genome, target, allowance, words, message, products=(60, 150)):
    """design_primers must stop with a ValueError whose text contains words."""
    try:
        design_primers(genome, target, allowance, *products)
    except ValueError as error:
        check(words in str(error), f"{message} (error: {error})")
        return
    raise AssertionError(f"{message}: no error raised")


def check_design(design, genome, target_start, target_end, allowance, products=(60, 150)):
    """Every rule a returned design must follow."""
    forward, reverse, probe = design["forward"], design["reverse"], design["probe"]
    return all([
        forward["seq"] == genome[forward["start"]:forward["end"]],
        reverse["seq"] == reverse_complement(genome[reverse["start"]:reverse["end"]]),
        probe["seq"] == genome[probe["start"]:probe["end"]],
        all(PRIMER_LENGTH[0] <= len(o["seq"]) <= PRIMER_LENGTH[1] for o in (forward, reverse)),
        PROBE_LENGTH[0] <= len(probe["seq"]) <= PROBE_LENGTH[1],
        target_start <= forward["start"] and reverse["end"] <= target_end,
        forward["end"] <= probe["start"] and probe["end"] <= reverse["start"],
        products[0] <= design["product_length"] <= products[1],
        design["product_length"] == reverse["end"] - forward["start"],
        all(OLIGO_GC[0] <= gc_fraction(o["seq"]) <= OLIGO_GC[1] for o in (forward, reverse, probe)),
        independent_check(design, genome, allowance),
    ])


def independent_check(design, genome, allowance):
    """run_model on the whole genome finds only the designed product."""
    result = run_model(design["forward"]["seq"], design["reverse"]["seq"], genome,
                       probe=design["probe"]["seq"], max_mismatches=allowance)
    products = result["products"]
    return (len(products) == 1 and products[0]["status"] == "product"
            and products[0]["start"] == design["product_start"]
            and products[0]["end"] == design["product_end"])


def no_overlap(designs):
    spans = sorted((d["product_start"], d["product_end"]) for d in designs)
    return all(end <= next_start for (_, end), (next_start, _) in zip(spans, spans[1:]))


# =============================================================================
# SECTION 2: CASES THAT SHOULD GIVE DESIGNS
# =============================================================================

def case_basic():
    print("\n--- Case 1: random 4,400-base genome, 400-base target ---")
    result = design_primers(GENOME, TARGET, 2)
    designs = result["designs"]
    check(result["target"] == {"start": 2000, "end": 2400, "strand": "+"},
          "target found at 2000-2399 on the + strand")
    check(len(designs) == 3, "3 designs")
    check(all(check_design(d, GENOME, 2000, 2400, 2) for d in designs),
          "every design follows every rule, and run_model agrees")
    check(no_overlap(designs), "the 3 products do not overlap")
    check(result["warnings"] == [], "no warnings (target GC about 45 %)")


def case_ranking():
    print("\n--- Case 2: ranking, best match first (worked out by hand) ---")
    # penalty = sum over the 3 oligos of |GC% - 50| + (product - shortest) / 10
    # with the default range 60-150:
    #   forward 9/20 = 45 % -> 5, reverse 10/20 = 50 % -> 0,
    #   probe 11/20 = 55 % -> 5, product 80 bp -> (80 - 60) / 10 = 2;  total 12
    forward = "AAAAAAAAAAAGGGGGGGGG"   # 9 G of 20
    reverse = "AAAAAAAAAAGGGGGGGGGG"   # 10 G of 20
    probe = "AAAAAAAAAGGGGGGGGGGG"     # 11 G of 20
    check(match_penalty(forward, reverse, probe, 80) == 12, "hand example: penalty 12")
    check(match_penalty(reverse, reverse, reverse, 60) == 0, "all 50 %, 60 bp: penalty 0")
    check(match_penalty(reverse, reverse, reverse, 150) == 9, "all 50 %, 150 bp: penalty 9")

    designs = design_primers(GENOME, TARGET, 2)["designs"]
    first = designs[0]
    check(first["match_penalty"] == 0 and first["product_length"] == 60,
          f"design 1: every oligo 50 % GC, 60 bp (the shortest allowed), penalty 0 "
          f"(forward {first['forward']['seq']})")
    penalties = [d["match_penalty"] for d in designs]
    check(penalties == sorted(penalties), f"designs come out best match first {penalties}")


def case_target_minus_strand():
    print("\n--- Case 3: target given as its reverse complement ---")
    plus = design_primers(GENOME, TARGET, 2)
    minus = design_primers(GENOME, reverse_complement(TARGET), 2)
    check(minus["target"]["strand"] == "-", "target found on the - strand")
    check(minus["designs"] == plus["designs"], "same designs as the + strand target")


def case_fasta_input():
    print("\n--- Case 4: genome and target pasted as FASTA, lowercase, line breaks ---")
    def as_fasta(name, seq):
        return f">{name}\n" + "\n".join(seq[i:i + 60].lower() for i in range(0, len(seq), 60))
    plain = design_primers(GENOME, TARGET, 2)
    fasta = design_primers(as_fasta("genome", GENOME), as_fasta("target", TARGET), 2)
    check(fasta["designs"] == plain["designs"], "same designs as plain uppercase input")


def case_n_in_target():
    print("\n--- Case 5: 30 N bases in the middle of the target ---")
    target = TARGET[:180] + "N" * 30 + TARGET[210:]
    genome = basic_genome(target)
    result = design_primers(genome, target, 2)
    oligos = [d[o]["seq"] for d in result["designs"] for o in ("forward", "reverse", "probe")]
    check(len(result["designs"]) > 0, f"{len(result['designs'])} designs still found")
    check(all("N" not in seq for seq in oligos), "no primer or probe contains N")


def case_gc_island():
    print("\n--- Case 6: a 40-base all-G/C stretch inside an ordinary target ---")
    target = TARGET[:100] + random_bases(40, seed=6, gc=1.0) + TARGET[140:]
    genome = basic_genome(target)
    result = design_primers(genome, target, 2)
    oligos = [d[o]["seq"] for d in result["designs"] for o in ("forward", "reverse", "probe")]
    check(len(result["designs"]) > 0, f"{len(result['designs'])} designs found")
    check(all(0.40 <= gc_fraction(seq) <= 0.60 for seq in oligos),
          "every primer and probe is 40-60 % GC (none from the G/C stretch)")


# =============================================================================
# SECTION 3: SPECIFICITY (the first model's job)
# =============================================================================

def case_two_copies():
    print("\n--- Case 7: the target appears twice in the genome ---")
    genome = basic_genome(TARGET) + TARGET
    result = design_primers(genome, TARGET, 2)
    check(result["designs"] == [], "no design: every primer binds in both copies")


def case_near_copy():
    print("\n--- Case 8: a second copy with every 12th base changed ---")
    # Any 18-24 bases in a row cover 1 or 2 changed bases (one in every 12),
    # so every primer matches the copy with at most 2 mismatches:
    # off-target at allowance 2. At allowance 0 the copy never counts.
    near_copy = "".join(SHIFT[b] if i % 12 == 0 else b for i, b in enumerate(TARGET))
    genome = basic_genome(TARGET) + near_copy
    at_2 = design_primers(genome, TARGET, 2)
    at_0 = design_primers(genome, TARGET, 0)
    check(at_2["designs"] == [], "allowance 2: no design (every primer also binds the copy)")
    check(len(at_0["designs"]) == 3, "allowance 0: 3 designs (the copy is 1-2 mismatches away)")
    check(all(check_design(d, genome, 2000, 2400, 0) for d in at_0["designs"]),
          "allowance 0 designs follow every rule, and run_model agrees")


# =============================================================================
# SECTION 4: GC CONTENT
# =============================================================================

def case_target_gc_too_high():
    print("\n--- Case 9: target 80 % GC ---")
    target = random_bases(TARGET_LENGTH, seed=9, gc=0.80)
    check_errors(basic_genome(target), target, 2, "GC content", "stops with a GC error")


def case_gc_checked_first():
    print("\n--- Case 10: GC-rich target that is not even in the genome ---")
    # If the GC check really runs first, the error is about GC, not about the
    # target being missing.
    target = random_bases(TARGET_LENGTH, seed=10, gc=0.80)
    check_errors(GENOME, target, 2, "GC content", "the GC error comes before the search")


def case_target_gc_medium():
    print("\n--- Case 11: target about 65 % GC ---")
    target = random_bases(TARGET_LENGTH, seed=11, gc=0.65)
    result = design_primers(basic_genome(target), target, 2)
    check(0.60 < gc_fraction(target) <= 0.70, f"target GC {gc_fraction(target):.1%} is medium")
    check(any("(medium" in w for w in result["warnings"]), "medium target warning given")
    check(len(result["designs"]) > 0, f"{len(result['designs'])} designs still returned")


def case_target_gc_low():
    print("\n--- Case 12: target about 33 % GC ---")
    target = random_bases(TARGET_LENGTH, seed=12, gc=0.33)
    result = design_primers(basic_genome(target), target, 2)
    check(gc_fraction(target) < 0.40, f"target GC {gc_fraction(target):.1%} is low")
    check(any("(low" in w for w in result["warnings"]), "low target warning, not an error")
    check(len(result["designs"]) > 0, f"{len(result['designs'])} designs still returned")


def case_no_gc_window():
    print("\n--- Case 13: target with no stretch at 40-60 % GC ---")
    # ATTA repeated: 0 % GC everywhere, so no primer can be chosen
    target = "ATTA" * 50
    result = design_primers(basic_genome(target), target, 2)
    check(result["possible"] == 0 and result["designs"] == [], "no possible design")
    check(any("40-60 % GC" in w for w in result["warnings"]), "the warning says GC is why")


def case_product_length():
    print("\n--- Case 14: the user chooses the product length range ---")
    short = design_primers(GENOME, TARGET, 2, 70, 90)
    check(all(check_design(d, GENOME, 2000, 2400, 2, (70, 90)) for d in short["designs"]),
          f"70-90 bp: every design follows every rule {[d['product_length'] for d in short['designs']]}")
    check(short["designs"][0]["product_length"] == 70, "70-90 bp: the best design is 70 bp (shortest)")
    check(short["warnings"] == [], "70-90 bp: no warning (inside 60-150)")
    long = design_primers(GENOME, TARGET, 2, 250, 350)
    check(all(check_design(d, GENOME, 2000, 2400, 2, (250, 350)) for d in long["designs"]),
          f"250-350 bp: every design follows every rule {[d['product_length'] for d in long['designs']]}")
    check(any("outside 60-150" in w for w in long["warnings"]), "250-350 bp: qPCR warning given")


# =============================================================================
# SECTION 5: BAD INPUT
# =============================================================================

def case_bad_input():
    print("\n--- Case 15: bad input ---")
    check_errors(GENOME, random_bases(TARGET_LENGTH, seed=14, gc=0.45), 2,
                 "not in the genome", "target not in the genome")
    check_errors(GENOME, TARGET, 3, "0 to 2", "allowance 3 is too high")
    check_errors(GENOME, TARGET, -1, "0 to 2", "allowance -1 is too low")
    check_errors(GENOME, TARGET, 2, "at least 54", "shortest product 53 is too short", (53, 150))
    check_errors(GENOME, TARGET, 2, "shorter than the shortest", "range 150-60 is backwards", (150, 60))
    check_errors(GENOME, TARGET[:59], 2, "shorter than the shortest product",
                 "a 59-base target cannot hold a 60 bp product")


# =============================================================================
# SECTION 6: RUN EVERY CASE
# =============================================================================

# =============================================================================
# SECTION 4: TEMPERATURE (version 4)
# -----------------------------------------------------------------------------
# With annealing at 55 C (the default) every design must pass the second
# model's qPCR rules (assay_temperature_check): primers 61-68 C and within
# 4 C of each other, probe 2 C above the higher primer, product Tm <= 90 C.
# =============================================================================

def temperature_ok(design, genome, annealing=55):
    tm = design["tm"]
    check_result = assay_temperature_check(design["forward"]["seq"], design["reverse"]["seq"],
                                           design["probe"]["seq"],
                                           genome[design["product_start"]:design["product_end"]],
                                           annealing_c=annealing)
    return check_result["ok"] and abs(check_result["forward"]["tm"] - tm["forward"]) < 1e-9


def case_temperature_basic():
    print("\n--- Case 16: with temperature, every design passes the qPCR rules ---")
    result = design_with_temperature(GENOME, TARGET, 2)
    designs = result["designs"]
    check(result["annealing_c"] == 55 and len(designs) == 3, "3 designs at the default 55 C")
    check(all(check_design(d, GENOME, 2000, 2400, 2) for d in designs),
          "every GC, length, position and specificity rule still holds")
    check(all(temperature_ok(d, GENOME) for d in designs),
          "every design passes assay_temperature_check "
          + str([(round(d['tm']['forward'], 1), round(d['tm']['reverse'], 1),
                  round(d['tm']['probe'], 1)) for d in designs]))


def case_temperature_penalty():
    print("\n--- Case 17: temperature ranking, worked out by hand ---")
    # annealing 55 C: primer target 63 C, probe target 70 C
    #   forward 62.0 -> 1.0, reverse 64.5 -> 1.5, probe 68.0 -> 2.0,
    #   product 80 bp -> (80 - 60) / 10 = 2.0;  total 6.5
    check(temperature_penalty(62.0, 64.5, 68.0, 80, 60, 55) == 6.5, "hand example: penalty 6.5")
    check(temperature_penalty(63.0, 63.0, 70.0, 60, 60, 55) == 0, "all on target, 60 bp: penalty 0")
    check(temperature_penalty(68.0, 68.0, 75.0, 60, 60, 60) == 0, "annealing 60 C: targets 68 and 75 C")
    designs = design_with_temperature(GENOME, TARGET, 2)["designs"]
    penalties = [d["match_penalty"] for d in designs]
    check(penalties == sorted(penalties), f"designs come out best first {penalties}")


def case_higher_annealing():
    print("\n--- Case 18: a higher annealing temperature asks for longer oligos ---")
    # 60 C annealing: primers are aimed at 68 C instead of 63 C, and at the
    # same GC a longer oligo has a higher Tm
    at_55 = design_with_temperature(GENOME, TARGET, 2, annealing_c=55)["designs"]
    at_60 = design_with_temperature(GENOME, TARGET, 2, annealing_c=60)["designs"]
    mean = lambda designs: sum(len(d[k]["seq"]) for d in designs for k in ("forward", "reverse")) / (2 * len(designs))
    check(at_60 and mean(at_60) > mean(at_55),
          f"mean primer length {mean(at_55):.1f} nt at 55 C < {mean(at_60):.1f} nt at 60 C")
    check(all(temperature_ok(d, GENOME, 60) for d in at_60), "every 60 C design passes the rules at 60 C")


def case_n_gene_temperature():
    print("\n--- Case 19: the N gene: version 3's low-Tm picks are gone ---")
    # Version 3 picked TGTCTGATAATGGACCCC (Tm 58.1 C at 500 nM), only 3 C above
    # the 55 C annealing step; Primer-BLAST would not pick it
    # (test_cases_web_specificity.txt, part C).
    reference = read_fasta(glob.glob(os.path.join(PROJECT, "data", "reference", "*"))[0])[1]
    gene = reference[28274 - 1:29533]
    old = design_primers(reference, gene, 2)["designs"]
    new = design_with_temperature(reference, gene, 2)["designs"]
    check(old[0]["forward"]["seq"] == "TGTCTGATAATGGACCCC", "version 3's first forward primer, for comparison")
    check(len(new) == 3 and all(temperature_ok(d, reference) for d in new),
          "version 4: 3 designs, all pass the qPCR rules")
    check(not {d["forward"]["seq"] for d in old} & {d["forward"]["seq"] for d in new},
          "none of version 3's forward primers is kept")
    blast = {"monovalent_mM": 50.0, "divalent_mM": 1.5, "dntp_mM": 0.6, "oligo_nM": 50.0}
    tms = [perfect_tm(d[k]["seq"], blast) for d in new for k in ("forward", "reverse")]
    check(all(57 <= tm <= 63 for tm in tms),
          f"at Primer-BLAST's conditions every primer is inside its 57-63 C window {[round(t, 1) for t in tms]}")


def case_annealing_errors():
    print("\n--- Case 20: annealing temperature out of range ---")
    for value in (40, 80):
        try:
            design_with_temperature(GENOME, TARGET, 2, annealing_c=value)
        except ValueError as error:
            check("annealing temperature" in str(error), f"{value} C: error ({error})")
        else:
            raise AssertionError(f"{value} C: no error")



def case_dimer_filter():
    print("\n--- Case 21: designs whose oligos form dimers are skipped ---")
    normal = design_with_temperature(GENOME, TARGET, 2)
    check(all(min(d["any"]["dG"] for d in design["dimers"]) >= -9000 for design in normal["designs"]),
          f"default limit -9 kcal/mol: every design passes ({normal['dimer_rejected']} skipped)")
    # a much stricter limit for this test only, so that some designs must be skipped
    saved = second_model_annealing_temp.DIMER_MAX_DG
    second_model_annealing_temp.DIMER_MAX_DG = -4000.0
    try:
        strict = design_with_temperature(GENOME, TARGET, 2)
    finally:
        second_model_annealing_temp.DIMER_MAX_DG = saved
    check(strict["dimer_rejected"] > 0, f"limit -4 kcal/mol: {strict['dimer_rejected']} designs skipped")
    check(all(min(d["any"]["dG"] for d in design["dimers"]) >= -4000 for design in strict["designs"]),
          "every design kept has no dimer more stable than -4 kcal/mol")


def case_hairpin_filter():
    print("\n--- Case 23: designs with a folding oligo are skipped ---")
    normal = design_with_temperature(GENOME, TARGET, 2)
    check(all(f is None or f["dG"] >= -3000 for d in normal["designs"] for f in d["hairpins"].values()),
          f"default limit -3 kcal/mol: no design kept has a stronger hairpin ({normal['hairpin_rejected']} skipped)")
    saved = second_model_annealing_temp.HAIRPIN_MAX_DG
    second_model_annealing_temp.HAIRPIN_MAX_DG = 0.0   # any fold more stable than 0: for this test only
    try:
        strict = design_with_temperature(GENOME, TARGET, 2)
    finally:
        second_model_annealing_temp.HAIRPIN_MAX_DG = saved
    check(strict["hairpin_rejected"] > 0, f"limit 0 kcal/mol: {strict['hairpin_rejected']} designs skipped")
    check(all(f is None or f["dG"] >= 0 for d in strict["designs"] for f in d["hairpins"].values()),
          "every design kept has no hairpin more stable than 0 kcal/mol")


def case_structure_switch():
    print("\n--- Case 24: Version 2 (no dimer or hairpin checks) vs Version 3 on the N gene ---")
    reference = read_fasta(glob.glob(os.path.join(PROJECT, "data", "reference", "*"))[0])[1]
    gene = reference[28274 - 1:29533]
    v2 = design_with_temperature(reference, gene, 2, structure_checks=False)
    v3 = design_with_temperature(reference, gene, 2)
    check(v2["designs"][0]["reverse"]["seq"] == "ATCTGTCAAGCAGCAGCAA" and "dimers" not in v2["designs"][0],
          "Version 2 keeps the CAG-repeat design (reverse ATCTGTCAAGCAGCAGCAA), unchecked for dimers")
    check(v3["dimer_rejected"] >= 1 and v3["designs"][0]["reverse"]["seq"] != "ATCTGTCAAGCAGCAGCAA",
          "Version 3 skips it: its reverse/probe dimer is -9.2 kcal/mol")


def case_user_conditions():
    print("\n--- Case 22: the user's salt and oligo amounts are used ---")
    mix = dict(QPCR_CONDITIONS, divalent_mM=5.0, primer_nM=250.0)
    result = design_with_temperature(GENOME, TARGET, 2, conditions=mix)
    check(result["conditions"]["divalent_mM"] == 5.0, "the result records the conditions used")
    design = result["designs"][0]
    check(abs(design["tm"]["forward"] - perfect_tm(design["forward"]["seq"], conditions_for("forward", mix))) < 1e-9,
          "the forward Tm is computed at 5 mM Mg2+ and 250 nM primer")
    check(all(temperature_ok_with(d, GENOME, mix) for d in result["designs"]),
          "every design passes the qPCR rules at those conditions")
    try:
        design_with_temperature(GENOME, TARGET, 2, conditions=dict(QPCR_CONDITIONS, divalent_mM=20.0))
    except ValueError as error:
        check("divalent" in str(error), f"20 mM Mg2+ refused ({error})")
    else:
        raise AssertionError("20 mM Mg2+ accepted")


def temperature_ok_with(design, genome, conditions):
    result = assay_temperature_check(design["forward"]["seq"], design["reverse"]["seq"], design["probe"]["seq"],
                                     genome[design["product_start"]:design["product_end"]],
                                     conditions=conditions)
    return result["ok"]


if __name__ == "__main__":
    case_basic()
    case_ranking()
    case_target_minus_strand()
    case_fasta_input()
    case_n_in_target()
    case_gc_island()
    case_two_copies()
    case_near_copy()
    case_target_gc_too_high()
    case_gc_checked_first()
    case_target_gc_medium()
    case_target_gc_low()
    case_no_gc_window()
    case_product_length()
    case_bad_input()
    case_temperature_basic()
    case_temperature_penalty()
    case_higher_annealing()
    case_n_gene_temperature()
    case_annealing_errors()
    case_dimer_filter()
    case_user_conditions()
    case_hairpin_filter()
    case_structure_switch()
    print("\nAll cases passed.")
