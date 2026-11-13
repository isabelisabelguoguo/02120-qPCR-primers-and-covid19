# =============================================================================
# TESTS FOR: PRIMER DESIGN   (version 1)
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
# Run it with:  python3 journal/test_primer_design.py
# =============================================================================

import random

from first_model_factor_specificity import reverse_complement, run_model
from primer_design import (
    GC_TOO_HIGH_ABOVE, PRIMER_LENGTH, PROBE_LENGTH, PRODUCT_LENGTH, design_primers, gc_fraction,
)


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


def check_errors(genome, target, allowance, words, message):
    """design_primers must stop with a ValueError whose text contains words."""
    try:
        design_primers(genome, target, allowance)
    except ValueError as error:
        check(words in str(error), f"{message} (error: {error})")
        return
    raise AssertionError(f"{message}: no error raised")


def check_design(design, genome, target_start, target_end, allowance):
    """Every rule a returned design must follow."""
    forward, reverse, probe = design["forward"], design["reverse"], design["probe"]
    return all([
        forward["seq"] == genome[forward["start"]:forward["end"]],
        reverse["seq"] == reverse_complement(genome[reverse["start"]:reverse["end"]]),
        probe["seq"] == genome[probe["start"]:probe["end"]],
        len(forward["seq"]) == PRIMER_LENGTH and len(reverse["seq"]) == PRIMER_LENGTH,
        len(probe["seq"]) == PROBE_LENGTH,
        target_start <= forward["start"] and reverse["end"] <= target_end,
        forward["end"] <= probe["start"] and probe["end"] <= reverse["start"],
        PRODUCT_LENGTH[0] <= design["product_length"] <= PRODUCT_LENGTH[1],
        design["product_length"] == reverse["end"] - forward["start"],
        all(gc_fraction(o["seq"]) <= GC_TOO_HIGH_ABOVE for o in (forward, reverse, probe)),
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


def case_hand_solved_first_design():
    print("\n--- Case 2: the first design, worked out by hand ---")
    # Best = product closest to 100 bp, then earliest start: forward at the
    # target's first base (2000), product 2000-2099, reverse window 2080-2099.
    # Gap between primers = 2020-2079 (60 bases); 24-base probe centred:
    # 2020 + (60 - 24) // 2 = 2038, so the probe is 2038-2061.
    first = design_primers(GENOME, TARGET, 2)["designs"][0]
    check(first["forward"]["start"] == 2000, "forward starts at 2000")
    check(first["reverse"]["start"] == 2080 and first["reverse"]["end"] == 2100,
          "reverse window 2080-2099")
    check(first["probe"]["start"] == 2038, "probe starts at 2038, centred in the gap")
    check(first["product_length"] == 100, "product is 100 bp")


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
    check(gc_fraction(target) <= 0.60, f"target GC {gc_fraction(target):.1%} is ok overall")
    check(len(result["designs"]) > 0, f"{len(result['designs'])} designs found")
    check(all(gc_fraction(seq) <= GC_TOO_HIGH_ABOVE for seq in oligos),
          "no primer or probe above 70 % GC")


# =============================================================================
# SECTION 3: SPECIFICITY (the first model's job)
# =============================================================================

def case_two_copies():
    print("\n--- Case 7: the target appears twice in the genome ---")
    genome = basic_genome(TARGET) + TARGET
    result = design_primers(genome, TARGET, 2)
    check(result["designs"] == [], "no design: every primer binds in both copies")


def case_near_copy():
    print("\n--- Case 8: a second copy with every 10th base changed ---")
    # Any 20-base primer covers exactly 2 changed bases of the copy, so it
    # binds the copy with 2 mismatches: off-target at allowance 2, but not at
    # allowance 1.
    near_copy = "".join(SHIFT[b] if i % 10 == 0 else b for i, b in enumerate(TARGET))
    genome = basic_genome(TARGET) + near_copy
    at_2 = design_primers(genome, TARGET, 2)
    at_1 = design_primers(genome, TARGET, 1)
    check(at_2["designs"] == [], "allowance 2: no design (every primer also binds the copy)")
    check(len(at_1["designs"]) == 3, "allowance 1: 3 designs (the copy is 2 mismatches away)")
    check(all(check_design(d, genome, 2000, 2400, 1) for d in at_1["designs"]),
          "allowance 1 designs follow every rule, and run_model agrees")


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
    check(any("target GC content" in w for w in result["warnings"]), "target warning given")
    check(len(result["designs"]) > 0, f"{len(result['designs'])} designs still returned")
    flagged = sum(len(d["medium_gc"]) for d in result["designs"])
    check(sum(w.startswith("design ") for w in result["warnings"]) == flagged,
          f"one warning per medium-GC oligo ({flagged})")


# =============================================================================
# SECTION 5: BAD INPUT
# =============================================================================

def case_bad_input():
    print("\n--- Case 12: bad input ---")
    check_errors(GENOME, random_bases(TARGET_LENGTH, seed=12, gc=0.45), 2,
                 "not in the genome", "target not in the genome")
    check_errors(GENOME, TARGET, 3, "0 to 2", "allowance 3 is too high")
    check_errors(GENOME, TARGET, -1, "0 to 2", "allowance -1 is too low")
    check_errors(GENOME, TARGET[:59], 2, "at least 60", "a 59-base target is too short")


# =============================================================================
# SECTION 6: RUN EVERY CASE
# =============================================================================

if __name__ == "__main__":
    case_basic()
    case_hand_solved_first_design()
    case_target_minus_strand()
    case_fasta_input()
    case_n_in_target()
    case_gc_island()
    case_two_copies()
    case_near_copy()
    case_target_gc_too_high()
    case_gc_checked_first()
    case_target_gc_medium()
    case_bad_input()
    print("\nAll cases passed.")
