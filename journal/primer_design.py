# =============================================================================
# PRIMER DESIGN: a qPCR primer pair + probe for a target gene   (version 1)
# =============================================================================
# 02-120 project: Do early SARS-CoV-2 qPCR primers still match later variants?
#
# The first model CHECKS primers someone else designed. This file goes the
# other way: given a genome and a target gene inside it, it DESIGNS a forward
# primer, a reverse primer and a probe, using only what the first model
# (first_model_factor_specificity.py) knows about: sequences, the mismatch
# allowance, binding sites on both strands, product length, and whether the
# probe lies inside the product.
#
# It does NOT use melting temperature or annealing temperature. It does use
# GC content (Section 2b), because G-C pairs hold three hydrogen bonds to A-T's
# two: a GC-rich primer binds too strongly, also at sites it only partly
# matches, and a GC-rich template folds on itself and copies poorly.
#
# Inputs:
#   genome          - the full genome sequence (a SARS-CoV-2 genome is ~30 kb),
#                     plain or FASTA
#   target          - the target gene's sequence; it must appear in the genome
#                     (on either strand), e.g. the N gene
#   max_mismatches  - the mismatch allowance, 0 to MAX_ALLOWED_MISMATCHES (2)
#
# Output, for each design (best first):
#   5'-forward-3', and 3'-reverse-5' (the reverse primer written backwards, so
#   it lines up under the genome the way it binds); where each primer binds on
#   the genome and its length; the product length; the probe and where it
#   binds.
#
# How it works:
#   0. check the target's GC content FIRST: too high -> stop with an error
#      and the reason; medium -> carry on, with a warning
#   1. find the target in the genome
#   2. every PRIMER_LENGTH-base window of the target made of known bases
#      (A/C/G/T only), with GC content not too high, is a possible forward
#      primer, and its reverse complement a possible reverse primer
#   3. pair each forward with each reverse downstream of it that gives a
#      product of PRODUCT_LENGTH, and put a PROBE_LENGTH probe in the middle
#      of the gap between them
#   4. walk the designs from the fewest medium-GC oligos, then the product
#      length closest to IDEAL_PRODUCT; keep
#      a design only if, on the WHOLE genome and within the mismatch allowance,
#        - each oligo has exactly one binding site, and
#        - run_model finds exactly one product, status "product", at the
#          designed place, with the probe inside it
#
# Run the demo with:  python3 journal/primer_design.py
#   (designs primers for the N gene of the reference genome)
# Or with your own files:
#   python3 journal/primer_design.py genome.fasta target.fasta 2
# =============================================================================

import glob
import os
import sys

from first_model_factor_specificity import (
    KNOWN_BASES, build_kmer_index, find_binding_sites, read_fasta, reverse_complement,
    run_model,
)


# =============================================================================
# SECTION 1: SETTINGS
# -----------------------------------------------------------------------------
# Lengths follow the CDC N1/N2/N3 oligos (primers 20-24 nt, probes 23-24 nt,
# products 67-72 bp).
#
# The mismatch allowance is capped at 2. The first model finds sites with
# k = 6 base seeds, and needs (oligo length) // (allowance + 1) >= 6, so a
# 20-base primer can use at most 2 (#15, #40). A higher allowance could miss
# off-target sites without warning.
# =============================================================================

K = 6
GC_MEDIUM_ABOVE = 0.60      # above this: still used, but flagged with a warning
GC_TOO_HIGH_ABOVE = 0.70    # above this: an error for the target, never used for an oligo
PRIMER_LENGTH = 20
PROBE_LENGTH = 24
PRODUCT_LENGTH = (60, 150)
IDEAL_PRODUCT = 100
MAX_ALLOWED_MISMATCHES = 2
DESIGNS_WANTED = 3


# =============================================================================
# SECTION 2: READING THE INPUT
# =============================================================================

def clean_sequence(text):
    """Accept plain sequence or FASTA text: drop '>' lines, spaces and line breaks."""
    lines = [line for line in text.splitlines() if not line.startswith(">")]
    return "".join("".join(lines).split()).upper()


def find_target(genome, target):
    """Return (start, end, strand) of the target in the genome, 0-based, end exclusive."""
    start = genome.find(target)
    if start >= 0:
        return start, start + len(target), "+"
    start = genome.find(reverse_complement(target))
    if start >= 0:
        return start, start + len(target), "-"
    raise ValueError("the target gene is not in the genome (checked both strands)")


# =============================================================================
# SECTION 2b: GC CONTENT
# -----------------------------------------------------------------------------
# GC content = the share of bases that are G or C. Three levels:
#   ok        - up to 60 %
#   medium    - above 60 %, up to 70 %: binds more strongly than it should;
#               still used, but flagged with a warning
#   too high  - above 70 %: binds too strongly, also where it only partly
#               matches, and the template folds on itself. A target this
#               GC-rich stops the design; an oligo this GC-rich is never used.
# Only unambiguous bases count, so missing data (N) does not lower the GC.
# =============================================================================

def gc_fraction(seq):
    known_bases = [base for base in seq if base in KNOWN_BASES]
    if not known_bases:
        return 0.0
    return sum(base in "GC" for base in known_bases) / len(known_bases)


def gc_level(seq):
    """Return 'ok', 'medium' or 'too high' for this sequence's GC content."""
    gc = gc_fraction(seq)
    if gc > GC_TOO_HIGH_ABOVE:
        return "too high"
    if gc > GC_MEDIUM_ABOVE:
        return "medium"
    return "ok"


def gc_percent(seq):
    return f"{100 * gc_fraction(seq):.1f} %"


def check_target_gc(target):
    """Stop with the reason if the target is too GC-rich; return a warning (or None) if medium."""
    level = gc_level(target)
    if level == "too high":
        raise ValueError(
            f"the target's GC content is {gc_percent(target)}, above "
            f"{round(100 * GC_TOO_HIGH_ABOVE)} %: primers here would bind too strongly, "
            f"also at partly matching sites, and the GC-rich template folds on itself "
            f"and copies poorly. Choose a less GC-rich target.")
    if level == "medium":
        return (f"target GC content is {gc_percent(target)} (medium, above "
                f"{round(100 * GC_MEDIUM_ABOVE)} %): primers may bind more strongly than "
                f"intended; check the designs before ordering")
    return None


# =============================================================================
# SECTION 3: EVERY POSSIBLE DESIGN IN THE TARGET
# -----------------------------------------------------------------------------
# A "window" is a stretch of the genome's + strand. The forward primer IS its
# window; the reverse primer is the reverse complement of its window (so it
# binds the + strand, and its 5' end is at the window's right end). The probe
# is written as its + strand window.
#
# Windows with N or ambiguous bases are skipped: the first model cannot say
# whether a primer binds there ("not assessable").
# =============================================================================

def known(seq):
    return all(base in KNOWN_BASES for base in seq)


def usable(seq):
    """Known bases only, and GC content not too high."""
    return known(seq) and gc_level(seq) != "too high"


def possible_designs(genome, region_start, region_end):
    """Every forward/reverse/probe design inside the target, best first."""
    starts = [s for s in range(region_start, region_end - PRIMER_LENGTH + 1)
              if usable(genome[s:s + PRIMER_LENGTH])]
    designs = []
    for forward_start in starts:
        for reverse_start in starts:
            product_end = reverse_start + PRIMER_LENGTH
            length = product_end - forward_start
            if not PRODUCT_LENGTH[0] <= length <= PRODUCT_LENGTH[1]:
                continue
            gap_start, gap_end = forward_start + PRIMER_LENGTH, reverse_start
            if gap_end - gap_start < PROBE_LENGTH:
                continue  # no room for the probe between the primers
            probe_start = gap_start + (gap_end - gap_start - PROBE_LENGTH) // 2
            probe = genome[probe_start:probe_start + PROBE_LENGTH]
            if not usable(probe):
                continue
            medium_gc = [name for name, seq in (
                ("forward", genome[forward_start:forward_start + PRIMER_LENGTH]),
                ("reverse", genome[reverse_start:product_end]),
                ("probe", probe)) if gc_level(seq) == "medium"]
            designs.append({
                "forward": {"seq": genome[forward_start:forward_start + PRIMER_LENGTH],
                            "start": forward_start, "end": forward_start + PRIMER_LENGTH},
                "reverse": {"seq": reverse_complement(genome[reverse_start:product_end]),
                            "start": reverse_start, "end": product_end},
                "probe": {"seq": probe, "start": probe_start,
                          "end": probe_start + PROBE_LENGTH},
                "product_start": forward_start,
                "product_end": product_end,
                "product_length": length,
                "medium_gc": medium_gc,
            })
    designs.sort(key=lambda d: (len(d["medium_gc"]), abs(d["product_length"] - IDEAL_PRODUCT),
                                d["product_start"]))
    return designs


# =============================================================================
# SECTION 4: SPECIFICITY CHECK WITH THE FIRST MODEL
# -----------------------------------------------------------------------------
#   4a. quick screen: each oligo must have exactly one binding site on the
#       whole genome, on either strand, within the mismatch allowance
#       (outcome "binds" or "not assessable" both count as a site). One shared
#       k-mer index; each oligo is checked once even if it is in many designs.
#   4b. full check: run_model on the whole genome must give exactly one
#       product, status "product", at the designed position.
# =============================================================================

def binds_once(oligo, genome, index, max_mismatches, cache):
    if oligo not in cache:
        sites = find_binding_sites("oligo", oligo, genome, index, K, max_mismatches)
        cache[oligo] = len(sites) == 1
    return cache[oligo]


def specific(design, genome, index, max_mismatches, cache):
    """True if the first model finds only the designed product on the whole genome."""
    oligos = (design["forward"]["seq"], design["reverse"]["seq"], design["probe"]["seq"])
    if not all(binds_once(oligo, genome, index, max_mismatches, cache) for oligo in oligos):
        return False
    result = run_model(oligos[0], oligos[1], genome, probe=oligos[2],
                       max_mismatches=max_mismatches, k=K)
    products = result["products"]
    return (len(products) == 1 and products[0]["status"] == "product"
            and products[0]["start"] == design["product_start"]
            and products[0]["end"] == design["product_end"])


# =============================================================================
# SECTION 5: THE WHOLE DESIGN IN ONE CALL
# =============================================================================

def design_primers(genome, target, max_mismatches=2, designs_wanted=DESIGNS_WANTED):
    """Design up to designs_wanted primer pairs + probes for target, best first."""
    genome, target = clean_sequence(genome), clean_sequence(target)
    if not 0 <= max_mismatches <= MAX_ALLOWED_MISMATCHES:
        raise ValueError(f"mismatch allowance must be 0 to {MAX_ALLOWED_MISMATCHES}")
    if len(target) < PRODUCT_LENGTH[0]:
        raise ValueError(f"the target must be at least {PRODUCT_LENGTH[0]} bases long")

    # -- the GC check comes first: a too GC-rich target stops here
    warnings = []
    target_warning = check_target_gc(target)
    if target_warning:
        warnings.append(target_warning)

    region_start, region_end, strand = find_target(genome, target)
    designs = possible_designs(genome, region_start, region_end)

    index = build_kmer_index(genome, K)
    cache, chosen, checked = {}, [], 0
    for design in designs:
        if len(chosen) == designs_wanted:
            break
        # each design's product must not overlap an earlier one's, so the
        # options cover different parts of the gene
        if any(design["product_start"] < d["product_end"] and d["product_start"] < design["product_end"]
               for d in chosen):
            continue
        checked += 1
        if specific(design, genome, index, max_mismatches, cache):
            chosen.append(design)

    for number, design in enumerate(chosen, 1):
        for name in design["medium_gc"]:
            seq = design[name]["seq"]
            warnings.append(f"design {number}: {name} GC content is {gc_percent(seq)} "
                            f"(medium): it may bind more strongly than intended")

    return {
        "warnings": warnings,
        "target": {"start": region_start, "end": region_end, "strand": strand},
        "possible": len(designs),
        "checked": checked,
        "designs": chosen,
    }


# =============================================================================
# SECTION 6: PRINTING (positions printed 1-based, like the CDC table)
# =============================================================================

def print_designs(result, max_mismatches):
    target = result["target"]
    print(f"Target: genome positions {target['start'] + 1}-{target['end']} "
          f"({target['strand']} strand), {target['end'] - target['start']} bases")
    print(f"{result['possible']} possible designs; {result['checked']} checked against "
          f"the whole genome (mismatch allowance {max_mismatches}).\n")
    for warning in result["warnings"]:
        print(f"WARNING: {warning}")
    if result["warnings"]:
        print()
    if not result["designs"]:
        print("No design binds only at its own site.")
        return
    for number, design in enumerate(result["designs"], 1):
        forward, reverse, probe = design["forward"], design["reverse"], design["probe"]
        flag = f"   (medium GC: {', '.join(design['medium_gc'])})" if design["medium_gc"] else ""
        print(f"Design {number}{flag}")
        print(f"  forward  5'-{forward['seq']}-3'   {len(forward['seq'])} nt, GC {gc_percent(forward['seq'])}")
        print(f"           binds {forward['start'] + 1}-{forward['end']} (starts at {forward['start'] + 1})")
        print(f"  reverse  3'-{reverse['seq'][::-1]}-5'   {len(reverse['seq'])} nt, GC {gc_percent(reverse['seq'])}")
        print(f"           binds {reverse['start'] + 1}-{reverse['end']} (5' end at {reverse['end']})")
        print(f"           order it as 5'-{reverse['seq']}-3'")
        print(f"  probe    5'-{probe['seq']}-3'   {len(probe['seq'])} nt, GC {gc_percent(probe['seq'])}")
        print(f"           binds {probe['start'] + 1}-{probe['end']}")
        print(f"  product  {design['product_start'] + 1}-{design['product_end']}, "
              f"{design['product_length']} bp, the only product on the whole genome")
        print()


# =============================================================================
# SECTION 7: RUN
# -----------------------------------------------------------------------------
# With no arguments: the N gene of the reference genome NC_045512.2
# (positions 28274-29533, the gene the CDC N1/N2/N3 assays target).
# =============================================================================

if __name__ == "__main__":
    if len(sys.argv) >= 3:
        genome = read_fasta(sys.argv[1])[1]
        target = read_fasta(sys.argv[2])[1]
        allowance = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    else:
        project = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        genome = read_fasta(glob.glob(os.path.join(project, "data", "reference", "*"))[0])[1]
        target = genome[28274 - 1:29533]
        allowance = 2
        print("Demo: N gene of the reference genome NC_045512.2\n")

    try:
        result = design_primers(genome, target, allowance)
    except ValueError as error:
        print(f"ERROR: {error}")
        sys.exit(1)
    print_designs(result, allowance)
