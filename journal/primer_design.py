# =============================================================================
# PRIMER DESIGN: a qPCR primer pair + probe for a target gene   (version 2)
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
# two. Keeping each primer and probe at 40-60 % GC gives stable binding: a
# GC-rich primer binds too strongly, also at sites it only partly matches; an
# AT-rich one binds too weakly. A GC-rich template also folds on itself and
# copies poorly.
#
# Inputs:
#   genome          - the full genome sequence (a SARS-CoV-2 genome is ~30 kb),
#                     plain or FASTA
#   target          - the target gene's sequence; it must appear in the genome
#                     (on either strand), e.g. the N gene
#   max_mismatches  - the mismatch allowance, 0 to MAX_ALLOWED_MISMATCHES (2)
#   ideal_product   - the product length to aim for (default 100 bp); designs
#                     are made within +/- 50 bp of it
#
# Output, for each design (best first):
#   5'-forward-3', and 3'-reverse-5' (the reverse primer written backwards, so
#   it lines up under the genome the way it binds); where each primer binds on
#   the genome and its length; the product length; the probe and where it
#   binds.
#
# How it works:
#   0. check the target's GC content FIRST: too high -> stop with an error
#      and the reason; medium or low -> carry on, with a warning
#   1. find the target in the genome
#   2. every position in the target can start a forward primer (or end a
#      reverse primer, or start a probe). Its length (primers 18-24, probes
#      18-30) is the one with GC closest to 50 %, among windows made of known
#      bases (A/C/G/T only) with 40-60 % GC
#   3. pair each forward with each reverse downstream of it that gives a
#      product within PRODUCT_SPREAD (50 bp) of ideal_product, and put the probe whose middle is closest
#      to the middle of the gap between them
#   4. rank the designs, best match first (Section 3b), and walk down the
#      ranking; keep a design only if, on the WHOLE genome and within the
#      mismatch allowance,
#        - each oligo has exactly one binding site, and
#        - run_model finds exactly one product, status "product", at the
#          designed place, with the probe inside it
#
# Run the demo with:  python3 journal/primer_design.py
#   (designs primers for the N gene of the reference genome)
# Or with your own files:
#   python3 journal/primer_design.py genome.fasta target.fasta 2 100
#   (the last two, allowance and product length, can be left out)
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
# Lengths: primers 18-24 bases (the usual qPCR range; the CDC primers are
# 18-24), probes 18-30 (a probe is often a little longer than the primers;
# the CDC probes are 23-24).
#
# A range alone would give every position 7 primer lengths and 13 probe
# lengths, so millions of designs and nothing to choose between them. Each
# position therefore gets ONE length: the one whose GC content is closest to
# IDEAL_GC (50 %), the shorter one on a tie (Section 3, choose_length).
#
# The mismatch allowance is capped at 2. The first model finds sites with
# k = 6 base seeds, and needs (oligo length) // (allowance + 1) >= 6, so the
# shortest oligo, 18 bases, can use at most 2 (#15, #40). A higher allowance
# could miss off-target sites without warning.
# =============================================================================

K = 6
OLIGO_GC = (0.40, 0.60)     # every primer and probe must be inside this range
IDEAL_GC = 0.50             # each oligo's length is chosen to get closest to this
TARGET_GC_LOW_BELOW = 0.40  # target below this: warning (AT-rich)
GC_MEDIUM_ABOVE = 0.60      # target above this: warning (GC-rich)
GC_TOO_HIGH_ABOVE = 0.70    # target above this: error, the design stops
PRIMER_LENGTH = (18, 24)
PROBE_LENGTH = (18, 30)
IDEAL_PRODUCT = 100         # default product length; the user can choose another
PRODUCT_SPREAD = 50         # designs are made within +/- this many bp of it
MIN_PRODUCT = 54            # room for an 18-base forward, probe and reverse
QPCR_PRODUCT = (60, 150)    # qPCR works best here; outside it gives a warning
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
# GC content = the share of bases that are G or C. Only unambiguous bases
# count, so missing data (N) does not lower the GC.
#
# Primers and probe: must be 40-60 % GC (OLIGO_GC), for stable hydrogen
# bonding. Anything outside is never used.
#
# Target gene, checked first, before anything else runs:
#   too high  - above 70 %: error and the reason; the design stops
#   medium    - above 60 %, up to 70 %: warning, primers may bind too strongly
#   low       - below 40 %: warning only, because SARS-CoV-2 as a whole is
#               about 38 % GC; 40-60 % stretches can usually still be found
#               inside an AT-rich gene
#   ok        - 40-60 %
# =============================================================================

def gc_fraction(seq):
    known_bases = [base for base in seq if base in KNOWN_BASES]
    if not known_bases:
        return 0.0
    return sum(base in "GC" for base in known_bases) / len(known_bases)


def gc_level(seq):
    """Return 'too high', 'medium', 'low' or 'ok' for a target's GC content."""
    gc = gc_fraction(seq)
    if gc > GC_TOO_HIGH_ABOVE:
        return "too high"
    if gc > GC_MEDIUM_ABOVE:
        return "medium"
    if gc < TARGET_GC_LOW_BELOW:
        return "low"
    return "ok"


def gc_percent(seq):
    return f"{100 * gc_fraction(seq):.1f} %"


def check_target_gc(target):
    """Stop with the reason if the target is too GC-rich; return a warning (or None)."""
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
    if level == "low":
        return (f"target GC content is {gc_percent(target)} (low, below "
                f"{round(100 * TARGET_GC_LOW_BELOW)} %): fewer stretches reach 40-60 % GC, "
                f"so there may be fewer designs")
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
    """Known bases only, and GC content 40-60 %."""
    return known(seq) and OLIGO_GC[0] <= gc_fraction(seq) <= OLIGO_GC[1]


def choose_length(genome, anchor, lengths, region_start, region_end, grow):
    """The one window for this anchor: the usable length with GC closest to 50 %.

    grow "right": the window is anchor .. anchor + length (forward primer, probe)
    grow "left":  the window is anchor - length .. anchor (reverse primer)
    Ties go to the shorter length. Returns (start, end), or None if no length works.
    """
    best = None
    for length in range(lengths[0], lengths[1] + 1):
        start, end = (anchor, anchor + length) if grow == "right" else (anchor - length, anchor)
        if start < region_start or end > region_end:
            break
        window = genome[start:end]
        if not usable(window):
            continue
        distance = abs(gc_fraction(window) - IDEAL_GC)
        if best is None or distance < best[0]:   # strictly closer, so ties keep the shorter
            best = (distance, start, end)
    return None if best is None else best[1:]


def all_windows(genome, anchors, lengths, region_start, region_end, grow):
    """{anchor: (start, end)} for every anchor that has a usable length."""
    windows = {}
    for anchor in anchors:
        window = choose_length(genome, anchor, lengths, region_start, region_end, grow)
        if window:
            windows[anchor] = window
    return windows


def centred_probe(probes, gap_start, gap_end):
    """The probe window inside the gap whose middle is closest to the gap's middle."""
    middle = (gap_start + gap_end) / 2
    best = None
    for start in range(gap_start, gap_end - PROBE_LENGTH[0] + 1):
        window = probes.get(start)
        if window is None or window[1] > gap_end:
            continue
        distance = abs((window[0] + window[1]) / 2 - middle)
        if best is None or distance < best[0]:
            best = (distance, window)
    return None if best is None else best[1]


def product_range(ideal_product):
    """(shortest, longest) product allowed for this ideal length."""
    return max(MIN_PRODUCT, ideal_product - PRODUCT_SPREAD), ideal_product + PRODUCT_SPREAD


def possible_designs(genome, region_start, region_end, ideal_product=IDEAL_PRODUCT):
    """Every forward/reverse/probe design inside the target, best first."""
    shortest, longest = product_range(ideal_product)
    # forward and probe: anchored at their 5' end (left), grow right
    forwards = all_windows(genome, range(region_start, region_end), PRIMER_LENGTH,
                           region_start, region_end, "right")
    probes = all_windows(genome, range(region_start, region_end), PROBE_LENGTH,
                         region_start, region_end, "right")
    # reverse: anchored at its 5' end, which is the window's RIGHT end; grows left
    reverses = all_windows(genome, range(region_start + 1, region_end + 1), PRIMER_LENGTH,
                           region_start, region_end, "left")

    designs = []
    for forward_start, (_, forward_end) in forwards.items():
        for product_end in range(forward_start + shortest, forward_start + longest + 1):
            if product_end not in reverses:
                continue
            reverse_start = reverses[product_end][0]
            if reverse_start - forward_end < PROBE_LENGTH[0]:
                continue  # no room for the probe between the primers
            probe_window = centred_probe(probes, forward_end, reverse_start)
            if probe_window is None:
                continue
            forward = genome[forward_start:forward_end]
            reverse_window = genome[reverse_start:product_end]
            probe = genome[probe_window[0]:probe_window[1]]
            designs.append({
                "forward": {"seq": forward, "start": forward_start, "end": forward_end},
                "reverse": {"seq": reverse_complement(reverse_window),
                            "start": reverse_start, "end": product_end},
                "probe": {"seq": probe, "start": probe_window[0], "end": probe_window[1]},
                "product_start": forward_start,
                "product_end": product_end,
                "product_length": product_end - forward_start,
                "match_penalty": match_penalty(forward, reverse_window, probe,
                                               product_end - forward_start, ideal_product),
            })
    designs.sort(key=lambda d: (d["match_penalty"], d["product_start"]))
    return designs


# =============================================================================
# SECTION 3b: RANKING: BEST MATCH FIRST
# -----------------------------------------------------------------------------
# Every design is copied from the target, so every oligo matches it
# perfectly. "Best match" means closest to the ideal design:
#
#   penalty = |forward GC - 50| + |reverse GC - 50| + |probe GC - 50|
#             (in percentage points)
#           + |product length - ideal product length| / 10
#
# Lower is better. 1 point = 1 GC percentage point on one oligo = 10 bp of
# product length. Ties go to the earlier start.
# =============================================================================

def match_penalty(forward, reverse_window, probe, product_length, ideal_product=IDEAL_PRODUCT):
    gc_points = sum(abs(100 * gc_fraction(seq) - 100 * IDEAL_GC)
                    for seq in (forward, reverse_window, probe))
    return round(gc_points + abs(product_length - ideal_product) / 10, 2)


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
                       max_mismatches=max_mismatches, k=K,
                       max_length=max(1000, design["product_length"]))
    products = result["products"]
    return (len(products) == 1 and products[0]["status"] == "product"
            and products[0]["start"] == design["product_start"]
            and products[0]["end"] == design["product_end"])


# =============================================================================
# SECTION 5: THE WHOLE DESIGN IN ONE CALL
# =============================================================================

def design_primers(genome, target, max_mismatches=2, ideal_product=IDEAL_PRODUCT,
                   designs_wanted=DESIGNS_WANTED):
    """Design up to designs_wanted primer pairs + probes for target, best first."""
    genome, target = clean_sequence(genome), clean_sequence(target)
    if not 0 <= max_mismatches <= MAX_ALLOWED_MISMATCHES:
        raise ValueError(f"mismatch allowance must be 0 to {MAX_ALLOWED_MISMATCHES}")
    if ideal_product < MIN_PRODUCT:
        raise ValueError(f"the product length must be at least {MIN_PRODUCT} bp "
                         f"(room for an 18-base forward primer, probe and reverse primer)")
    if len(target) < ideal_product:
        raise ValueError(f"the target ({len(target)} bases) is shorter than the "
                         f"{ideal_product} bp product asked for")

    # -- the GC check comes first: a too GC-rich target stops here
    warnings = []
    target_warning = check_target_gc(target)
    if target_warning:
        warnings.append(target_warning)
    if not QPCR_PRODUCT[0] <= ideal_product <= QPCR_PRODUCT[1]:
        warnings.append(f"a {ideal_product} bp product is outside {QPCR_PRODUCT[0]}-"
                        f"{QPCR_PRODUCT[1]} bp, where qPCR works best")

    region_start, region_end, strand = find_target(genome, target)
    designs = possible_designs(genome, region_start, region_end, ideal_product)

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

    if not designs:
        warnings.append("no design: the target has no forward, reverse and probe "
                        f"stretches with {round(100 * OLIGO_GC[0])}-{round(100 * OLIGO_GC[1])} % GC "
                        "that fit together into a {}-{} bp product".format(*product_range(ideal_product)))

    return {
        "warnings": warnings,
        "target": {"start": region_start, "end": region_end, "strand": strand},
        "ideal_product": ideal_product,
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
    shortest, longest = product_range(result["ideal_product"])
    print(f"Product length: aiming for {result['ideal_product']} bp "
          f"(designs from {shortest} to {longest} bp)")
    print(f"{result['possible']} possible designs; {result['checked']} checked against "
          f"the whole genome (mismatch allowance {max_mismatches}).\n")
    for warning in result["warnings"]:
        print(f"WARNING: {warning}")
    if result["warnings"]:
        print()
    if not result["designs"]:
        if result["possible"]:
            print("No design binds only at its own site.")
        return
    for number, design in enumerate(result["designs"], 1):
        forward, reverse, probe = design["forward"], design["reverse"], design["probe"]
        print(f"Design {number}   (match penalty {design['match_penalty']}, lower is better)")
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
        product = int(sys.argv[4]) if len(sys.argv) > 4 else IDEAL_PRODUCT
    else:
        project = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        genome = read_fasta(glob.glob(os.path.join(project, "data", "reference", "*"))[0])[1]
        target = genome[28274 - 1:29533]
        allowance = 2
        product = IDEAL_PRODUCT
        print("Demo: N gene of the reference genome NC_045512.2\n")

    try:
        result = design_primers(genome, target, allowance, product)
    except ValueError as error:
        print(f"ERROR: {error}")
        sys.exit(1)
    print_designs(result, allowance)
