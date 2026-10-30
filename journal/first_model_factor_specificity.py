# =============================================================================
# FIRST MODEL: FACTOR: SPECIFICITY   (version 3)
# =============================================================================
# 02-120 project: Do early SARS-CoV-2 qPCR primers still match later variants?
#
# This preliminary model answers one question only: given a primer pair (and
# optionally a probe) and a target sequence, which PCR products could the
# pair amplify, and would the probe detect them?
#
# What it does:
#   - indexes the target's k-mers so primers are not slid across every position
#   - finds primer binding sites on BOTH strands
#   - allows a user-chosen number of mismatches
#   - reports each mismatch by its distance from the primer's 3' end AND the
#     base change (primer base vs genome base)
#   - flags mismatches close to the 3' end as "3' at risk"
#   - accepts degenerate primer bases (Y = C or T, R = A or G, N = any, ...)
#   - pairs sites that face each other at a plausible distance into products
#   - checks that the qPCR probe binds inside each product
#   - for every primer that fails, reports its closest near-miss site,
#     including a single insertion or deletion
#   - gives every primer one of THREE outcomes:
#       binds           - a site within the mismatch allowance
#       mismatched      - no such site: too many real A/C/G/T differences
#       not assessable  - the only candidate site touches missing data
#
# Positions are stored 0-based inside the code (Python style) and PRINTED
# 1-based, the way the CDC table and genome browsers count.
#
# What it does NOT do yet (left for later models):
#   - melting temperature / annealing temperature (nearest-neighbor model)
#   - primer dimers
#   - reading any real genome or primer files
#
# Everything is tested on FAKE primers and a FAKE target sequence built by
# hand, so every expected answer can be worked out on paper first
# ("Cases solvable by hand" in the proposal).
# See journal/test_cases_specificity.txt for every case and its outcome.
#
# Run it with:  python3 journal/first_model_factor_specificity.py
# =============================================================================
#
# -----------------------------------------------------------------------------
# HANDLING BAD SEQUENCING DATA
# -----------------------------------------------------------------------------
# The genomes this project uses are finished consensus sequences from NCBI:
# reads were already aligned and assembled by whoever uploaded them, so read
# alignment tools (Salmon, STAR, TopHat) are not needed here. Poor sequencing
# still leaves two traces in a finished genome:
#
#   1. Unknown bases: regions that did not sequence well are written as N, and
#      uncertain bases as IUPAC codes (Y, R, ...).
#   2. Missing ends: the genome starts or stops partway into the reference.
#
# Version 1 counted an N as a mismatch, so missing data could look like a
# mutation and a primer could be wrongly reported as lost. Since version 2 the
# model never counts unknown or missing bases as mismatches. A site whose KNOWN
# bases fit the mismatch allowance, but which overlaps N/ambiguous bases or
# runs off either end of the sequence, is reported as "not assessable" instead
# of "binds" or "mismatched".
#
# GC content is deliberately NOT used to judge data quality: a FASTA genome has
# no quality scores, and SARS-CoV-2 has genuinely AT-rich regions that a GC
# filter would wrongly throw away.
# -----------------------------------------------------------------------------


# =============================================================================
# SECTION 1: DNA HELPERS (complements, degenerate bases, known bases)
# -----------------------------------------------------------------------------
# A primer that binds the minus strand appears on the plus strand as its
# reverse complement, so we need a way to compute it. Ambiguity codes
# (N, Y, R, ...) are complemented too.
#
# Degenerate primers: some published primers and probes contain IUPAC codes
# (the CDC N3 probe has a Y). A primer Y matches a genome C or T.
# In the TARGET, any letter other than A/C/G/T is unknown data (see the
# "HANDLING BAD SEQUENCING DATA" note above).
# =============================================================================

COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}

# which genome bases each primer letter accepts
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

KNOWN_BASES = "ACGT"   # anything else in a target counts as unknown data


def reverse_complement(seq):
    """Return the reverse complement of a DNA string (5'->3')."""
    return "".join(COMPLEMENT[base] for base in reversed(seq.upper()))


def clean_oligo(name, seq):
    """Uppercase a primer/probe and make sure every letter is a DNA/IUPAC code."""
    seq = seq.upper()
    bad = sorted(set(seq) - set(IUPAC))
    if bad:
        raise ValueError(f"{name}: not DNA letters: {''.join(bad)}")
    return seq


def base_matches(genome_base, primer_base):
    """True if the genome base is one the primer letter accepts."""
    return genome_base in IUPAC[primer_base]


def expand_degenerate(seq):
    """Every plain A/C/G/T sequence a degenerate sequence stands for (YA -> CA, TA)."""
    options = [""]
    for base in seq:
        options = [option + choice for option in options for choice in IUPAC[base]]
    return options


# =============================================================================
# SECTION 2: K-MER INDEX OF THE TARGET SEQUENCE
# -----------------------------------------------------------------------------
# Instead of comparing a primer to every position of the target, we record
# where every length-k word (k-mer) occurs. Looking up a primer's k-mers then
# tells us the only places worth checking. A degenerate k-mer is looked up as
# every plain k-mer it stands for.
# =============================================================================

def build_kmer_index(target, k):
    """Map each k-mer in target to the list of positions (0-based) where it starts."""
    index = {}
    for position in range(len(target) - k + 1):
        kmer = target[position:position + k]
        index.setdefault(kmer, []).append(position)
    return index


def seed_hits(index, seed):
    """All target positions where a (possibly degenerate) seed occurs exactly."""
    hits = []
    for plain in expand_degenerate(seed):
        hits += index.get(plain, [])
    return hits


# =============================================================================
# SECTION 3: COMPARING A SITE TO THE PRIMER (mismatches vs unknown bases)
# -----------------------------------------------------------------------------
# Polymerase extends from the primer's 3' end, so a mismatch there matters
# far more than one at the 5' end. Each mismatch is reported as its distance
# from the 3' end: 1 = the very last (3'-terminal) base of the primer.
#
# Positions are split into two kinds:
#   mismatch - the target has a real base (A/C/G/T) the primer does not accept
#   unknown  - the target has N/an ambiguity code, or the position lies
#              beyond the start or end of the sequence (no data at all)
#
# Each mismatch is also described as a base change, both bases written in the
# primer's own 5'->3' direction: "primer C, genome G" means the primer has a C
# where a perfectly matching primer would need a G.
# =============================================================================

def compare_site(target, start, pattern):
    """Return (mismatch offsets, unknown offsets), left to right on the plus strand."""
    mismatches, unknowns = [], []
    for offset, primer_base in enumerate(pattern):
        position = start + offset
        if position < 0 or position >= len(target) or target[position] not in KNOWN_BASES:
            unknowns.append(offset)
        elif not base_matches(target[position], primer_base):
            mismatches.append(offset)
    return mismatches, unknowns


def distances_from_3prime(offsets, length, strand):
    """Convert plus-strand offsets into distances from the primer's 3' end.

    On the plus strand the primer reads left to right, so its 3' end is the
    RIGHT end of the site. On the minus strand the primer reads right to left
    along the plus strand, so its 3' end is the LEFT end of the site.
    """
    if strand == "+":
        return sorted(length - offset for offset in offsets)
    return sorted(offset + 1 for offset in offsets)


def describe_changes(target, pattern, pairs, strand):
    """Describe each mismatched (pattern offset, target position) pair as a base change."""
    changes = []
    for offset, position in pairs:
        primer_base, genome_base = pattern[offset], target[position]
        if strand == "-":   # turn both back into the primer's own direction
            primer_base, genome_base = COMPLEMENT[primer_base], COMPLEMENT[genome_base]
        changes.append({
            "from_3prime": distances_from_3prime([offset], len(pattern), strand)[0],
            "primer": primer_base,
            "genome": genome_base,
        })
    return sorted(changes, key=lambda change: change["from_3prime"])


# *****************************************************************************
# *****************************************************************************
# ***                                                                       ***
# ***   >>> USER SETTING: 3' END RISK WINDOW <<<                            ***
# ***                                                                       ***
# ***   A mismatch this many bases (or fewer) from the 3' end flags the     ***
# ***   site and its product as "3' at risk" (default: 3).                  ***
# ***   The site still counts as binding; the flag is a warning.            ***
# ***                                                                       ***
# *****************************************************************************
# *****************************************************************************
THREE_PRIME_WINDOW = 3


# =============================================================================
# SECTION 4: FINDING PRIMER BINDING SITES ON BOTH STRANDS
# -----------------------------------------------------------------------------
# Seeding rule (pigeonhole principle): if a primer is split into
# (max_mismatches + 1) pieces, at least one piece must match EXACTLY at any
# true site. So we look up the first k bases of each piece in the index, and
# only then check the full primer at those candidate positions.
#
# Seeds cannot land on unknown bases, so two more groups of candidates are
# checked directly (they are few, so this stays fast):
#   - every placement that overlaps an N/ambiguous base
#   - every placement that hangs off the start or the end of the sequence
#
#   "+" strand site: the primer itself appears on the plus strand.
#                    It extends to the RIGHT.
#   "-" strand site: the primer's reverse complement appears on the plus
#                    strand. It extends to the LEFT.
# =============================================================================

# A "not assessable" site needs at least this many KNOWN bases that agree with
# the primer; otherwise a run of Ns would make every primer "not assessable".
MIN_KNOWN_BASES = 10


def candidate_starts(pattern, target, index, k, max_mismatches):
    """Return every start position worth comparing the full pattern against."""
    pieces = max_mismatches + 1
    piece_length = len(pattern) // pieces
    candidates = set()

    # -- 4a: exact seed hits from the k-mer index
    for piece in range(pieces):
        offset = piece * piece_length
        for hit in seed_hits(index, pattern[offset:offset + k]):
            candidates.add(hit - offset)

    # -- 4b: placements overlapping an N/ambiguous base
    for position, base in enumerate(target):
        if base not in KNOWN_BASES:
            candidates.update(range(position - len(pattern) + 1, position + 1))

    # -- 4c: placements hanging off either end of the sequence
    candidates.update(range(-(len(pattern) - MIN_KNOWN_BASES), 0))
    candidates.update(range(len(target) - len(pattern) + 1,
                            len(target) - MIN_KNOWN_BASES + 1))
    return candidates


def find_binding_sites(primer_name, primer, target, index, k, max_mismatches):
    """Return every site that binds or is not assessable for this primer."""
    if len(primer) // (max_mismatches + 1) < k:
        raise ValueError(
            f"{primer_name}: {len(primer)} bases is too short to seed with k={k} "
            f"while allowing {max_mismatches} mismatches; lower k or the tolerance."
        )

    sites = []
    for strand, pattern in (("+", primer), ("-", reverse_complement(primer))):
        for start in sorted(candidate_starts(pattern, target, index, k, max_mismatches)):
            mismatches, unknowns = compare_site(target, start, pattern)
            known = len(pattern) - len(unknowns)
            if len(mismatches) > max_mismatches or known < MIN_KNOWN_BASES:
                continue  # too different (or too little data) to be a site

            # -- 4d: decide the outcome for this site
            if not unknowns:
                outcome, reason = "binds", ""
            elif start < 0:
                outcome, reason = "not assessable", "runs off the start of the sequence"
            elif start + len(pattern) > len(target):
                outcome, reason = "not assessable", "runs off the end of the sequence"
            else:
                outcome, reason = "not assessable", "overlaps N/ambiguous bases"

            changes = describe_changes(target, pattern, [(o, start + o) for o in mismatches], strand)
            sites.append({
                "primer": primer_name,
                "strand": strand,
                "start": start,                    # 0-based, plus-strand coords
                "end": start + len(pattern),       # exclusive
                "outcome": outcome,
                "reason": reason,
                "mismatches": len(mismatches),
                "mismatch_from_3prime": [change["from_3prime"] for change in changes],
                "changes": changes,
                "three_prime_risk": any(c["from_3prime"] <= THREE_PRIME_WINDOW for c in changes),
                "unknown_from_3prime": distances_from_3prime(unknowns, len(pattern), strand),
            })
    return sites


# =============================================================================
# SECTION 5: NEAR-MISS SEARCH FOR PRIMERS THAT FAIL (incl. insertions/deletions)
# -----------------------------------------------------------------------------
# When a primer is "mismatched", the useful scientific answer is WHY: where
# did it almost bind, and what changed there? This section finds the closest
# site, allowing at most ONE insertion or deletion (up to MAX_INDEL bases)
# plus substitutions. Differences are counted as "events":
#   each substitution = 1 event, one insertion/deletion of any length = 1 event
# so a 6-base deletion (like Alpha's spike 69-70 deletion) is ONE event,
# not 6 mismatches.
#
#   deletion  - the genome is missing bases the primer has
#   insertion - the genome has extra bases the primer does not have
#
# Only fully sequenced placements (no N, not off the ends) are considered.
# Two rules stop coincidences in a long genome from being reported:
#   - at most NEAR_MISS_MAX_EVENTS events, or NEAR_MISS_MAX_EVENTS_WITH_INDEL
#     if one of them is an insertion/deletion (a gap can sit in many places,
#     which gives coincidences many more chances)
#   - at least MIN_MATCHED_BASES primer bases must match the genome exactly
#     (a deletion leaves fewer bases to compare, which makes chance matches
#     easy: on a random 30,000-base sequence a 20-base primer found a fake
#     "6-base deletion + 1 mismatch" before this rule was added)
# So a 16-base primer can show deletions of up to 2 bases, while the 20-24
# base CDC primers can show a 6-base deletion like Alpha's 69-70 deletion.
#
# Limit: the search starts from exact k-base (5-base) stretches the primer
# shares with the genome. If the differences cut the primer into pieces all
# shorter than 5 bases, the near-miss is not found ("none close"). A 16-base
# primer is always found with up to 2 events, a 20-base primer with up to 3.
# =============================================================================

MAX_INDEL = 6
NEAR_MISS_MAX_EVENTS = 4
NEAR_MISS_MAX_EVENTS_WITH_INDEL = 2   # the indel + at most 1 substitution
MIN_MATCHED_BASES = 14


def gapped_pairs(length, left_start, gap_at, gap_length, kind):
    """Pair each primer offset with the genome position it faces; also return the site end.

    kind "none":      no gap, primer offset o faces left_start + o
    kind "deletion":  primer offsets gap_at..gap_at+gap_length-1 face nothing
    kind "insertion": gap_length extra genome bases sit before primer offset gap_at
    """
    if kind == "none":
        return [(o, left_start + o) for o in range(length)], left_start + length
    left = [(o, left_start + o) for o in range(gap_at)]
    if kind == "deletion":
        right = [(o, left_start + o - gap_length) for o in range(gap_at + gap_length, length)]
        return left + right, left_start + length - gap_length
    right = [(o, left_start + o + gap_length) for o in range(gap_at, length)]
    return left + right, left_start + length + gap_length


def placements(anchor, length):
    """Every (left start, gap position, gap length, kind) worth trying around one seed anchor.

    A seed may sit on either side of the gap, so the left part is tried both at
    the anchor and shifted by the gap length.
    """
    yield anchor, 0, 0, "none"
    for gap_length in range(1, MAX_INDEL + 1):
        for gap_at in range(1, length):
            for left_start in {anchor, anchor - gap_length}:
                yield left_start, gap_at, gap_length, "insertion"
            if gap_at + gap_length < length:
                for left_start in {anchor, anchor + gap_length}:
                    yield left_start, gap_at, gap_length, "deletion"


def find_near_miss(primer_name, primer, target, index, k):
    """Return the closest placement of a failed primer, or None if nothing is close."""
    best, best_key = None, None
    for strand, pattern in (("+", primer), ("-", reverse_complement(primer))):
        # -- 5a: anchors from every exact k-mer the primer shares with the target
        anchors = set()
        for offset in range(len(pattern) - k + 1):
            for hit in seed_hits(index, pattern[offset:offset + k]):
                anchors.add(hit - offset)

        # -- 5b: score every placement around every anchor
        for anchor in anchors:
            for left_start, gap_at, gap_length, kind in placements(anchor, len(pattern)):
                pairs, end = gapped_pairs(len(pattern), left_start, gap_at, gap_length, kind)
                if left_start < 0 or end > len(target):
                    continue
                if any(base not in KNOWN_BASES for base in target[left_start:end]):
                    continue
                mismatched = [(o, p) for o, p in pairs if not base_matches(target[p], pattern[o])]
                if len(pairs) - len(mismatched) < MIN_MATCHED_BASES:
                    continue  # too few matching bases to be more than a coincidence
                events = len(mismatched) + (kind != "none")
                limit = NEAR_MISS_MAX_EVENTS if kind == "none" else NEAR_MISS_MAX_EVENTS_WITH_INDEL
                if events > limit:
                    continue  # too different to be more than a coincidence
                key =(events, gap_length)   # fewer events first, then the smaller gap
                if best_key is not None and key >= best_key:
                    continue

                # -- 5c: describe the insertion/deletion, if any
                indel = None
                if kind == "deletion":
                    span = distances_from_3prime(range(gap_at, gap_at + gap_length), len(pattern), strand)
                    indel = {"kind": "deletion", "length": gap_length, "from_3prime": (span[0], span[-1])}
                elif kind == "insertion":
                    span = distances_from_3prime([gap_at - 1, gap_at], len(pattern), strand)
                    indel = {"kind": "insertion", "length": gap_length, "from_3prime": (span[0], span[-1])}

                best_key = key
                best = {
                    "primer": primer_name,
                    "strand": strand,
                    "start": left_start,
                    "end": end,
                    "events": events,
                    "changes": describe_changes(target, pattern, mismatched, strand),
                    "indel": indel,
                }
    return best


# =============================================================================
# SECTION 6: PAIRING SITES INTO PCR PRODUCTS
# -----------------------------------------------------------------------------
# A product forms when a "+" site (extending right) sits upstream of a "-"
# site (extending left) and the distance between their outer ends is within
# a plausible product length. The forward and reverse primers can pair with
# each other or with themselves, because that is what happens in a real tube.
# The probe is not a primer and never makes products.
#
# Sites that are "not assessable" are paired too, so the model can tell
# "no product possible" apart from "a product might exist, but data is
# missing". Section 7 decides which is which.
# =============================================================================

def pair_into_products(sites, min_length, max_length):
    """Return every possible product from '+' sites facing downstream '-' sites."""
    usable = [s for s in sites if s["primer"] in ("forward", "reverse")]
    plus_sites = [s for s in usable if s["strand"] == "+"]
    minus_sites = [s for s in usable if s["strand"] == "-"]
    products = []
    for left in plus_sites:
        for right in minus_sites:
            if right["end"] <= left["start"]:
                continue  # the two primers point away from each other
            length = right["end"] - left["start"]
            if min_length <= length <= max_length and right["start"] >= left["start"]:
                products.append({
                    "left": left,
                    "right": right,
                    "start": left["start"],
                    "end": right["end"],
                    "length": length,
                    "three_prime_risk": left["three_prime_risk"] or right["three_prime_risk"],
                })
    products.sort(key=lambda p: (p["start"], p["length"]))
    return products


# =============================================================================
# SECTION 7: OUTCOMES (per primer, per product, per assay)
# -----------------------------------------------------------------------------
# Per primer (forward, reverse, probe):
#   binds          - at least one site binds
#   not assessable - no site binds, but one touches missing data
#   mismatched     - no site within the mismatch allowance at all
#
# Probe inside a product: the best outcome among probe sites lying fully
# inside that product ("mismatched" if there are none).
#
# Per product (status):
#   product        - both primer sites bind (and the probe binds, if given)
#   no signal      - both primers bind but the probe does not: DNA is copied,
#                    yet a qPCR machine would see no fluorescence
#   not assessable - missing data decides it
#
# Per assay, the best product status wins, in this order:
#   product > not assessable > no signal > no product (nothing could form)
# =============================================================================

def oligo_outcome(name, sites):
    """Return 'binds', 'not assessable' or 'mismatched' for one primer or probe."""
    outcomes = {s["outcome"] for s in sites if s["primer"] == name}
    if "binds" in outcomes:
        return "binds"
    if "not assessable" in outcomes:
        return "not assessable"
    return "mismatched"


def probe_inside(product, sites):
    """Return the probe's outcome inside one product."""
    inside = [s for s in sites if s["primer"] == "probe"
              and s["start"] >= product["start"] and s["end"] <= product["end"]]
    return oligo_outcome("probe", inside)


def product_status(product, probe_given):
    """Return 'product', 'no signal' or 'not assessable' for one possible product."""
    if probe_given and product["probe"] == "mismatched":
        return "no signal"
    parts = [product["left"]["outcome"], product["right"]["outcome"]]
    if probe_given:
        parts.append(product["probe"])
    if "not assessable" in parts:
        return "not assessable"
    return "product"


def assay_outcome(products):
    """Return the assay's overall outcome from its products' statuses."""
    statuses = {p["status"] for p in products}
    for status in ("product", "not assessable", "no signal"):
        if status in statuses:
            return status
    return "no product"


# =============================================================================
# SECTION 8: THE WHOLE MODEL IN ONE CALL
# -----------------------------------------------------------------------------
# Primer pair (+ optional probe) and target sequence in; one result out with
# sites, products, outcomes, near-misses and warnings.
# =============================================================================

# *****************************************************************************
# *****************************************************************************
# ***                                                                       ***
# ***   >>> USER SETTING: MISMATCH ALLOWANCE <<<                            ***
# ***                                                                       ***
# ***   max_mismatches = how many mismatches a primer may have at a site    ***
# ***                    and still count as binding (default: 2)            ***
# ***                                                                       ***
# ***   Limit: (primer length) // (max_mismatches + 1) must be >= k (5),    ***
# ***          e.g. a 16-base primer allows at most 2 mismatches.           ***
# ***                                                                       ***
# *****************************************************************************
# *****************************************************************************
def run_model(forward, reverse, target, probe=None, max_mismatches=2, k=5,
              min_length=20, max_length=1000):
    """Predict every product the primer pair could amplify from target, and explain it."""
    target = target.upper()
    index = build_kmer_index(target, k)
    oligos = {"forward": clean_oligo("forward", forward),
              "reverse": clean_oligo("reverse", reverse)}
    if probe:
        oligos["probe"] = clean_oligo("probe", probe)

    # -- 8a: binding sites for every oligo
    sites = []
    for name, seq in oligos.items():
        sites += find_binding_sites(name, seq, target, index, k, max_mismatches)

    # -- 8b: products, probe check and status
    products = pair_into_products(sites, min_length, max_length)
    for product in products:
        if probe:
            product["probe"] = probe_inside(product, sites)
        product["status"] = product_status(product, probe is not None)

    # -- 8c: outcomes and near-misses for anything that failed
    outcomes = {name: oligo_outcome(name, sites) for name in oligos}
    near_misses = {name: find_near_miss(name, seq, target, index, k)
                   for name, seq in oligos.items() if outcomes[name] == "mismatched"}

    # -- 8d: warnings worth reading even when a product forms
    warnings = [f"3' at risk: {p['left']['primer']} x {p['right']['primer']} product "
                f"at {p['start'] + 1}-{p['end']}"
                for p in products if p["status"] == "product" and p["three_prime_risk"]]

    return {
        "sites": sites,
        "products": products,
        "outcomes": outcomes,
        "near_misses": near_misses,
        "assay": assay_outcome(products),
        "warnings": warnings,
    }


# =============================================================================
# SECTION 9: PRINTING A READABLE REPORT (positions printed 1-based)
# =============================================================================

def format_changes(changes):
    return ", ".join(f"{c['from_3prime']} from 3': primer {c['primer']}, genome {c['genome']}"
                     for c in changes)


def format_indel(indel):
    low, high = indel["from_3prime"]
    if indel["kind"] == "deletion":
        return (f"deletion of {indel['length']} bases "
                f"(primer bases {low}-{high} from 3' have no genome partner)")
    return (f"insertion of {indel['length']} bases "
            f"(extra genome bases between primer bases {low} and {high} from 3')")


def print_report(title, result):
    """Print sites, products, outcomes, near-misses and warnings."""
    print(f"\n--- {title} ---")
    print(f"sites: {len(result['sites'])}")
    for s in result["sites"]:
        line = (f"  {s['primer']:<8} strand {s['strand']}  {s['start'] + 1:>4}-{s['end']:<4} "
                f"{s['outcome']:<15} mismatches={s['mismatches']}")
        if s["changes"]:
            line += f" ({format_changes(s['changes'])})"
        if s["three_prime_risk"]:
            line += "  [3' AT RISK]"
        if s["reason"]:
            line += f"  [{s['reason']}; unknown from 3' end={s['unknown_from_3prime']}]"
        print(line)

    print(f"products: {len(result['products'])}")
    for p in result["products"]:
        line = (f"  {p['left']['primer']} x {p['right']['primer']}: "
                f"{p['start'] + 1}-{p['end']}  length {p['length']} bp  status: {p['status']}")
        if "probe" in p:
            line += f"  probe: {p['probe']}"
        print(line)

    for name, miss in result["near_misses"].items():
        if miss is None:
            print(f"near-miss for {name}: none close")
            continue
        details = [format_changes(miss["changes"])] if miss["changes"] else []
        if miss["indel"]:
            details.append(format_indel(miss["indel"]))
        print(f"near-miss for {name}: strand {miss['strand']} {miss['start'] + 1}-{miss['end']}, "
              f"{miss['events']} event(s): {'; '.join(details)}")

    for warning in result["warnings"]:
        print(f"WARNING: {warning}")

    outcomes = ", ".join(f"{name} {outcome}" for name, outcome in result["outcomes"].items())
    print(f"outcome: {outcomes}, assay {result['assay']}")


# =============================================================================
# SECTION 10: FAKE PRIMERS AND FAKE TARGET BUILDING BLOCKS
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
# SECTION 11: CASES SOLVABLE BY HAND: CORE (version 1)
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
def case_mismatches():
    seq = mutate(mutate(amplicon(), 15, "G"), 36, "T")
    target = filler(10) + seq + filler(10)

    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 3a: mismatches, tolerance 0", result)
    check(product_lengths(result) == [], "no product when no mismatches are allowed")
    check(result["outcomes"]["forward"] == "mismatched", "forward is 'mismatched'")
    check(result["assay"] == "no product", "assay is 'no product'")
    check(result["near_misses"]["forward"]["changes"] ==
          [{"from_3prime": 1, "primer": "C", "genome": "G"}],
          "forward near-miss: 1 from 3' end, primer C, genome G")
    check(result["near_misses"]["reverse"]["changes"] ==
          [{"from_3prime": 1, "primer": "G", "genome": "A"}],
          "reverse near-miss: 1 from 3' end, primer G, genome A")

    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 3b: mismatches, tolerance 1", result)
    check(product_lengths(result) == [52], "one product when 1 mismatch is allowed")
    check(result["products"][0]["left"]["mismatch_from_3prime"] == [1],
          "forward mismatch reported at the 3'-terminal base")
    check(result["products"][0]["right"]["mismatch_from_3prime"] == [1],
          "reverse mismatch reported at the 3'-terminal base")
    check(result["products"][0]["three_prime_risk"], "product is flagged 3' at risk")
    check(len(result["warnings"]) == 1, "one 3' at risk warning is printed")


# -- Case 3c: mismatch far from the 3' end -----------------------------------
# Change the forward site's FIRST base (offset 0): distance 16 from the 3' end.
def case_mismatch_5prime():
    target = filler(10) + mutate(amplicon(), 0, "C") + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 3c: mismatch at the 5' end", result)
    check(result["products"][0]["left"]["changes"] ==
          [{"from_3prime": 16, "primer": "G", "genome": "C"}],
          "forward mismatch 16 from 3' end: primer G, genome C")
    check(not result["products"][0]["three_prime_risk"], "a 5' mismatch is NOT 3' at risk")


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
# SECTION 12: CASES SOLVABLE BY HAND: BAD SEQUENCING DATA (version 2)
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
def case_mismatched_beats_unknown():
    seq = amplicon()
    for offset in (40, 44, 48):
        seq = change_base(seq, offset)
    seq = mutate(seq, 3, "N")
    target = filler(10) + seq + filler(10)
    result = run_model(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 10: reverse truly mismatched, N elsewhere", result)
    check(result["outcomes"]["reverse"] == "mismatched", "reverse is 'mismatched'")
    check(result["outcomes"]["forward"] == "not assessable", "forward is 'not assessable'")
    check(result["assay"] == "no product", "assay is 'no product': the N cannot rescue it")
    check(result["near_misses"]["reverse"] is None,
          "near-miss not found: no exact 5-base stretch left (known limit)")


# =============================================================================
# SECTION 13: CASES SOLVABLE BY HAND: EDGE CASES (version 3)
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
    check(site_of(result, "forward")["changes"] == [{"from_3prime": 8, "primer": "Y", "genome": "A"}],
          "Y does not match genome A: 8 from 3' end, primer Y, genome A")


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
# SECTION 14: RUN EVERY CASE
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
    print("\nAll hand-solvable cases passed.")
