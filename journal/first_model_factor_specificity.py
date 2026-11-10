# =============================================================================
# FIRST MODEL: FACTOR: SPECIFICITY   (version 4)
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
#   - reports each mismatch by its distance from the primer's 3' end, the
#     base change (primer base vs genome base) and how severe that kind of
#     mismatch is (severe or mild)
#   - grades mismatches close to the 3' end as "3' risk: high / low / possible"
#     and marks a product "at risk" when a primer's 3' risk is high
#   - accepts degenerate bases in primers (Y = C or T, R = A or G, N = any, ...)
#     and reads ambiguity codes in the genome (genome Y = C or T)
#   - lets a site with a 1-base insertion/deletion in the primer's 5' half
#     still bind (it counts as 1 mismatch)
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
#   - melting temperature / annealing temperature: see second_model_annealing_temp.py
#   - primer dimers
#
# Files that go with this model (all in journal/):
#   test_first_model_factor_specificity.py - fake, hand-solvable test cases
#                                            plus search-completeness tests
#   cdc_reference_check.py                 - CDC N1/N2/N3 (+ RNase P control)
#                                            on the real reference genome
#   test_cases_specificity.txt             - every test case and its outcome
#   modification documentation.md          - every problem found and fixed
#
# Run the tests with:     python3 journal/test_first_model_factor_specificity.py
# Run the CDC check with: python3 journal/cdc_reference_check.py
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
#
# Ambiguity codes in the GENOME are read for what they say, not thrown away:
#   every base the genome letter could be is accepted by the primer -> match
#   none of them is accepted                                        -> mismatch
#   some are, some are not                                          -> unknown
# e.g. primer C vs genome Y (C or T) is unknown, primer C vs genome R (A or G)
# is a certain mismatch, primer Y vs genome Y is a certain match.
# Anything that is not a DNA letter at all (e.g. "-") is unknown.
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


def genome_vs_primer(genome_base, primer_base):
    """Return 'match', 'mismatch' or 'unknown' for one genome letter facing one primer letter."""
    could_be = IUPAC.get(genome_base)
    if could_be is None:
        return "unknown"   # not a DNA letter at all
    accepted = IUPAC[primer_base]
    if all(base in accepted for base in could_be):
        return "match"
    if not any(base in accepted for base in could_be):
        return "mismatch"
    return "unknown"


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
#
# k = 6 (the default in run_model). There are 4^6 = 4,096 possible 6-mers,
# so in a 30,000-base genome a random seed hits about 7 places (k = 5: about
# 29), leaving fewer wrong candidates to check. The price: every primer must
# be at least k x (max_mismatches + 1) bases long, 18 bases with the default
# 2 mismatches, which is the short end of real qPCR primers (18-30 bases).
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
#   mismatch - the genome letter certainly differs from what the primer accepts
#   unknown  - the genome has N/an ambiguity code that might or might not
#              match, or the position lies beyond the start or end of the
#              sequence (no data at all)
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
        if position < 0 or position >= len(target):
            unknowns.append(offset)
            continue
        verdict = genome_vs_primer(target[position], primer_base)
        if verdict == "mismatch":
            mismatches.append(offset)
        elif verdict == "unknown":
            unknowns.append(offset)
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
            "severity": mismatch_severity(primer_base, genome_base),
        })
    return sorted(changes, key=lambda change: change["from_3prime"])


# =============================================================================
# SECTION 3b: HOW SEVERE IS A MISMATCH? (graded, not all-or-nothing)
# -----------------------------------------------------------------------------
# Not every mismatch hurts PCR equally. In a mismatch, the primer base sits
# opposite a TEMPLATE base (the complement of the genome base shown in the
# report). Measured effects near the 3' end (Stadhouders et al. 2010,
# J Mol Diagn 12:109-117; Kwok et al. 1990, Nucleic Acids Res 18:999-1005):
#   severe - primer-template A-A, G-A, A-G and C-C: these block extension
#            strongly (large rise in Ct, or no product)
#   mild   - every other pairing, e.g. G-T, T-G, T-T, T-C: the polymerase
#            often extends through them with only a small delay
# A degenerate letter counts as severe if ANY base it stands for is severe.
#
# Severity changes the 3' risk grade below; the mismatch allowance still
# counts every mismatch as 1, so the search guarantees in Section 4 hold.
# =============================================================================

SEVERE_PAIRS = {("A", "A"), ("G", "A"), ("A", "G"), ("C", "C")}   # (primer, template)


def mismatch_severity(primer_base, genome_base):
    """'severe' or 'mild' for a primer base facing the template opposite genome_base."""
    template = COMPLEMENT[genome_base]
    for primer_option in IUPAC[primer_base]:
        for template_option in IUPAC[template]:
            if (primer_option, template_option) in SEVERE_PAIRS:
                return "severe"
    return "mild"


# *****************************************************************************
# *****************************************************************************
# ***                                                                       ***
# ***   >>> USER SETTING: 3' END RISK WINDOW <<<                            ***
# ***                                                                       ***
# ***   Mismatches this many bases (or fewer) from the 3' end give the      ***
# ***   site and its product a 3' risk grade (default: 3):                  ***
# ***     high     - a severe mismatch in the window                        ***
# ***     low      - only mild mismatches in the window                     ***
# ***     possible - no mismatch, but an unknown base in the window         ***
# ***   The site still counts as binding; a HIGH grade turns its product    ***
# ***   "at risk" (ranked below a clean "product").                         ***
# ***                                                                       ***
# *****************************************************************************
# *****************************************************************************
THREE_PRIME_WINDOW = 3
RISK_ORDER = [None, "possible", "low", "high"]   # least to most worrying


def three_prime_risk(changes, unknown_from_3prime):
    """Grade the 3' end of one site: 'high', 'low', 'possible' or None."""
    near = [c for c in changes if c["from_3prime"] <= THREE_PRIME_WINDOW]
    if any(c["severity"] == "severe" for c in near):
        return "high"
    if near:
        return "low"
    if any(d <= THREE_PRIME_WINDOW for d in unknown_from_3prime):
        return "possible"
    return None


def worst_risk(*risks):
    """The most worrying of several 3' risk grades."""
    return max(risks, key=RISK_ORDER.index)


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
            # If the site would still fit the allowance even if EVERY unknown
            # base were a mismatch, the missing data cannot change the answer.
            if not unknowns:
                outcome, reason = "binds", ""
            elif len(mismatches) + len(unknowns) <= max_mismatches:
                outcome, reason = "binds", "binds even if every unknown base is a mismatch"
            elif start < 0:
                outcome, reason = "not assessable", "runs off the start of the sequence"
            elif start + len(pattern) > len(target):
                outcome, reason = "not assessable", "runs off the end of the sequence"
            else:
                outcome, reason = "not assessable", "overlaps N/ambiguous bases"

            changes = describe_changes(target, pattern, [(o, start + o) for o in mismatches], strand)
            unknown_from_3prime = distances_from_3prime(unknowns, len(pattern), strand)
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
                "indel": None,
                "three_prime_risk": three_prime_risk(changes, unknown_from_3prime),
                "unknown_from_3prime": unknown_from_3prime,
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
# So an 18-base primer can show deletions of up to 4 bases, while 20-30 base
# primers (like the CDC ones) can show a 6-base deletion like Alpha's 69-70
# deletion.
#
# Limit: the search starts from exact k-base (6-base) stretches the primer
# shares with the genome. If the differences cut the primer into pieces all
# shorter than 6 bases, the near-miss is not found ("none close"). An 18- or
# 20-base primer is always found with up to 2 events, a 24-base primer with
# up to 3, a 30-base primer with up to 4.
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


def placements(anchor, length, max_indel=MAX_INDEL):
    """Every (left start, gap position, gap length, kind) worth trying around one seed anchor.

    A seed may sit on either side of the gap, so the left part is tried both at
    the anchor and shifted by the gap length.
    """
    yield anchor, 0, 0, "none"
    for gap_length in range(1, max_indel + 1):
        for gap_at in range(1, length):
            for left_start in {anchor, anchor - gap_length}:
                yield left_start, gap_at, gap_length, "insertion"
            if gap_at + gap_length < length:
                for left_start in {anchor, anchor + gap_length}:
                    yield left_start, gap_at, gap_length, "deletion"


def describe_indel(kind, gap_at, gap_length, length, strand):
    """Describe a gap by its kind, length and distance from the primer's 3' end."""
    if kind == "deletion":
        span = distances_from_3prime(range(gap_at, gap_at + gap_length), length, strand)
    else:
        span = distances_from_3prime([gap_at - 1, gap_at], length, strand)
    return {"kind": kind, "length": gap_length, "from_3prime": (span[0], span[-1])}


def anchors_for(pattern, index, k):
    """Start positions suggested by every exact k-mer the pattern shares with the target."""
    anchors = set()
    for offset in range(len(pattern) - k + 1):
        for hit in seed_hits(index, pattern[offset:offset + k]):
            anchors.add(hit - offset)
    return anchors


def find_near_miss(primer_name, primer, target, index, k):
    """Return the closest placement of a failed primer, or None if nothing is close."""
    best, best_key = None, None
    for strand, pattern in (("+", primer), ("-", reverse_complement(primer))):
        # -- 5a: anchors from every exact k-mer the primer shares with the target,
        # -- 5b: then score every placement around every anchor
        for anchor in anchors_for(pattern, index, k):
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
                key = (events, gap_length)   # fewer events first, then the smaller gap
                if best_key is not None and key >= best_key:
                    continue

                # -- 5c: describe the insertion/deletion, if any
                indel = None
                if kind != "none":
                    indel = describe_indel(kind, gap_at, gap_length, len(pattern), strand)

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
# SECTION 5b: BINDING THROUGH A SMALL INSERTION/DELETION IN THE 5' HALF
# -----------------------------------------------------------------------------
# A 1-base insertion or deletion near a primer's 5' end usually still lets it
# bind and extend: the 3' end, where polymerase starts, is intact. So a site
# whose only gap is
#   - at most MAX_BINDING_INDEL bases long, and
#   - in the primer's 5' half (more than half the primer length from the 3' end)
# still "binds". The gap counts as 1 mismatch toward the allowance, so with
# 0 mismatches allowed no gapped site binds. A gap in the 3' half never binds;
# the near-miss search (Section 5) reports it instead.
#
# Only fully sequenced placements are used, and a gapped site is skipped where
# the same primer already has an ungapped site on that strand.
# Limit: like the near-miss search, this starts from exact 6-base stretches.
# =============================================================================

MAX_BINDING_INDEL = 1


def find_indel_sites(primer_name, primer, target, index, k, max_mismatches, existing):
    """Return extra 'binds' sites that need one small gap in the primer's 5' half."""
    if max_mismatches < 1:
        return []
    found = []
    for strand, pattern in (("+", primer), ("-", reverse_complement(primer))):
        length = len(pattern)

        # -- 5b-1: collect every allowed gapped placement
        candidates = []
        for anchor in anchors_for(pattern, index, k):
            for left_start, gap_at, gap_length, kind in placements(anchor, length, MAX_BINDING_INDEL):
                if kind == "none":
                    continue
                indel = describe_indel(kind, gap_at, gap_length, length, strand)
                if indel["from_3prime"][0] <= length // 2:
                    continue  # gap is in the 3' half: not allowed to bind
                pairs, end = gapped_pairs(length, left_start, gap_at, gap_length, kind)
                if left_start < 0 or end > len(target):
                    continue
                if any(base not in KNOWN_BASES for base in target[left_start:end]):
                    continue
                mismatched = [(o, p) for o, p in pairs if not base_matches(target[p], pattern[o])]
                events = len(mismatched) + 1
                if events <= max_mismatches:
                    candidates.append((events, left_start, end, indel, mismatched))

        # -- 5b-2: keep the best placement per location, skipping ungapped sites
        taken = [(s["start"], s["end"]) for s in existing
                 if s["primer"] == primer_name and s["strand"] == strand]
        for events, left_start, end, indel, mismatched in sorted(candidates, key=lambda c: (c[0], c[1])):
            if any(left_start < t_end and t_start < end for t_start, t_end in taken):
                continue
            taken.append((left_start, end))
            changes = describe_changes(target, pattern, mismatched, strand)
            found.append({
                "primer": primer_name,
                "strand": strand,
                "start": left_start,
                "end": end,
                "outcome": "binds",
                "reason": "",
                "mismatches": events,              # the gap counts as 1
                "mismatch_from_3prime": [c["from_3prime"] for c in changes],
                "changes": changes,
                "indel": indel,
                "three_prime_risk": three_prime_risk(changes, []),
                "unknown_from_3prime": [],
            })
    return found


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
                    "three_prime_risk": worst_risk(left["three_prime_risk"], right["three_prime_risk"]),
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
#   at risk        - as "product", but a primer has a SEVERE mismatch within
#                    THREE_PRIME_WINDOW bases of its 3' end (3' risk HIGH):
#                    it binds, but extension may fail or be delayed (later Ct)
#   no signal      - both primers bind but the probe does not: DNA is copied,
#                    yet a qPCR machine would see no fluorescence
#   not assessable - missing data decides it
# Only the primers' 3' risk counts: a TaqMan probe is never extended, so a
# mismatch at its 3' end does not block anything.
#
# Per assay, the best product status wins, in the order of STATUS_ORDER:
#   product > at risk > not assessable > no signal > no product (nothing could form)
# The same order can rank primer pairs for a user: "at risk" pairs go below
# every pair that gives a clean product.
# =============================================================================

STATUS_ORDER = ["product", "at risk", "not assessable", "no signal", "no product"]   # best to worst

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
    """Return 'product', 'at risk', 'no signal' or 'not assessable' for one possible product."""
    if probe_given and product["probe"] == "mismatched":
        return "no signal"
    parts = [product["left"]["outcome"], product["right"]["outcome"]]
    if probe_given:
        parts.append(product["probe"])
    if "not assessable" in parts:
        return "not assessable"
    if product["three_prime_risk"] == "high":
        return "at risk"
    return "product"


def assay_outcome(products):
    """Return the assay's overall outcome from its products' statuses."""
    statuses = {p["status"] for p in products}
    for status in STATUS_ORDER:
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
# ***   Limit: (primer length) // (max_mismatches + 1) must be >= k (6),    ***
# ***          e.g. an 18-base primer allows at most 2 mismatches.          ***
# ***                                                                       ***
# *****************************************************************************
# *****************************************************************************
def run_model(forward, reverse, target, probe=None, max_mismatches=2, k=6,
              min_length=20, max_length=1000):
    """Predict every product the primer pair could amplify from target, and explain it."""
    target = target.upper()
    index = build_kmer_index(target, k)
    oligos = {"forward": clean_oligo("forward", forward),
              "reverse": clean_oligo("reverse", reverse)}
    if probe:
        oligos["probe"] = clean_oligo("probe", probe)

    # -- 8a: binding sites for every oligo (plain sites, then small-gap sites)
    sites = []
    for name, seq in oligos.items():
        plain = find_binding_sites(name, seq, target, index, k, max_mismatches)
        sites += plain + find_indel_sites(name, seq, target, index, k, max_mismatches, plain)

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
    warnings = [f"3' risk {p['three_prime_risk']}: {p['left']['primer']} x {p['right']['primer']} "
                f"product at {p['start'] + 1}-{p['end']}"
                for p in products if p["status"] in ("product", "at risk") and p["three_prime_risk"]]

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
    return ", ".join(f"{c['from_3prime']} from 3': primer {c['primer']}, genome {c['genome']} "
                     f"({c['severity']})" for c in changes)


def format_indel(indel):
    low, high = indel["from_3prime"]
    size = "1 base" if indel["length"] == 1 else f"{indel['length']} bases"
    if indel["kind"] == "deletion":
        where = f"primer base {low}" if low == high else f"primer bases {low}-{high}"
        return f"deletion of {size} ({where} from 3' has no genome partner)"
    return f"insertion of {size} (extra genome bases between primer bases {low} and {high} from 3')"


def print_report(title, result):
    """Print sites, products, outcomes, near-misses and warnings."""
    print(f"\n--- {title} ---")
    print(f"sites: {len(result['sites'])}")
    for s in result["sites"]:
        line = (f"  {s['primer']:<8} strand {s['strand']}  {s['start'] + 1:>4}-{s['end']:<4} "
                f"{s['outcome']:<15} mismatches={s['mismatches']}")
        if s["changes"]:
            line += f" ({format_changes(s['changes'])})"
        if s["indel"]:
            line += f" ({format_indel(s['indel'])})"
        if s["three_prime_risk"]:
            line += f"  [3' RISK: {s['three_prime_risk'].upper()}]"
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
# SECTION 10: READING A GENOME FILE (FASTA)
# -----------------------------------------------------------------------------
# A FASTA file is one ">" header line followed by the sequence split over many
# lines. The lines are joined into one string.
# =============================================================================

def read_fasta(path):
    """Return (header, sequence) from a single-record FASTA file."""
    header, lines = "", []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            elif line:
                lines.append(line)
    return header, "".join(lines).upper()
