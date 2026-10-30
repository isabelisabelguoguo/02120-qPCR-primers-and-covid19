# =============================================================================
# FIRST MODEL: FACTOR: SPECIFICITY   (version 2)
# =============================================================================
# 02-120 project: Do early SARS-CoV-2 qPCR primers still match later variants?
#
# This preliminary model answers one question only: given a primer pair and a
# target sequence, which PCR products could the pair amplify?
#
# What it does:
#   - indexes the target's k-mers so primers are not slid across every position
#   - finds primer binding sites on BOTH strands
#   - allows a user-chosen number of mismatches
#   - reports each mismatch by its distance from the primer's 3' end
#   - pairs sites that face each other at a plausible distance into products
#   - gives every primer one of THREE outcomes (new in version 2):
#       binds           - a site within the mismatch allowance
#       mismatched      - no such site: too many real A/C/G/T differences
#       not assessable  - the only candidate site touches missing data
#
# What it does NOT do yet (left for later models):
#   - melting temperature / annealing temperature (nearest-neighbor model)
#   - primer dimers
#   - insertions/deletions (a deletion shows up as several mismatches)
#   - reading any real genome or primer files
#
# Everything is tested on FAKE primers and a FAKE target sequence built by
# hand, so every expected answer can be worked out on paper first
# ("Cases solvable by hand" in the proposal).
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
# mutation and a primer could be wrongly reported as lost. Version 2 never
# counts unknown or missing bases as mismatches. A site whose KNOWN bases fit
# the mismatch allowance, but which overlaps N/ambiguous bases or runs off
# either end of the sequence, is reported as "not assessable" instead of
# "binds" or "mismatched".
#
# GC content is deliberately NOT used to judge data quality: a FASTA genome has
# no quality scores, and SARS-CoV-2 has genuinely AT-rich regions that a GC
# filter would wrongly throw away.
# -----------------------------------------------------------------------------


# =============================================================================
# SECTION 1: DNA HELPERS (reverse complement, known vs unknown bases)
# -----------------------------------------------------------------------------
# A primer that binds the minus strand appears on the plus strand as its
# reverse complement, so we need a way to compute it. Ambiguity codes
# (N, Y, R, ...) are complemented too, so real genomes containing them work.
# =============================================================================

COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}

KNOWN_BASES = "ACGT"   # anything else in a target counts as unknown data


def reverse_complement(seq):
    """Return the reverse complement of a DNA string (5'->3')."""
    return "".join(COMPLEMENT[base] for base in reversed(seq.upper()))


# =============================================================================
# SECTION 2: K-MER INDEX OF THE TARGET SEQUENCE
# -----------------------------------------------------------------------------
# Instead of comparing a primer to every position of the target, we record
# where every length-k word (k-mer) occurs. Looking up a primer's k-mers then
# tells us the only places worth checking.
# =============================================================================

def build_kmer_index(target, k):
    """Map each k-mer in target to the list of positions (0-based) where it starts."""
    index = {}
    for position in range(len(target) - k + 1):
        kmer = target[position:position + k]
        index.setdefault(kmer, []).append(position)
    return index


# =============================================================================
# SECTION 3: COMPARING A SITE TO THE PRIMER (mismatches vs unknown bases)
# -----------------------------------------------------------------------------
# Polymerase extends from the primer's 3' end, so a mismatch there matters
# far more than one at the 5' end. Each mismatch is reported as its distance
# from the 3' end: 1 = the very last (3'-terminal) base of the primer.
#
# Positions are split into two kinds:
#   mismatch - the target has a real base (A/C/G/T) that differs
#   unknown  - the target has N/an ambiguity code, or the position lies
#              beyond the start or end of the sequence (no data at all)
# =============================================================================

def compare_site(target, start, pattern):
    """Return (mismatch offsets, unknown offsets), left to right on the plus strand."""
    mismatches, unknowns = [], []
    for offset, primer_base in enumerate(pattern):
        position = start + offset
        if position < 0 or position >= len(target) or target[position] not in KNOWN_BASES:
            unknowns.append(offset)
        elif target[position] != primer_base:
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
        for hit in index.get(pattern[offset:offset + k], []):
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
    primer = primer.upper()
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

            sites.append({
                "primer": primer_name,
                "strand": strand,
                "start": start,                    # 0-based, plus-strand coords
                "end": start + len(pattern),       # exclusive
                "outcome": outcome,
                "reason": reason,
                "mismatches": len(mismatches),
                "mismatch_from_3prime": distances_from_3prime(mismatches, len(pattern), strand),
                "unknown_from_3prime": distances_from_3prime(unknowns, len(pattern), strand),
            })
    return sites


# =============================================================================
# SECTION 5: PAIRING SITES INTO PCR PRODUCTS
# -----------------------------------------------------------------------------
# A product forms when a "+" site (extending right) sits upstream of a "-"
# site (extending left) and the distance between their outer ends is within
# a plausible product length. Any primer can pair with any other primer, or
# with itself, because that is what happens in a real tube.
# Only sites whose outcome is "binds" are used: a product is never predicted
# from missing data.
# =============================================================================

def pair_into_products(sites, min_length, max_length):
    """Return every product formed by a '+' site facing a downstream '-' site."""
    binding = [s for s in sites if s["outcome"] == "binds"]
    plus_sites = [s for s in binding if s["strand"] == "+"]
    minus_sites = [s for s in binding if s["strand"] == "-"]
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
                })
    products.sort(key=lambda p: (p["start"], p["length"]))
    return products


# =============================================================================
# SECTION 6: THE THREE OUTCOMES (per primer and per assay)
# -----------------------------------------------------------------------------
# Per primer:
#   binds          - at least one site binds
#   not assessable - no site binds, but one touches missing data
#   mismatched     - no site within the mismatch allowance at all
# Per assay (the primer pair together):
#   product        - at least one product is predicted
#   not assessable - no product, and a primer is not assessable
#   no product     - no product, and the data says so clearly
# =============================================================================

def primer_outcome(primer_name, sites):
    """Return 'binds', 'not assessable' or 'mismatched' for one primer."""
    outcomes = {s["outcome"] for s in sites if s["primer"] == primer_name}
    if "binds" in outcomes:
        return "binds"
    if "not assessable" in outcomes:
        return "not assessable"
    return "mismatched"


def assay_outcome(sites, products):
    """Return 'product', 'not assessable' or 'no product' for the primer pair."""
    if products:
        return "product"
    if "not assessable" in (primer_outcome("forward", sites), primer_outcome("reverse", sites)):
        return "not assessable"
    return "no product"


# =============================================================================
# SECTION 7: THE WHOLE MODEL IN ONE CALL
# -----------------------------------------------------------------------------
# Primer pair + target sequence in, predicted sites and products out.
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
def predict_products(forward, reverse, target, max_mismatches=2, k=5,
                     min_length=20, max_length=1000):
    """Predict every PCR product the primer pair could amplify from target."""
    target = target.upper()
    index = build_kmer_index(target, k)
    sites = []
    sites += find_binding_sites("forward", forward, target, index, k, max_mismatches)
    sites += find_binding_sites("reverse", reverse, target, index, k, max_mismatches)
    return sites, pair_into_products(sites, min_length, max_length)


def print_report(title, sites, products):
    """Print the binding sites, products and outcomes in a readable form."""
    print(f"\n--- {title} ---")
    print(f"sites: {len(sites)}")
    for s in sites:
        line = (f"  {s['primer']:<8} strand {s['strand']}  {s['start']:>4}-{s['end']:<4} "
                f"{s['outcome']:<15} mismatches={s['mismatches']} "
                f"from 3' end={s['mismatch_from_3prime']}")
        if s["reason"]:
            line += f"  [{s['reason']}; unknown from 3' end={s['unknown_from_3prime']}]"
        print(line)
    print(f"products: {len(products)}")
    for p in products:
        print(f"  {p['left']['primer']} x {p['right']['primer']}: "
              f"{p['start']}-{p['end']}  length {p['length']} bp")
    print(f"outcome: forward {primer_outcome('forward', sites)}, "
          f"reverse {primer_outcome('reverse', sites)}, "
          f"assay {assay_outcome(sites, products)}")


# =============================================================================
# SECTION 8: FAKE PRIMERS AND FAKE TARGET BUILDING BLOCKS
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


def filler(n):
    """n bases of background sequence that no primer can bind."""
    return "A" * n


def amplicon():
    """forward site + insert + reverse site = 16 + 20 + 16 = 52 bp product."""
    return FAKE_FORWARD + INSERT + REVERSE_SITE


def mutate(seq, offset, new_base):
    """Return seq with the base at offset replaced (used to plant mismatches or Ns)."""
    return seq[:offset] + new_base + seq[offset + 1:]


# =============================================================================
# SECTION 9: CASES SOLVABLE BY HAND
# -----------------------------------------------------------------------------
# Each case states its expected answer (worked out on paper), runs the model,
# and checks the model agrees. A failed check stops the program.
# =============================================================================

def check(condition, message):
    if not condition:
        raise AssertionError(message)
    print(f"  PASS: {message}")


# -- Case 1: exact match, one product ----------------------------------------
# target = 10 filler + 52 bp amplicon + 10 filler
# forward binds "+" at 10..26, reverse binds "-" at 46..62, product 10..62 = 52 bp
def case_exact_match():
    target = filler(10) + amplicon() + filler(10)
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 1: exact match", sites, products)
    check(len(products) == 1, "exactly one product")
    check(products[0]["length"] == 52, "product is 52 bp (16 + 20 + 16)")
    check((products[0]["start"], products[0]["end"]) == (10, 62), "product spans 10-62")


# -- Case 2: reverse-strand matches ------------------------------------------
# Reverse-complement the whole Case 1 target. The product is the same 52 bp,
# but now the REVERSE primer is the "+" site and the FORWARD primer is "-".
# target length 72, so the product spans 72-62=10 .. 72-10=62.
def case_reverse_strand():
    target = reverse_complement(filler(10) + amplicon() + filler(10))
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 2: whole target on the other strand", sites, products)
    check(len(products) == 1, "exactly one product")
    check(products[0]["length"] == 52, "product is still 52 bp")
    check(products[0]["left"]["primer"] == "reverse", "reverse primer is now the '+' site")
    check(products[0]["right"]["primer"] == "forward", "forward primer is now the '-' site")


# -- Case 3: introduced mismatches -------------------------------------------
# Plant one mismatch in the forward site at its LAST base (3'-terminal),
# i.e. offset 15 of the amplicon -> distance 1 from the 3' end.
# Plant one mismatch in the reverse site at its FIRST plus-strand base,
# i.e. amplicon offset 36. The reverse primer's 3' end is the LEFT end of a
# "-" site, so this is also distance 1 from its 3' end.
def case_mismatches():
    seq = amplicon()
    seq = mutate(seq, 15, "G")   # forward site last base C -> G
    seq = mutate(seq, 36, "T")   # reverse site first base C -> T
    target = filler(10) + seq + filler(10)

    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 3a: mismatches, tolerance 0", sites, products)
    check(len(products) == 0, "no product when no mismatches are allowed")
    check(primer_outcome("forward", sites) == "mismatched", "forward is 'mismatched'")
    check(assay_outcome(sites, products) == "no product", "assay is 'no product'")

    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 3b: mismatches, tolerance 1", sites, products)
    check(len(products) == 1, "one product when 1 mismatch is allowed")
    check(products[0]["left"]["mismatch_from_3prime"] == [1],
          "forward mismatch reported at the 3'-terminal base")
    check(products[0]["right"]["mismatch_from_3prime"] == [1],
          "reverse mismatch reported at the 3'-terminal base")


# -- Case 3c: mismatch far from the 3' end -----------------------------------
# Change the forward site's FIRST base (offset 0): distance 16 from the 3' end.
def case_mismatch_5prime():
    target = filler(10) + mutate(amplicon(), 0, "C") + filler(10)
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 3c: mismatch at the 5' end", sites, products)
    check(products[0]["left"]["mismatch_from_3prime"] == [16],
          "forward mismatch reported 16 bases from the 3' end")


# -- Case 4: multiple products -----------------------------------------------
# Two copies of the amplicon: 10 filler + amplicon + 30 filler + amplicon + 10.
# Copy 1 at 10..62, copy 2 at 92..144. Products with max_length=1000:
#   copy1 F x copy1 R = 52 bp, copy2 F x copy2 R = 52 bp,
#   copy1 F x copy2 R = 144 - 10 = 134 bp  (the long "outer" product)
# copy2 F x copy1 R faces away -> no product.
def case_multiple_products():
    target = filler(10) + amplicon() + filler(30) + amplicon() + filler(10)
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 4: two copies of the target", sites, products)
    check(len(sites) == 4, "4 binding sites (2 forward, 2 reverse)")
    check(sorted(p["length"] for p in products) == [52, 52, 134],
          "products of 52, 52 and 134 bp")

    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target,
                                       max_mismatches=0, max_length=100)
    print_report("Case 4b: same, max product length 100", sites, products)
    check(sorted(p["length"] for p in products) == [52, 52],
          "the 134 bp product is dropped when the length limit is 100")


# -- Case 5: primers facing away from each other -----------------------------
# Put the reverse site FIRST and the forward site after it:
# they point away from each other, so nothing is amplified.
def case_facing_away():
    target = filler(10) + REVERSE_SITE + INSERT + FAKE_FORWARD + filler(10)
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 5: primers facing away", sites, products)
    check(len(sites) == 2, "both primers still bind")
    check(len(products) == 0, "but no product forms")


# -- Case 6: no binding site at all ------------------------------------------
def case_no_site():
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, filler(200), max_mismatches=2)
    print_report("Case 6: target without the primer sites", sites, products)
    check(len(sites) == 0 and len(products) == 0, "no sites and no products")
    check(assay_outcome(sites, products) == "no product", "assay is 'no product'")


# =============================================================================
# SECTION 10: BAD-DATA CASES SOLVABLE BY HAND (new in version 2)
# -----------------------------------------------------------------------------
# Each case imitates a trace that poor sequencing leaves in a finished genome
# and checks the model says "not assessable" instead of "mismatched".
# =============================================================================

# -- Case 7: Ns inside a primer site -----------------------------------------
# Replace amplicon offsets 5 and 6 (inside the forward site) with N.
# On the "+" strand, offset 5 is 16-5 = 11 and offset 6 is 16-6 = 10 bases
# from the 3' end. With tolerance 0 the Ns must NOT count as mismatches:
# forward is "not assessable" (0 mismatches, 2 unknown), reverse still binds,
# no product is predicted, and the assay is "not assessable".
def case_n_in_site():
    seq = mutate(mutate(amplicon(), 5, "N"), 6, "N")
    target = filler(10) + seq + filler(10)
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 7: Ns inside the forward site", sites, products)
    forward_site = [s for s in sites if s["primer"] == "forward"][0]
    check(forward_site["mismatches"] == 0, "Ns are not counted as mismatches")
    check(forward_site["unknown_from_3prime"] == [10, 11], "Ns reported 10 and 11 bases from the 3' end")
    check(primer_outcome("forward", sites) == "not assessable", "forward is 'not assessable'")
    check(primer_outcome("reverse", sites) == "binds", "reverse still binds")
    check(assay_outcome(sites, products) == "not assessable", "assay is 'not assessable'")


# -- Case 8: sequence ends partway through a primer site ---------------------
# target = 10 filler + forward site + insert + first 12 bases of the reverse
# site. The reverse site would span 46..62 but the sequence stops at 58, so
# 12 of its 16 bases are present (>= 10 required) and all of them match.
def case_truncated_end():
    target = filler(10) + FAKE_FORWARD + INSERT + REVERSE_SITE[:12]
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 8: sequence stops inside the reverse site", sites, products)
    reverse_site = [s for s in sites if s["primer"] == "reverse"][0]
    check(reverse_site["reason"] == "runs off the end of the sequence", "reason is 'runs off the end'")
    check(primer_outcome("reverse", sites) == "not assessable", "reverse is 'not assessable'")
    check(assay_outcome(sites, products) == "not assessable", "assay is 'not assessable'")


# -- Case 9: sequence starts partway through a primer site -------------------
# Drop the first 4 bases of the amplicon: 12 of the forward site's 16 bases
# remain at the very start of the sequence, so it "runs off the start".
def case_truncated_start():
    target = amplicon()[4:] + filler(10)
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=0)
    print_report("Case 9: sequence starts inside the forward site", sites, products)
    forward_site = [s for s in sites if s["primer"] == "forward"][0]
    check(forward_site["start"] == -4, "forward site would have started 4 bases before the sequence")
    check(forward_site["reason"] == "runs off the start of the sequence", "reason is 'runs off the start'")
    check(assay_outcome(sites, products) == "not assessable", "assay is 'not assessable'")


# -- Case 10: real mutations next to missing data are still 'mismatched' ----
# Three real mismatches in the reverse site (amplicon offsets 40, 44, 48)
# plus an N in the forward site. Tolerance 1: the reverse primer has more
# real differences than allowed, so it is "mismatched", not "not assessable",
# even though the target also contains missing data elsewhere.
def case_mismatched_beats_unknown():
    seq = amplicon()
    for offset in (40, 44, 48):
        seq = mutate(seq, offset, "A" if seq[offset] != "A" else "C")
    seq = mutate(seq, 3, "N")
    target = filler(10) + seq + filler(10)
    sites, products = predict_products(FAKE_FORWARD, FAKE_REVERSE, target, max_mismatches=1)
    print_report("Case 10: reverse truly mismatched, N elsewhere", sites, products)
    check(primer_outcome("reverse", sites) == "mismatched", "reverse is 'mismatched'")
    check(primer_outcome("forward", sites) == "not assessable", "forward is 'not assessable'")


# =============================================================================
# SECTION 11: RUN EVERY CASE
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
    print("\nAll hand-solvable cases passed.")
