# =============================================================================
# FIRST MODEL: FACTOR: SPECIFICITY
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
#
# What it does NOT do yet (left for later models):
#   - melting temperature / annealing temperature (nearest-neighbor model)
#   - primer dimers
#   - reading any real genome or primer files
#
# Everything is tested on FAKE primers and a FAKE target sequence built by
# hand, so every expected answer can be worked out on paper first
# ("Cases solvable by hand" in the proposal).
#
# Run it with:  python3 journal/first_model_factor_specificity.py
# =============================================================================


# =============================================================================
# SECTION 1: DNA HELPERS (reverse complement)
# -----------------------------------------------------------------------------
# A primer that binds the minus strand appears on the plus strand as its
# reverse complement, so we need a way to compute it.
# =============================================================================

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


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
# SECTION 3: COUNTING MISMATCHES AND THEIR DISTANCE FROM THE 3' END
# -----------------------------------------------------------------------------
# Polymerase extends from the primer's 3' end, so a mismatch there matters
# far more than one at the 5' end. Each mismatch is reported as its distance
# from the 3' end: 1 = the very last (3'-terminal) base of the primer.
# =============================================================================

def mismatch_offsets(site_seq, plus_strand_pattern):
    """Return the 0-based offsets (left to right on the plus strand) where they differ."""
    return [i for i, (a, b) in enumerate(zip(site_seq, plus_strand_pattern)) if a != b]


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
#   "+" strand site: the primer itself appears on the plus strand.
#                    It extends to the RIGHT.
#   "-" strand site: the primer's reverse complement appears on the plus
#                    strand. It extends to the LEFT.
# =============================================================================

def find_binding_sites(primer_name, primer, target, index, k, max_mismatches):
    """Return every place the primer binds with at most max_mismatches mismatches."""
    primer = primer.upper()
    pieces = max_mismatches + 1
    piece_length = len(primer) // pieces
    if piece_length < k:
        raise ValueError(
            f"{primer_name}: {len(primer)} bases is too short to seed with k={k} "
            f"while allowing {max_mismatches} mismatches; lower k or the tolerance."
        )

    sites = []
    for strand, pattern in (("+", primer), ("-", reverse_complement(primer))):
        # -- 4a: collect candidate start positions from exact seed hits
        candidates = set()
        for piece in range(pieces):
            offset = piece * piece_length
            seed = pattern[offset:offset + k]
            for hit in index.get(seed, []):
                start = hit - offset
                if 0 <= start <= len(target) - len(pattern):
                    candidates.add(start)

        # -- 4b: verify each candidate over the full primer length
        for start in sorted(candidates):
            site_seq = target[start:start + len(pattern)]
            offsets = mismatch_offsets(site_seq, pattern)
            if len(offsets) <= max_mismatches:
                sites.append({
                    "primer": primer_name,
                    "strand": strand,
                    "start": start,                    # 0-based, plus-strand coords
                    "end": start + len(pattern),       # exclusive
                    "mismatches": len(offsets),
                    "mismatch_from_3prime": distances_from_3prime(offsets, len(pattern), strand),
                })
    return sites


# =============================================================================
# SECTION 5: PAIRING SITES INTO PCR PRODUCTS
# -----------------------------------------------------------------------------
# A product forms when a "+" site (extending right) sits upstream of a "-"
# site (extending left) and the distance between their outer ends is within
# a plausible product length. Any primer can pair with any other primer, or
# with itself, because that is what happens in a real tube.
# =============================================================================

def pair_into_products(sites, min_length, max_length):
    """Return every product formed by a '+' site facing a downstream '-' site."""
    plus_sites = [s for s in sites if s["strand"] == "+"]
    minus_sites = [s for s in sites if s["strand"] == "-"]
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
# SECTION 6: THE WHOLE MODEL IN ONE CALL
# -----------------------------------------------------------------------------
# Primer pair + target sequence in, predicted products out.
# =============================================================================

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
    """Print the binding sites and products in a readable form."""
    print(f"\n--- {title} ---")
    print(f"binding sites: {len(sites)}")
    for s in sites:
        print(f"  {s['primer']:<8} strand {s['strand']}  {s['start']:>4}-{s['end']:<4} "
              f"mismatches={s['mismatches']} from 3' end={s['mismatch_from_3prime']}")
    print(f"products: {len(products)}")
    for p in products:
        print(f"  {p['left']['primer']} x {p['right']['primer']}: "
              f"{p['start']}-{p['end']}  length {p['length']} bp")


# =============================================================================
# SECTION 7: FAKE PRIMERS AND FAKE TARGET BUILDING BLOCKS
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
    """Return seq with the base at offset replaced (used to plant mismatches)."""
    return seq[:offset] + new_base + seq[offset + 1:]


# =============================================================================
# SECTION 8: CASES SOLVABLE BY HAND
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


# =============================================================================
# SECTION 9: RUN EVERY CASE
# =============================================================================

if __name__ == "__main__":
    case_exact_match()
    case_reverse_strand()
    case_mismatches()
    case_mismatch_5prime()
    case_multiple_products()
    case_facing_away()
    case_no_site()
    print("\nAll hand-solvable cases passed.")
