# =============================================================================
# SECOND MODEL: FACTOR: ANNEALING TEMPERATURE   (version 3)
# =============================================================================
# 02-120 project: Do early SARS-CoV-2 qPCR primers still match later variants?
#
# The first model (first_model_factor_specificity.py) answers "is there a site
# within the mismatch allowance?" with yes/no. But a primer or probe only works
# if it is actually STUCK to its site at the annealing temperature of the PCR.
# A mismatch weakens the duplex and lowers its melting temperature (Tm); if
# the Tm drops close to (or below) the annealing temperature, most copies of
# the oligo float free and the product (or the signal) is delayed or lost.
#
# What it does:
#   - runs the first model to find every binding site and product
#   - for every site, builds the duplex: the oligo (5'->3') on top, the
#     template base it actually faces below (mismatches included)
#   - adds up nearest-neighbor dH and dS (SantaLucia 1998), using mismatch
#     parameters for internal and terminal mismatches (Allawi & SantaLucia;
#     Peyret; SantaLucia & Peyret)
#   - corrects for salt (Na+/K+, Mg2+, dNTPs) and computes the Tm exactly the
#     way Primer3 does (checked to 0.01 C on all 11 CDC oligos)
#   - at the ANNEALING temperature, computes the fraction of template sites
#     that have an oligo bound
#   - reports each site's Tm, how much the variant lowered it (dTm), its margin
#     above the annealing temperature, and the fraction bound
#   - adds a temperature-aware product status: a product whose primer is
#     mostly unbound at the annealing temperature is "weak binding"; one whose
#     probe is mostly unbound gives "no signal"
#   - (version 2) checks whether a design's temperatures suit qPCR at all
#     (Section 9): primer Tm 6-13 C above annealing, primers within 4 C,
#     probe 2 C above the primers, product melts at the 95 C step; each
#     oligo at the CDC's own concentration (500 nM primers, 125 nM probe)
#   - (version 3) finds primer dimers (Section 10): every pair of oligos
#     lined up at every offset, the paired region's nearest-neighbor energy
#     summed; salt and oligo amounts can be set by the user (QPCR_CONDITIONS)
#
# What it does NOT do yet:
#   - insertions/deletions: no bulge-loop parameters yet (Tm "not computed")
#   - two mismatches next to each other: no tandem-mismatch parameters
#   - hairpins (an oligo folding on itself needs loop parameters)
#   - a full kinetic/efficiency model: the fraction bound is an equilibrium
#     estimate, not a predicted Ct
#
# Files that go with this model (all in journal/):
#   test_second_model_annealing_temp.py   - hand-solvable cases + Primer3 check
#   variant_check_temperature.py          - CDC N1/N2/N3 Tm on every variant
#
# Run the tests with:     python3 journal/test_second_model_annealing_temp.py
# Run the variant table:  python3 journal/variant_check_temperature.py
# =============================================================================

import math
import os

from first_model_factor_specificity import (
    COMPLEMENT, clean_oligo, genome_vs_primer, reverse_complement, run_model,
)


# =============================================================================
# SECTION 1: THE PARAMETER TABLES (data/thermodynamics/)
# -----------------------------------------------------------------------------
# A DNA duplex is held together less by single base pairs than by the
# STACKING of each pair on its neighbor, so the energy is summed over
# neighboring pairs ("nearest neighbors"). Each table row is one stack:
#
#   key "XY/ZW" = 5'-XY-3' on the top strand facing 3'-ZW-5' on the bottom
#   e.g. "CA/GT" = top C then A, bottom G then T (both pairs Watson-Crick)
#
# The same stack read from the other strand is the key reversed ("TG/AC"
# reversed is "CA/GT"), so each stack is stored only once.
#
#   matched stacks     - santalucia_1998_nn.tsv (10 stacks + initiation)
#   internal mismatch  - one of the two pairs is wrong (51 stacks)
#   terminal mismatch  - the wrong pair is the very last pair of the duplex
#
# dH in kcal/mol, dS in cal/(K*mol), all measured in 1 M NaCl.
# =============================================================================

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
THERMO_DIR = os.path.join(PROJECT, "data", "thermodynamics")


def read_nn_table(filename):
    """Return {key: (dH, dS)} from one parameter table."""
    table = {}
    with open(os.path.join(THERMO_DIR, filename)) as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("sequence") or not line.strip():
                continue
            key, dh, ds = line.split("\t")
            table[key] = (float(dh), float(ds))
    return table


MATCHED = read_nn_table("santalucia_1998_nn.tsv")
INTERNAL_MISMATCH = read_nn_table("allawi_santalucia_internal_mismatch_nn.tsv")
TERMINAL_MISMATCH = read_nn_table("santalucia_peyret_terminal_mismatch_nn.tsv")

R = 1.9872   # gas constant, cal/(K*mol)


def stack_value(top, bottom, table):
    """(dH, dS) for one 2-base stack, read from either strand, or None if not in the table."""
    key = f"{top}/{bottom}"
    if key in table:
        return table[key]
    return table.get(key[::-1])


# =============================================================================
# SECTION 2: REACTION CONDITIONS AND THE ANNEALING TEMPERATURE
# -----------------------------------------------------------------------------
# Salt shields the negative phosphate backbones from each other, so the
# duplex is MORE stable in salt than in water. The tables are for 1 M NaCl;
# PCR buffer has ~50 mM K+ and a few mM Mg2+.
#
#   Mg2+ is converted to an equivalent amount of monovalent ions
#   (von Ahsen et al. 2001): Na_eq = mono + 120 * sqrt(Mg - dNTP).
#   dNTPs grab Mg2+ one-to-one, so only the free Mg2+ counts.
#
#   The entropy is then corrected (SantaLucia 1998):
#   dS += 0.368 * (N - 1) * ln(Na_eq in mol/L)     N = oligo length
#
# These are the same rules Primer3 uses by default, and the same conditions
# as data/primers/primer3_reference_tm.tsv.
#
# The CDC protocol (2019-nCoV RT-PCR panel) cycles 95 C for 3 s, then
# 55 C for 30 s (annealing + extension in one step), 45 cycles, so 55 C is
# the temperature at which primers AND the probe must stay bound.
# (CDC uses 500 nM primers and 125 nM probe; the 250 nM below matches the
# Primer3 reference values, which is a difference of about 1 C in Tm.)
# =============================================================================

# *****************************************************************************
# *****************************************************************************
# ***                                                                       ***
# ***   >>> USER SETTINGS: TEMPERATURE <<<                                  ***
# ***                                                                       ***
# ***   ANNEALING_TEMP_C   = the annealing step of the PCR (CDC: 55 C)      ***
# ***   MIN_FRACTION_BOUND = an oligo counts as holding on at the           ***
# ***                        annealing temperature if at least this         ***
# ***                        fraction of its sites is bound (default 0.5)   ***
# ***   CONDITIONS         = salt and oligo concentrations                  ***
# ***                                                                       ***
# *****************************************************************************
# *****************************************************************************
ANNEALING_TEMP_C = 55.0
MIN_FRACTION_BOUND = 0.5
CONDITIONS = {
    "monovalent_mM": 50.0,   # K+ / Na+
    "divalent_mM": 3.0,      # Mg2+
    "dntp_mM": 0.8,          # dNTPs (each binds one Mg2+)
    "oligo_nM": 250.0,       # primer or probe concentration
}


def equivalent_sodium_molar(conditions):
    """Monovalent-equivalent salt in mol/L (von Ahsen et al. 2001)."""
    free_mg = conditions["divalent_mM"] - conditions["dntp_mM"]
    sodium_mM = conditions["monovalent_mM"] + (120 * math.sqrt(free_mg) if free_mg > 0 else 0)
    return sodium_mM / 1000


# =============================================================================
# SECTION 3: ENERGY OF ONE DUPLEX
# -----------------------------------------------------------------------------
# The duplex is written as two aligned strings:
#   top    = the oligo, 5'->3'
#   bottom = the template base facing each oligo base, 3'->5'
# A position is paired when bottom = complement(top).
#
# Steps:
#   1. A mismatch at either END of the duplex is a terminal mismatch: its
#      stack comes from the terminal table and the duplex is then treated as
#      one base shorter on that side.
#   2. Every stack between the two end pairs comes from the matched table
#      (both pairs right) or the internal-mismatch table (one pair wrong).
#      Two wrong pairs side by side have no parameters: not computed.
#   3. Initiation: one term per END of the duplex, depending on whether the
#      end pair is G-C or A-T (SantaLucia 1998).
#   4. A self-complementary oligo (a perfect palindrome) gets the symmetry
#      correction.
# =============================================================================

def duplex_energy(top, bottom):
    """Return (dH, dS, None) at 1 M NaCl, or (None, None, reason) if it cannot be computed."""
    n = len(top)
    paired = [bottom[i] == COMPLEMENT[top[i]] for i in range(n)]
    dh = ds = 0.0
    first, last = 0, n - 1

    # -- step 1: terminal mismatches
    if not paired[0]:
        if n < 3 or not paired[1]:
            return None, None, "two mismatches side by side at the 5' end"
        # read from the other strand so the mismatch is at the key's 3' end
        value = TERMINAL_MISMATCH.get(f"{bottom[1]}{bottom[0]}/{top[1]}{top[0]}")
        dh, ds = dh + value[0], ds + value[1]
        first = 1
    if not paired[n - 1]:
        if n < 3 or not paired[n - 2]:
            return None, None, "two mismatches side by side at the 3' end"
        value = TERMINAL_MISMATCH.get(f"{top[n - 2:]}/{bottom[n - 2:]}")
        dh, ds = dh + value[0], ds + value[1]
        last = n - 2

    # -- step 2: stacks between the end pairs
    for i in range(first, last):
        if not paired[i] and not paired[i + 1]:
            return None, None, "two mismatches side by side (no tandem-mismatch parameters)"
        table = MATCHED if paired[i] and paired[i + 1] else INTERNAL_MISMATCH
        value = stack_value(top[i:i + 2], bottom[i:i + 2], table)
        dh, ds = dh + value[0], ds + value[1]

    # -- step 3: initiation at each end pair
    for end in (top[first], top[last]):
        value = MATCHED["init_w_term_GC"] if end in "GC" else MATCHED["init_w_term_AT"]
        dh, ds = dh + value[0], ds + value[1]

    # -- step 4: symmetry
    if all(paired) and top == reverse_complement(top):
        dh, ds = dh + MATCHED["symmetry_correction"][0], ds + MATCHED["symmetry_correction"][1]
    return dh, ds, None


# =============================================================================
# SECTION 4: MELTING TEMPERATURE AND FRACTION BOUND
# -----------------------------------------------------------------------------
# Tm (Primer3's definition, two strands at equal concentration C):
#
#   Tm = 1000 * dH / (dS + R * ln(C / 4)) - 273.15        (dS salt-corrected)
#
# (C instead of C/4 for a self-complementary oligo.)
#
# In a real PCR the primer is in huge excess over the template, so what
# matters at the annealing temperature T is the fraction of TEMPLATE sites
# that have an oligo on them:
#
#   dG = 1000 * dH - T * dS             (cal/mol, T in kelvin)
#   K  = exp(-dG / (R * T))             (binding constant)
#   fraction bound = K*C / (1 + K*C)    (C = oligo concentration)
#
# Half the sites are bound when K*C = 1, i.e. at
# 1000 * dH / (dS + R * ln C), a little ABOVE the Primer3 Tm because the
# oligo does not run out. Each 1 C drop in Tm costs roughly 1 C of margin
# above the annealing temperature.
# =============================================================================

def salt_corrected_entropy(ds, length, conditions):
    return ds + 0.368 * (length - 1) * math.log(equivalent_sodium_molar(conditions))


def melting_temperature(dh, ds, length, conditions, self_complementary=False):
    """Tm in C from 1 M NaCl dH and dS, as Primer3 computes it."""
    ds = salt_corrected_entropy(ds, length, conditions)
    concentration = conditions["oligo_nM"] * 1e-9
    if not self_complementary:
        concentration /= 4
    return 1000 * dh / (ds + R * math.log(concentration)) - 273.15


def fraction_bound(dh, ds, length, conditions, temp_c):
    """Fraction of template sites with an oligo bound at temp_c (oligo in excess)."""
    ds = salt_corrected_entropy(ds, length, conditions)
    kelvin = temp_c + 273.15
    dg = 1000 * dh - kelvin * ds
    kc = math.exp(-dg / (R * kelvin)) * conditions["oligo_nM"] * 1e-9
    return kc / (1 + kc)


# =============================================================================
# SECTION 5: ONE BINDING SITE FROM THE FIRST MODEL
# -----------------------------------------------------------------------------
# The genome at the site is turned into the oligo's own direction (reverse
# complement for a "-" site); the template base facing each oligo base is its
# complement. Degenerate oligo letters (e.g. Y in the N3 probe) take the
# matching base: a degenerate oligo is a mix, and the version that matches
# is the one that binds best.
#
# Tm is NOT computed (reason given) when the site has an insertion/deletion,
# runs off the sequence, faces an unknown genome base, or has a mismatch at a
# degenerate oligo letter.
# =============================================================================

def site_duplex(oligo, target, site):
    """Return (top, bottom, None) for an oligo on its site, or (None, None, reason)."""
    if site["indel"]:
        return None, None, "insertion/deletion site (no bulge parameters yet)"
    if site["start"] < 0 or site["end"] > len(target):
        return None, None, "site runs off the sequence"
    genome = target[site["start"]:site["end"]]
    if site["strand"] == "-":
        genome = reverse_complement(genome)

    top, bottom = [], []
    for primer_base, genome_base in zip(oligo, genome):
        verdict = genome_vs_primer(genome_base, primer_base)
        if genome_base not in "ACGT" or verdict == "unknown":
            return None, None, "site has unknown or ambiguous genome bases"
        if verdict == "match":
            top.append(genome_base)
        elif primer_base in "ACGT":
            top.append(primer_base)
        else:
            return None, None, "mismatch at a degenerate oligo letter"
        bottom.append(COMPLEMENT[genome_base])
    return "".join(top), "".join(bottom), None


def perfect_tm(oligo, conditions):
    """Tm of the oligo on a perfectly matching site (degenerate letters: not computed)."""
    if any(base not in "ACGT" for base in oligo):
        return None
    bottom = "".join(COMPLEMENT[base] for base in oligo)
    dh, ds, _ = duplex_energy(oligo, bottom)
    return melting_temperature(dh, ds, len(oligo), conditions,
                               self_complementary=oligo == reverse_complement(oligo))


def site_temperature(oligo, target, site, annealing_c, conditions):
    """Tm, dTm, margin and fraction bound for one site (or the reason it is not computed)."""
    top, bottom, reason = site_duplex(oligo, target, site)
    if reason is None:
        dh, ds, reason = duplex_energy(top, bottom)
    if reason is not None:
        return {"tm": None, "reason": reason}

    palindrome = top == reverse_complement(top) and bottom == "".join(COMPLEMENT[b] for b in top)
    tm = melting_temperature(dh, ds, len(top), conditions, self_complementary=palindrome)
    reference = perfect_tm(oligo, conditions)
    bound = fraction_bound(dh, ds, len(top), conditions, annealing_c)
    return {
        "tm": tm,
        "dTm": None if reference is None else tm - reference,   # change caused by mismatches
        "margin": tm - annealing_c,                              # how far above annealing
        "fraction_bound": bound,
        "holds": bound >= MIN_FRACTION_BOUND,
        "dH": dh,
        "dS": ds,
        "reason": "",
    }


# =============================================================================
# SECTION 6: PRODUCT STATUS AT THE ANNEALING TEMPERATURE
# -----------------------------------------------------------------------------
# Starts from the first model's status and only ever makes it worse:
#   product / at risk -> "weak binding" if either primer site is bound less
#                        than MIN_FRACTION_BOUND at the annealing temperature
#                        (the product may be much delayed or not form)
#   product / at risk -> "no signal" if the primers hold but no probe site
#                        inside the product holds (DNA copied, no light)
# A site whose Tm could not be computed is given the benefit of the doubt
# (the first model's status stands) and is listed under "notes".
#
# STATUS_ORDER_AT_TA adds "weak binding" just below "at risk":
#   product > at risk > weak binding > not assessable > no signal > no product
# =============================================================================

STATUS_ORDER_AT_TA = ["product", "at risk", "weak binding", "not assessable",
                      "no signal", "no product"]


def holds(site):
    """False only if the Tm was computed and the site is mostly unbound."""
    temperature = site["temperature"]
    return temperature["tm"] is None or temperature["holds"]


def status_at_annealing(product, sites, probe_given):
    status = product["status"]
    if status not in ("product", "at risk"):
        return status
    if not (holds(product["left"]) and holds(product["right"])):
        return "weak binding"
    if probe_given:
        probes = [s for s in sites if s["primer"] == "probe" and s["outcome"] == "binds"
                  and s["start"] >= product["start"] and s["end"] <= product["end"]]
        if probes and not any(holds(s) for s in probes):
            return "no signal"
    return status


def assay_at_annealing(products):
    statuses = {p["status_at_ta"] for p in products}
    for status in STATUS_ORDER_AT_TA:
        if status in statuses:
            return status
    return "no product"


# =============================================================================
# SECTION 7: THE WHOLE MODEL IN ONE CALL
# =============================================================================

def run_temperature_model(forward, reverse, target, probe=None, annealing_c=ANNEALING_TEMP_C,
                          conditions=CONDITIONS, **first_model_options):
    """First model + Tm and fraction bound at the annealing temperature for every site."""
    result = run_model(forward, reverse, target, probe=probe, **first_model_options)
    target = target.upper()
    oligos = {"forward": clean_oligo("forward", forward), "reverse": clean_oligo("reverse", reverse)}
    if probe:
        oligos["probe"] = clean_oligo("probe", probe)

    for site in result["sites"]:
        site["temperature"] = site_temperature(oligos[site["primer"]], target, site,
                                               annealing_c, conditions)
    for product in result["products"]:
        product["status_at_ta"] = status_at_annealing(product, result["sites"], probe is not None)

    result["annealing_c"] = annealing_c
    result["perfect_tm"] = {name: perfect_tm(seq, conditions) for name, seq in oligos.items()}
    result["assay_at_ta"] = assay_at_annealing(result["products"])
    result["notes"] = [f"{s['primer']} site {s['start'] + 1}-{s['end']}: Tm not computed "
                       f"({s['temperature']['reason']})"
                       for s in result["sites"] if s["temperature"]["tm"] is None]
    return result


# =============================================================================
# SECTION 8: PRINTING A READABLE REPORT
# =============================================================================

def format_temperature(temperature):
    if temperature["tm"] is None:
        return f"Tm not computed ({temperature['reason']})"
    dtm = "" if temperature["dTm"] is None else f" (dTm {temperature['dTm']:+.1f})"
    return (f"Tm {temperature['tm']:.1f} C{dtm}, margin {temperature['margin']:+.1f} C, "
            f"bound {100 * temperature['fraction_bound']:.1f}%")


def print_temperature_report(title, result):
    print(f"\n--- {title} (annealing at {result['annealing_c']:.0f} C) ---")
    perfect = ", ".join(f"{name} {'-' if tm is None else f'{tm:.1f}'}"
                        for name, tm in result["perfect_tm"].items())
    print(f"perfect-match Tm: {perfect}")
    for s in result["sites"]:
        print(f"  {s['primer']:<8} strand {s['strand']}  {s['start'] + 1:>5}-{s['end']:<5} "
              f"{s['outcome']:<15} mismatches={s['mismatches']}  {format_temperature(s['temperature'])}")
    for p in result["products"]:
        print(f"  product {p['start'] + 1}-{p['end']} ({p['length']} bp): "
              f"first model {p['status']}, at annealing {p['status_at_ta']}")
    for note in result["notes"]:
        print(f"NOTE: {note}")
    print(f"assay: first model {result['assay']}, at annealing {result['assay_at_ta']}")


# =============================================================================
# SECTION 9: TEMPERATURE RULES FOR A qPCR ASSAY   (added in version 2)
# -----------------------------------------------------------------------------
# Sections 1-8 check ONE assay against a genome. This section asks the
# design question: are these oligos' temperatures right for qPCR at all?
# Each qPCR cycle has two steps (CDC protocol):
#
#   95 C, 3 s   DENATURATION: the two DNA strands melt apart (unwind)
#   55 C, 30 s  ANNEALING + EXTENSION: primers and probe bind, the polymerase
#               copies, and the probe is cut (that is the signal)
#
# For that to work every cycle:
#   1. PRIMERS MUST HOLD AT THE ANNEALING TEMPERATURE. A primer whose Tm is
#      close to 55 C is only partly bound, so fewer templates are copied
#      each cycle. Every CDC primer has a Tm at least 6.6 C above 55 C (so
#      > 99 % bound). Rule: primer Tm >= annealing + 6 C.
#      But not too far above: a primer that binds very strongly also holds
#      on to sites with a few mismatches, so it is less specific.
#      Rule: primer Tm <= annealing + 13 C (CDC's highest is +12.4).
#   2. BOTH PRIMERS ALIKE. One annealing temperature serves both, so their
#      Tms must be close: Rule: within 4 C (CDC N1: 3.7 C; Primer3's default
#      is 3 C at its own conditions).
#   3. PROBE ABOVE THE PRIMERS. The probe must already be sitting on the
#      template when the primer is extended through it, or nothing is cut
#      and there is no light. Rule: probe Tm >= the higher primer Tm + 2 C.
#      (TaqMan guides say 5-10 C; CDC N1-N3 are +5.0 to +6.3, but CDC's
#      RNase P probe is only +2.6, and it works, so 2 C is the floor.)
#   4. THE PRODUCT MUST UNWIND AT THE DENATURATION STEP. The copied product
#      is double-stranded; if its own Tm is near 95 C it does not fully melt
#      apart, so it cannot be copied again. GC-rich and long products melt
#      highest. Rule: product Tm <= denaturation - 5 C.
#
# Concentrations: CDC uses 500 nM of each primer and 125 nM probe (vs the
# 250 nM of the Primer3 check), so each oligo's Tm uses its own amount.
#
# Product Tm (a long duplex; Primer3's formula for products):
#   Tm = 81.5 + 16.6 * log10(Na_eq) + 0.41 * GC% - 600 / length
# with Na_eq the same Mg2+-corrected salt as in Section 2.
#
# The Tm rules here (1-3) are for the perfect-match site. How a VARIANT's
# mismatches lower those Tms is what Sections 5-6 compute.
# =============================================================================

# *****************************************************************************
# ***   >>> USER SETTINGS: qPCR TEMPERATURE RULES <<<                       ***
# *****************************************************************************
DENATURATION_C = 95.0
PRIMER_NM = 500.0                     # CDC: 500 nM each primer
PROBE_NM = 125.0                      # CDC: 125 nM probe
PRIMER_TM_ABOVE_ANNEALING = (6.0, 13.0)   # primer Tm from annealing + 6 to + 13 C
PRIMER_TM_MAX_DIFFERENCE = 4.0
PROBE_ABOVE_PRIMERS_C = 2.0
PRODUCT_BELOW_DENATURATION_C = 5.0


# The qPCR conditions the designer uses. The user can change any of them;
# primer_nM and probe_nM replace oligo_nM (one amount for every oligo),
# which stays for the Primer3 check.
QPCR_CONDITIONS = dict(CONDITIONS, primer_nM=PRIMER_NM, probe_nM=PROBE_NM)
CONDITION_LIMITS = {                   # what the designer accepts (typical qPCR mixes)
    "monovalent_mM": (10.0, 200.0),
    "divalent_mM": (0.0, 10.0),
    "dntp_mM": (0.0, 2.0),
    "primer_nM": (50.0, 2000.0),
    "probe_nM": (50.0, 1000.0),
}


def conditions_for(role, conditions=CONDITIONS):
    """The reaction conditions with this oligo's own concentration
    (conditions' primer_nM / probe_nM if given, else the CDC amounts)."""
    if role == "probe":
        amount = conditions.get("probe_nM", PROBE_NM)
    else:
        amount = conditions.get("primer_nM", PRIMER_NM)
    return dict(conditions, oligo_nM=amount)


def check_conditions(conditions):
    """Raise ValueError for a missing or out-of-range condition."""
    for name, (low, high) in CONDITION_LIMITS.items():
        value = conditions.get(name)
        if not isinstance(value, (int, float)) or value != value or not low <= value <= high:
            raise ValueError(f"{name.replace('_', ' ')} must be {low:g}-{high:g}")


def product_tm(product, conditions=CONDITIONS):
    """Tm of the double-stranded product (Primer3's long-sequence formula)."""
    product = product.upper()
    gc_percent = 100 * sum(base in "GC" for base in product) / len(product)
    return (81.5 + 16.6 * math.log10(equivalent_sodium_molar(conditions))
            + 0.41 * gc_percent - 600 / len(product))


def oligo_at_annealing(oligo, role, annealing_c=ANNEALING_TEMP_C, conditions=CONDITIONS):
    """Perfect-match Tm and fraction bound at the annealing temperature, at the role's amount."""
    own = conditions_for(role, conditions)
    if any(base not in "ACGT" for base in oligo):
        return {"tm": None, "fraction_bound": None}
    dh, ds, _ = duplex_energy(oligo, "".join(COMPLEMENT[b] for b in oligo))
    return {"tm": perfect_tm(oligo, own),
            "fraction_bound": fraction_bound(dh, ds, len(oligo), own, annealing_c)}


def assay_temperature_check(forward, reverse, probe=None, product=None,
                            annealing_c=ANNEALING_TEMP_C, conditions=CONDITIONS, dimers=True):
    """Check rules 1-4 (and primer dimers, Section 10); returns each oligo's Tm
    and fraction bound, the product Tm, the dimers, and a list of problems
    (empty = the assay suits qPCR)."""
    oligos = {"forward": forward.upper(), "reverse": reverse.upper()}
    if probe:
        oligos["probe"] = probe.upper()
    result = {name: oligo_at_annealing(seq, name, annealing_c, conditions)
              for name, seq in oligos.items()}
    problems = []
    low, high = (annealing_c + PRIMER_TM_ABOVE_ANNEALING[0], annealing_c + PRIMER_TM_ABOVE_ANNEALING[1])

    primer_tms = []
    for name in ("forward", "reverse"):
        tm = result[name]["tm"]
        if tm is None:
            continue
        primer_tms.append(tm)
        if tm < low:
            problems.append(f"{name} primer Tm {tm:.2f} C is too low: it should be at least "
                            f"{low:.1f} C ({PRIMER_TM_ABOVE_ANNEALING[0]:.0f} C above the {annealing_c:.0f} C "
                            f"annealing step)")
        elif tm > high:
            problems.append(f"{name} primer Tm {tm:.2f} C is too high (above {high:.1f} C): "
                            f"it would also hold on to sites with mismatches")
    if len(primer_tms) == 2 and abs(primer_tms[0] - primer_tms[1]) > PRIMER_TM_MAX_DIFFERENCE:
        problems.append(f"the primers' Tms differ by {abs(primer_tms[0] - primer_tms[1]):.1f} C "
                        f"(at most {PRIMER_TM_MAX_DIFFERENCE:.0f} C): one annealing temperature "
                        f"cannot suit both")
    if probe and result["probe"]["tm"] is not None and primer_tms:
        needed = max(primer_tms) + PROBE_ABOVE_PRIMERS_C
        if result["probe"]["tm"] < needed:
            problems.append(f"probe Tm {result['probe']['tm']:.2f} C is too low: it should be at least "
                            f"{needed:.2f} C ({PROBE_ABOVE_PRIMERS_C:.0f} C above the higher primer), "
                            f"so it is bound before the primers are extended")
    if product:
        tm = product_tm(product, conditions)
        result["product_tm"] = tm
        if tm > DENATURATION_C - PRODUCT_BELOW_DENATURATION_C:
            problems.append(f"product Tm {tm:.1f} C is too close to the {DENATURATION_C:.0f} C "
                            f"denaturation step: the product may not fully unwind")
    if dimers and all(set(seq) <= set("ACGT") for seq in oligos.values()):
        result["dimers"], dimer_problems = dimer_check(oligos, conditions_for("forward", conditions))
        problems.extend(dimer_problems)
    result["problems"] = problems
    result["ok"] = not problems
    return result



# =============================================================================
# SECTION 10: PRIMER DIMERS   (added in version 3)
# -----------------------------------------------------------------------------
# Two oligos in the same tube can stick to EACH OTHER instead of the
# template. If a primer's 3' end is in such a duplex, the polymerase extends
# it and copies the other oligo: a short "primer-dimer" product that uses up
# primers and, with SYBR-type dyes, gives false signal.
#
# Method (as in the proposal): line the two oligos up antiparallel at EVERY
# relative offset (oligo A 5'->3' on top, oligo B 3'->5' below), and at each
# offset take every paired region: a stretch that starts and ends with a
# Watson-Crick pair. Its energy is summed with the same nearest-neighbor
# tables as everything else (duplex_energy): matched stacks, plus single
# internal mismatches inside the region (two side by side end it). The most
# stable region over all offsets is the dimer.
#
#   dG(T) = 1000 * dH - (T + 273.15) * dS     (cal/mol, dS salt-corrected)
#   reported at 37 C, the usual reference temperature for dimers.
#
# Compared with Primer3's homodimer dG for the CDC oligos
# (data/primers/primer3_reference_tm.tsv), these values are 1-4 kcal/mol
# LESS negative: Primer3 also counts bases dangling past the paired region
# and loops. So the limits below were set on THIS method, using the CDC
# assays (test case T14):
#   most stable dimer of any kind: N2 probe with itself, -5.9 kcal/mol
#   most stable dimer with a primer's 3' end paired: N2 forward on the
#   probe, -3.1 kcal/mol
# Rules (a design breaking either is flagged):
#   any dimer                  dG(37 C) >= -9 kcal/mol (the common guideline)
#   a primer's 3' end paired   dG(37 C) >= -5 kcal/mol (it can be extended)
# A probe's 3' end does not count: TaqMan probes are blocked there.
# =============================================================================

DIMER_TEMP_C = 37.0
DIMER_MAX_DG = -9000.0          # cal/mol
DIMER_3END_MAX_DG = -5000.0     # cal/mol


def dimer_scan(a, b, conditions=CONDITIONS, temp_c=DIMER_TEMP_C, extendable=(True, True)):
    """The most stable paired region between oligos a and b (both 5'->3').

    extendable = (a's 3' end counts, b's 3' end counts) for the 3'-end search.
    Returns {"any": region, "three_prime": region or None}, each region
    {"dG", "dH", "dS", "top", "bottom", "a_start", "b_start", "length"}.
    """
    below = b[::-1]                          # b written 3'->5' under a
    best_any = best_end = None
    log_salt = math.log(equivalent_sodium_molar(conditions))
    kelvin = temp_c + 273.15
    for offset in range(-(len(below) - 1), len(a)):
        # a[i] faces below[i - offset]
        first, last = max(0, offset), min(len(a), len(below) + offset)
        paired = [below[i - offset] == COMPLEMENT.get(a[i]) for i in range(first, last)]
        for start in range(first, last):
            if not paired[start - first]:
                continue
            # grow the region one base at a time; the stacks are summed left
            # to right as in duplex_energy, so the numbers are the same
            dh = ds = 0.0
            for k in range(start, last - 1):
                this, nxt = paired[k - first], paired[k + 1 - first]
                if not this and not nxt:
                    break                    # two mismatches side by side end the region
                table = MATCHED if this and nxt else INTERNAL_MISMATCH
                value = stack_value(a[k:k + 2], below[k - offset:k - offset + 2], table)
                dh, ds = dh + value[0], ds + value[1]
                if not nxt:
                    continue                 # a region must end with a pair
                end = k + 1
                top, bottom = a[start:end + 1], below[start - offset:end - offset + 1]
                region_dh, region_ds = dh, ds
                for base in (top[0], top[-1]):
                    init = MATCHED["init_w_term_GC"] if base in "GC" else MATCHED["init_w_term_AT"]
                    region_dh, region_ds = region_dh + init[0], region_ds + init[1]
                if all(paired[start - first:end - first + 1]) and top == reverse_complement(top):
                    region_dh += MATCHED["symmetry_correction"][0]
                    region_ds += MATCHED["symmetry_correction"][1]
                length = end - start + 1
                dg = 1000 * region_dh - kelvin * (region_ds + 0.368 * (length - 1) * log_salt)
                if best_any is not None and dg >= best_any["dG"] and not (
                        (extendable[0] and end == len(a) - 1) or (extendable[1] and start - offset == 0)):
                    continue
                region = {"dG": dg, "dH": region_dh, "dS": region_ds, "top": top, "bottom": bottom,
                          "a_start": start, "b_start": len(b) - 1 - (end - offset), "length": length}
                if best_any is None or dg < best_any["dG"]:
                    best_any = region
                a_end = extendable[0] and end == len(a) - 1
                b_end = extendable[1] and start - offset == 0
                if (a_end or b_end) and (best_end is None or dg < best_end["dG"]):
                    best_end = region
    return {"any": best_any, "three_prime": best_end}


def dimer_check(oligos, conditions=CONDITIONS, cache=None):
    """Every self- and cross-dimer among {name: seq}; returns (dimers, problems).

    dimers: [{"pair": "forward/reverse", "any": ..., "three_prime": ...}, ...]
    cache: an optional dict, reused across calls with the SAME conditions, so a
    primer that appears in many designs is scanned once per partner.
    """
    names = list(oligos)
    dimers, problems = [], []
    for x, first in enumerate(names):
        for second in names[x:]:
            extendable = (first != "probe", second != "probe")
            key = (oligos[first], oligos[second], extendable)
            if cache is not None and key in cache:
                scan = cache[key]
            else:
                scan = dimer_scan(oligos[first], oligos[second], conditions, extendable=extendable)
                if cache is not None:
                    cache[key] = scan
            pair = first if first == second else f"{first}/{second}"
            dimers.append({"pair": pair, **scan})
            worst = scan["any"]
            if worst and worst["dG"] < DIMER_MAX_DG:
                problems.append(f"{pair} dimer: {worst['length']} bp, dG {worst['dG'] / 1000:.1f} kcal/mol "
                                f"(more stable than {DIMER_MAX_DG / 1000:.0f}): the oligos stick to each other")
            end = scan["three_prime"]
            if end and end["dG"] < DIMER_3END_MAX_DG:
                problems.append(f"{pair} dimer at a primer's 3' end: {end['length']} bp, dG "
                                f"{end['dG'] / 1000:.1f} kcal/mol (more stable than "
                                f"{DIMER_3END_MAX_DG / 1000:.0f}): it can be extended into a primer-dimer")
    return dimers, problems
