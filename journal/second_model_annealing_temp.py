# =============================================================================
# SECOND MODEL: FACTOR: ANNEALING TEMPERATURE   (version 1)
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
#
# What it does NOT do yet:
#   - insertions/deletions: no bulge-loop parameters yet (Tm "not computed")
#   - two mismatches next to each other: no tandem-mismatch parameters
#   - hairpins and primer dimers
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
