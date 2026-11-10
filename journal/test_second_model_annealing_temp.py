# =============================================================================
# TESTS FOR: SECOND MODEL: FACTOR: ANNEALING TEMPERATURE   (version 1)
# =============================================================================
# The model lives in second_model_annealing_temp.py; this file only tests it.
# The numbers in Cases T1-T4 can be added up by hand from the three tables in
# data/thermodynamics/ (every stack used is listed in the comments). Case T5
# reproduces published values: Primer3's Tm for the 11 CDC oligos.
#
# Fake primers and targets are the 18-30 base ones from the first model's
# tests (test_first_model_factor_specificity.py).
#
# Run it with:  python3 journal/test_second_model_annealing_temp.py
# =============================================================================

import csv
import math
import os

from first_model_factor_specificity import COMPLEMENT
from second_model_annealing_temp import (
    CONDITIONS, PROJECT, duplex_energy, equivalent_sodium_molar, fraction_bound,
    melting_temperature, perfect_tm, print_temperature_report, run_temperature_model,
)
from test_first_model_factor_specificity import (
    FORWARD_A, FORWARD_B, PROBE_A, REVERSE_A, amplicon, check, filler, mutate,
)

ONE_MOLAR = {"monovalent_mM": 1000.0, "divalent_mM": 0.0, "dntp_mM": 0.0, "oligo_nM": 250.0}


def template_for(oligo):
    """The perfectly matching template, written 3'->5' under the oligo."""
    return "".join(COMPLEMENT[base] for base in oligo)


def close(a, b, tolerance=0.01):
    return abs(a - b) <= tolerance


# -- Case T1: adding up the stacks by hand -------------------------------------
# FORWARD_B = 5'-CAGTGCTTACGGATCGAG-3' (18 bases, 17 stacks). Each stack is
# looked up in santalucia_1998_nn.tsv, read from the other strand when needed
# (e.g. AG/TC is stored as CT/GA):
#   CA -8.5/-22.7   AG -7.8/-21.0   GT -8.4/-22.4   TG -8.5/-22.7
#   GC -9.8/-24.4   CT -7.8/-21.0   TT -7.9/-22.2   TA -7.2/-21.3
#   AC -8.4/-22.4   CG -10.6/-27.2  GG -8.0/-19.9   GA -8.2/-22.2
#   AT -7.2/-20.4   TC -8.2/-22.2   CG -10.6/-27.2  GA -8.2/-22.2
#   AG -7.8/-21.0
#   stacks:     dH = -143.1   dS = -382.4
#   initiation: 5' end C (G-C pair) +0.1/-2.8, 3' end G (G-C pair) +0.1/-2.8
#   total:      dH = -142.9 kcal/mol, dS = -388.0 cal/(K*mol)
def case_stack_sum():
    print("\n--- Case T1: stack sum by hand, 18-base primer ---")
    dh, ds, reason = duplex_energy(FORWARD_B, template_for(FORWARD_B))
    check(reason is None, "a perfect duplex can be computed")
    check(close(dh, -142.9), f"dH = -142.9 kcal/mol (got {dh:.2f})")
    check(close(ds, -388.0), f"dS = -388.0 cal/(K*mol) (got {ds:.2f})")


# -- Case T2: from dH and dS to Tm by hand ------------------------------------
# a) 1 M NaCl, 250 nM oligo: C/4 = 62.5e-9 M, ln(62.5e-9) = -16.588
#    Tm = -142900 / (-388.0 + 1.9872 * -16.588) - 273.15
#       = -142900 / -420.96 - 273.15 = 339.46 - 273.15 = 66.31 C
# b) PCR buffer: 50 mM K+, 3 mM Mg2+, 0.8 mM dNTP
#    Na_eq = 50 + 120 * sqrt(3 - 0.8) = 50 + 177.99 = 227.99 mM
#    dS = -388.0 + 0.368 * 17 * ln(0.22799) = -388.0 - 9.25 = -397.25
#    Tm = -142900 / (-397.25 - 32.96) - 273.15 = 332.16 - 273.15 = 59.01 C
#    Less salt than 1 M, so the duplex is less stable: Tm drops 7.3 C.
def case_tm_by_hand():
    print("\n--- Case T2: Tm by hand, 1 M NaCl and PCR buffer ---")
    check(close(equivalent_sodium_molar(CONDITIONS) * 1000, 227.99), "Na_eq = 227.99 mM")
    check(close(melting_temperature(-142.9, -388.0, 18, ONE_MOLAR), 66.31), "Tm in 1 M NaCl = 66.31 C")
    check(close(melting_temperature(-142.9, -388.0, 18, CONDITIONS), 59.01), "Tm in PCR buffer = 59.01 C")


# -- Case T3: one internal mismatch -------------------------------------------
# FORWARD_B offset 8 is A. The genome now has C there, so the template base
# facing the primer's A is G (an A-G mismatch, "severe" in the first model).
# Two stacks touch the mismatch; both now come from the internal-mismatch table:
#   removed  TA/AT -7.2/-21.3 and AC/TG (= GT/CA) -8.4/-22.4  -> -15.6/-43.7
#   added    TA/AG +0.7/+0.7  and AC/GG (= GG/CA) +0.5/+3.2   ->  +1.2/ +3.9
#   change   dH +16.8, dS +47.6
#   new      dH = -126.1, dS = -340.4 -> Tm 56.43 C, dTm = -2.59 C
def case_internal_mismatch():
    print("\n--- Case T3: one internal mismatch (A-G at offset 8) ---")
    bottom = mutate(template_for(FORWARD_B), 8, "G")
    dh, ds, _ = duplex_energy(FORWARD_B, bottom)
    check(close(dh, -126.1) and close(ds, -340.4), f"dH = -126.1, dS = -340.4 (got {dh:.2f}, {ds:.2f})")
    tm = melting_temperature(dh, ds, 18, CONDITIONS)
    check(close(tm, 56.43), f"Tm = 56.43 C (got {tm:.2f})")
    check(close(tm - 59.01, -2.59, 0.02), "the mismatch lowers Tm by 2.59 C")


# -- Case T4: a severe mismatch at the 3' end barely changes Tm ---------------
# Same primer, last base G now faces template A (G-A: severe; the first model
# calls this 3' risk HIGH -> "at risk"). The last stack becomes a TERMINAL
# mismatch and the duplex's 3' end pair is now the A one base in:
#   removed  AG/TC (= CT/GA) -7.8/-21.0, initiation at G +0.1/-2.8
#   added    AG/TA (terminal table) -2.5/-5.9, initiation at A +2.3/+4.1
#   new      dH = -135.4, dS = -366.0 -> Tm 58.54 C, dTm = only -0.47 C
# The primer still sticks almost as well, but its last base is not paired, so
# the polymerase may not extend it. That is why the first model's 3' risk and
# this model's Tm answer different questions and both are kept.
def case_terminal_mismatch():
    print("\n--- Case T4: severe mismatch at the very 3' end ---")
    bottom = mutate(template_for(FORWARD_B), 17, "A")
    dh, ds, _ = duplex_energy(FORWARD_B, bottom)
    check(close(dh, -135.4) and close(ds, -366.0), f"dH = -135.4, dS = -366.0 (got {dh:.2f}, {ds:.2f})")
    tm = melting_temperature(dh, ds, 18, CONDITIONS)
    check(close(tm, 58.54), f"Tm = 58.54 C, only 0.47 C lower (got {tm:.2f})")


# -- Case T4b: two mismatches side by side -------------------------------------
# Offsets 8 and 9 both wrong: no table has a stack with two wrong pairs, so
# the Tm is not computed (and the reason says why).
def case_tandem_mismatch():
    print("\n--- Case T4b: two mismatches side by side ---")
    bottom = mutate(mutate(template_for(FORWARD_B), 8, "G"), 9, "A")
    dh, _, reason = duplex_energy(FORWARD_B, bottom)
    check(dh is None and "side by side" in reason, f"not computed: {reason}")


# -- Case T4c: every internal single mismatch lowers Tm ------------------------
# For every fake primer/probe, every position except the two end bases, and
# every wrong base: 3 x (n - 2) duplexes per oligo, Tm always lower.
def case_mismatch_always_lowers():
    print("\n--- Case T4c: every internal single mismatch lowers Tm ---")
    count = 0
    for oligo in (FORWARD_A, FORWARD_B, REVERSE_A, PROBE_A):
        perfect = template_for(oligo)
        reference = perfect_tm(oligo, CONDITIONS)
        for offset in range(1, len(oligo) - 1):
            for wrong in "ACGT":
                if wrong == perfect[offset]:
                    continue
                dh, ds, _ = duplex_energy(oligo, mutate(perfect, offset, wrong))
                if melting_temperature(dh, ds, len(oligo), CONDITIONS) >= reference:
                    raise AssertionError(f"{oligo} offset {offset} {wrong} did not lower Tm")
                count += 1
    # internal positions: 20-2 + 18-2 + 24-2 + 26-2 = 80, times 3 wrong bases
    check(count == 3 * 80, f"all {count} internal mismatches lower Tm")


# -- Case T4d: fraction bound around the annealing temperature -----------------
# Same 18-base primer (dH -142.9, dS -397.25 after salt). Half the template
# sites are bound when K*C = 1, i.e. at
#   T_half = -142900 / (-397.25 + 1.9872 * ln(250e-9)) = 334.30 K = 61.15 C
# a little above its Tm (59.01 C) because the primer is in excess.
#   at 61.15 C -> 50%;  at 55 C (CDC annealing) -> 98.3%;  at 95 C -> ~0%
def case_fraction_bound():
    print("\n--- Case T4d: fraction bound vs temperature ---")
    check(close(fraction_bound(-142.9, -388.0, 18, CONDITIONS, 61.15), 0.5, 0.005), "50% bound at 61.15 C")
    check(close(fraction_bound(-142.9, -388.0, 18, CONDITIONS, 55), 0.983, 0.001), "98.3% bound at 55 C")
    check(fraction_bound(-142.9, -388.0, 18, CONDITIONS, 95) < 1e-6, "essentially nothing bound at 95 C")


# -- Case T5: reproduce Primer3 for the 11 CDC oligos -------------------------
# primer3_reference_tm.tsv was made with primer3-py 2.3.1 at the same
# conditions. Every Tm must match to 0.01 C (the N3 probe has a degenerate Y,
# so Primer3 gave no value for it).
def case_primer3_reference():
    print("\n--- Case T5: Primer3 Tm for the CDC oligos ---")
    path = os.path.join(PROJECT, "data", "primers", "primer3_reference_tm.tsv")
    with open(path) as handle:
        rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    compared = 0
    for row in rows:
        if row["primer3_tm_C"].startswith("NA"):
            continue
        ours = perfect_tm(row["sequence_5to3"], CONDITIONS)
        check(close(ours, float(row["primer3_tm_C"])),
              f"{row['assay']} {row['oligo']}: {ours:.2f} C, Primer3 {row['primer3_tm_C']} C")
        compared += 1
    check(compared == 11, "all 11 oligos with a Primer3 value were compared")


# -- Case T6: product status at the annealing temperature ---------------------
# Set A on target = 10 filler + 70 bp amplicon + 10 filler (product 11-80).
# Perfect-match Tm: forward 64.1, reverse 66.9, probe 69.1 C.
# a) 55 C (CDC): all far above, everything bound -> "product"
# b) 68 C: above the forward primer's Tm; forward only 21% bound
#    -> "weak binding" (the first model still says "product")
def case_status_at_annealing():
    target = filler(10) + amplicon() + filler(10)
    result = run_temperature_model(FORWARD_A, REVERSE_A, target, probe=PROBE_A)
    print_temperature_report("Case T6a: set A at the CDC annealing temperature", result)
    check(result["assay_at_ta"] == "product", "at 55 C: product")

    result = run_temperature_model(FORWARD_A, REVERSE_A, target, probe=PROBE_A, annealing_c=68)
    print_temperature_report("Case T6b: set A annealed too hot", result)
    check(result["assay"] == "product", "first model alone: product")
    check(result["assay_at_ta"] == "weak binding", "at 68 C: weak binding (forward primer melts)")


# -- Case T7: the probe melts, the primers hold -> no signal -------------------
# Two mild G-T mismatches in the probe site (target 38 and 47 = probe offsets
# 8 and 17, genome G -> A). The probe still BINDS for the first model (2
# mismatches allowed), but its Tm falls from 69.1 to 60.1 C. At 64 C:
#   forward 81% bound, reverse 98% bound, probe 14% bound -> "no signal"
def case_probe_melts():
    target = filler(10) + amplicon() + filler(10)
    target = mutate(mutate(target, 38, "A"), 47, "A")
    result = run_temperature_model(FORWARD_A, REVERSE_A, target, probe=PROBE_A, annealing_c=64)
    print_temperature_report("Case T7: probe with 2 mismatches at 64 C", result)
    probe = [s for s in result["sites"] if s["primer"] == "probe"][0]
    check(probe["mismatches"] == 2 and probe["outcome"] == "binds", "first model: probe binds (2 mismatches)")
    check(close(probe["temperature"]["tm"], 60.08), f"probe Tm 60.08 C (got {probe['temperature']['tm']:.2f})")
    check(result["assay"] == "product" and result["assay_at_ta"] == "no signal",
          "first model product, at 64 C no signal")


# -- Case T8: an insertion/deletion site is not computed ----------------------
# The first model's Case 22 (allowance 1): forward site loses amplicon offset 2
# (5' half, so it still binds with a 1-base deletion; product 69 bp). No bulge
# parameters yet, so Tm is not computed, a note says so, and the status stands.
# (With allowance 2 the first model reads the same site without a gap, as 2
# mismatches side by side, which is also not computed.)
def case_indel_not_computed():
    target = filler(10) + amplicon()[:2] + amplicon()[3:] + filler(10)
    result = run_temperature_model(FORWARD_A, REVERSE_A, target, probe=PROBE_A, max_mismatches=1)
    print_temperature_report("Case T8: deletion in the forward site", result)
    forward = [s for s in result["sites"] if s["primer"] == "forward"][0]
    check(forward["indel"] is not None and forward["temperature"]["tm"] is None,
          "forward site has a deletion: Tm not computed")
    check(any("insertion/deletion" in note for note in result["notes"]), "a note explains why")
    check(result["assay_at_ta"] == result["assay"] == "product", "status unchanged: product")


if __name__ == "__main__":
    case_stack_sum()
    case_tm_by_hand()
    case_internal_mismatch()
    case_terminal_mismatch()
    case_tandem_mismatch()
    case_mismatch_always_lowers()
    case_fraction_bound()
    case_primer3_reference()
    case_status_at_annealing()
    case_probe_melts()
    case_indel_not_computed()
    print("\nAll cases passed.")
