# =============================================================================
# VARIANT CHECK: SECOND MODEL: FACTOR: ANNEALING TEMPERATURE
# =============================================================================
# The same question as variant_check.py, one step further: the CDC N1, N2
# and N3 oligos still BIND the variants (first model), but how much did each
# variant mutation weaken them, and are they still stuck at the CDC
# annealing temperature (55 C)?
#
# For every oligo on every genome the table shows its Tm at the site inside
# the product, the change from a perfect match (dTm), and the percentage of
# sites bound at 55 C. The last column is the product status at 55 C.
#
# FLAG: the "Beta" genome (MZ068155) was first saved as Alpha (see #42).
# One isolate per variant.
#
# Run it with:  python3 journal/variant_check_temperature.py
# =============================================================================

from cdc_reference_check import MAX_MISMATCHES, PRIMER_TABLE, read_assays
from second_model_annealing_temp import ANNEALING_TEMP_C, run_temperature_model
from variant_check import load_genomes


def describe(site):
    """'Tm (dTm) bound%' for one site, or why it was not computed."""
    if site is None:
        return "-"
    temperature = site["temperature"]
    if temperature["tm"] is None:
        return "not computed"
    dtm = "" if temperature["dTm"] is None or abs(temperature["dTm"]) < 0.05 else f" ({temperature['dTm']:+.1f})"
    return f"{temperature['tm']:.1f}{dtm} {100 * temperature['fraction_bound']:.0f}%"


def probe_site(result, product):
    """The best-bound probe site inside a product."""
    inside = [s for s in result["sites"] if s["primer"] == "probe" and s["outcome"] == "binds"
              and s["start"] >= product["start"] and s["end"] <= product["end"]]
    computed = [s for s in inside if s["temperature"]["tm"] is not None]
    if computed:
        return max(computed, key=lambda s: s["temperature"]["fraction_bound"])
    return inside[0] if inside else None


def summarize(assay, oligos, label, genome):
    result = run_temperature_model(oligos["forward"]["sequence_5to3"], oligos["reverse"]["sequence_5to3"],
                                   genome, probe=oligos["probe"]["sequence_5to3"],
                                   max_mismatches=MAX_MISMATCHES)
    products = [p for p in result["products"] if p["status"] in ("product", "at risk")]
    row = {"assay": assay, "genome": label, "F": "-", "R": "-", "P": "-", "at 55 C": result["assay_at_ta"]}
    if products:
        product = products[0]
        row["F"], row["R"] = describe(product["left"]), describe(product["right"])
        row["P"] = describe(probe_site(result, product))
    return row


if __name__ == "__main__":
    assays = read_assays(PRIMER_TABLE)
    genomes = load_genomes()
    print(f"Annealing at {ANNEALING_TEMP_C:.0f} C, mismatch allowance {MAX_MISMATCHES}. "
          f"Each oligo: Tm in C (change from a perfect match) and % of sites bound at "
          f"{ANNEALING_TEMP_C:.0f} C.\n")
    columns = ("assay", "genome", "F", "R", "P", "at 55 C")
    widths = (6, 10, 18, 18, 18, 14)
    print("".join(c.ljust(w) for c, w in zip(columns, widths)))
    for name in ("N1", "N2", "N3"):
        for label, genome in genomes:
            row = summarize(name, assays[name], label, genome)
            print("".join(row[c].ljust(w) for c, w in zip(columns, widths)))
        print()
