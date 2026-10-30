# =============================================================================
# VARIANT CHECK: FIRST MODEL: FACTOR: SPECIFICITY (quick first pass)
# =============================================================================
# The project's actual question: do the CDC N1, N2 and N3 assays still match
# the Alpha, Delta and Omicron genomes? Each assay is run on the reference and
# on every genome in data/variants/, and one summary table is printed.
#
# For each oligo the table shows its fewest mismatches at any site (0 = perfect
# match, "-" = no site within the allowance), plus the assay outcome, product
# size and 3' risk.
#
# This is a premature first pass: one isolate per variant, and no melting
# temperature yet.
#
# Run it with:  python3 journal/variant_check.py
# =============================================================================

import glob
import os

from cdc_reference_check import MAX_MISMATCHES, PRIMER_TABLE, PROJECT, REFERENCE, read_assays
from first_model_factor_specificity import read_fasta, run_model


# =============================================================================
# SECTION 1: GENOMES
# -----------------------------------------------------------------------------
# The reference first, then each variant named by the start of its file name
# (e.g. "Alpha_B.1.1.7_MZ068155.fasta" -> "Alpha").
# =============================================================================

def load_genomes():
    """Return [(label, sequence)], reference first."""
    genomes = [("Reference", read_fasta(REFERENCE)[1])]
    for path in sorted(glob.glob(os.path.join(PROJECT, "data", "variants", "*.fasta"))):
        genomes.append((os.path.basename(path).split("_")[0], read_fasta(path)[1]))
    return genomes


# =============================================================================
# SECTION 2: ONE ROW PER ASSAY PER GENOME
# =============================================================================

def fewest_mismatches(result, oligo):
    """Fewest mismatches at any binding site for this oligo, or '-' if none."""
    counts = [s["mismatches"] for s in result["sites"]
              if s["primer"] == oligo and s["outcome"] == "binds"]
    return str(min(counts)) if counts else "-"


def summarize(assay, oligos, label, genome):
    result = run_model(oligos["forward"]["sequence_5to3"], oligos["reverse"]["sequence_5to3"],
                       genome, probe=oligos["probe"]["sequence_5to3"], max_mismatches=MAX_MISMATCHES)
    products = [p for p in result["products"] if p["status"] == "product"]
    return {
        "assay": assay,
        "genome": label,
        "F": fewest_mismatches(result, "forward"),
        "R": fewest_mismatches(result, "reverse"),
        "P": fewest_mismatches(result, "probe"),
        "outcome": result["assay"],
        "size": ", ".join(str(p["length"]) for p in products) or "-",
        "risk": ", ".join(p["three_prime_risk"] for p in products if p["three_prime_risk"]) or "-",
    }


# =============================================================================
# SECTION 3: RUN
# =============================================================================

if __name__ == "__main__":
    assays = read_assays(PRIMER_TABLE)
    genomes = load_genomes()
    print(f"Mismatch allowance: {MAX_MISMATCHES}. F/R/P = fewest mismatches for "
          f"forward/reverse/probe ('-' = no site).\n")
    columns = ("assay", "genome", "F", "R", "P", "outcome", "size", "risk")
    widths = (6, 10, 3, 3, 3, 15, 8, 10)
    print("".join(c.ljust(w) for c, w in zip(columns, widths)))
    for name in ("N1", "N2", "N3"):
        for label, genome in genomes:
            row = summarize(name, assays[name], label, genome)
            print("".join(row[c].ljust(w) for c, w in zip(columns, widths)))
        print()
