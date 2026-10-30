# =============================================================================
# CDC REFERENCE CHECK: FIRST MODEL: FACTOR: SPECIFICITY
# =============================================================================
# The first check in the proposal ("Reproduction of a published value"):
# the CDC panel's three SARS-CoV-2 assays (N1, N2, N3) must each give exactly
# ONE product on the reference genome NC_045512.2, at the amplicon sizes
# published in Lu et al. (2020), with the probe binding inside it.
#
# The human RNase P assay (RP) is a negative control: it targets a human gene,
# so on the virus genome it must give NO product.
#
# Published positions are offset by 16 bases from NC_045512.2, so the
# positions are checked against the "verified" columns of the primer table,
# which were located on NC_045512.2 itself.
#
# Run it with:  python3 journal/cdc_reference_check.py
# =============================================================================

import csv
import glob
import os

from first_model_factor_specificity import print_report, read_fasta, run_model


# =============================================================================
# SECTION 1: FILE LOCATIONS
# -----------------------------------------------------------------------------
# Paths are built from this file's location, so the script works from any
# folder. The reference genome is the only file in data/reference/.
# =============================================================================

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRIMER_TABLE = os.path.join(PROJECT, "data", "primers", "CDC primer sequences.tsv")
REFERENCE = glob.glob(os.path.join(PROJECT, "data", "reference", "*"))[0]

MAX_MISMATCHES = 2   # the model's default allowance


# =============================================================================
# SECTION 2: READING THE PRIMER TABLE
# -----------------------------------------------------------------------------
# One row per oligo (forward, reverse, probe); rows are grouped by assay.
# =============================================================================

def read_assays(path):
    """Return {assay: {"forward": row, "reverse": row, "probe": row}}."""
    assays = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row = {key.strip(): value.strip() for key, value in row.items()}
            assays.setdefault(row["assay"], {})[row["oligo"]] = row
    return assays


# =============================================================================
# SECTION 3: CHECKS THAT KEEP GOING
# -----------------------------------------------------------------------------
# Unlike the test file, a failed check here does not stop the program: every
# assay is checked and a summary is printed at the end.
# =============================================================================

failures = []


def check(condition, message):
    print(f"  {'PASS' if condition else 'FAIL'}: {message}")
    if not condition:
        failures.append(message)


# =============================================================================
# SECTION 4: THE CHECK FOR EACH ASSAY
# =============================================================================

def check_viral_assay(name, oligos, genome):
    """N1/N2/N3: one product, published size, verified position, probe binds."""
    result = run_model(oligos["forward"]["sequence_5to3"], oligos["reverse"]["sequence_5to3"],
                       genome, probe=oligos["probe"]["sequence_5to3"], max_mismatches=MAX_MISMATCHES)
    print_report(f"{name} on the reference genome", result)

    products = [p for p in result["products"] if p["status"] == "product"]
    published = int(oligos["forward"]["pub_amplicon_bp"])
    first = int(oligos["forward"]["verified_start_NC045512.2"])
    last = int(oligos["reverse"]["verified_start_NC045512.2"])   # reverse primer's 5' end

    check(len(products) == 1, f"{name}: exactly one product")
    if not products:
        return
    product = products[0]
    check((product["start"] + 1, product["end"]) == (first, last),
          f"{name}: product at {product['start'] + 1}-{product['end']}, "
          f"verified primer positions {first}-{last}")
    check(product["length"] == published,
          f"{name}: product {product['length']} bp, published {published} bp")
    check(product["probe"] == "binds", f"{name}: probe binds inside the product")
    check(product["left"]["mismatches"] == 0 and product["right"]["mismatches"] == 0,
          f"{name}: both primers match the reference perfectly")


def check_negative_control(name, oligos, genome):
    """RP: a human gene, so no product on the virus genome."""
    result = run_model(oligos["forward"]["sequence_5to3"], oligos["reverse"]["sequence_5to3"],
                       genome, probe=oligos["probe"]["sequence_5to3"], max_mismatches=MAX_MISMATCHES)
    print_report(f"{name} (negative control) on the reference genome", result)
    check(result["assay"] == "no product", f"{name}: no product on the virus genome")


# =============================================================================
# SECTION 5: RUN
# =============================================================================

if __name__ == "__main__":
    header, genome = read_fasta(REFERENCE)
    print(f"Reference: {header}")
    print(f"Length: {len(genome)} bases, mismatch allowance: {MAX_MISMATCHES}")

    assays = read_assays(PRIMER_TABLE)
    for name in ("N1", "N2", "N3"):
        check_viral_assay(name, assays[name], genome)
    check_negative_control("RP", assays["RP"], genome)

    print()
    if failures:
        print(f"{len(failures)} check(s) failed:")
        for message in failures:
            print(f"  - {message}")
    else:
        print("All CDC reference checks passed.")
