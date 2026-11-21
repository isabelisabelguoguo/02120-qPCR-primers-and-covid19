# =============================================================================
# COMPARE: OMICRON SPIKE GENE: website Version 1 vs Version 2 vs Primer-BLAST
# =============================================================================
# Target: the spike (S) gene of Omicron BA.1 (data/variants/OL672836.1),
# genome positions 21497-25309 (3,813 bases; found by matching the first and
# last 30 bases of the reference spike, NC_045512.2 21563-25384; Omicron's
# spike deletions and its 9-base insertion make it 9 bases shorter).
#
# Designers, all with the websites' default settings (mismatch allowance 2,
# products 60-150 bp):
#   Version 1   design_primers(annealing_c=None): GC content only
#   Version 2   design_primers(annealing_c=55): temperature + primer dimers
#   NCBI        Primer-BLAST (Primer3 + BLAST), our length/GC/product
#               settings, its own Tm defaults, specificity vs the human
#               genome. Saved in data/ncbi_blast/primer_blast_omicron_spike.txt
#
# For every design:
#   - Tms at qPCR conditions and Version 2's rules (assay_temperature_check)
#   - the assay status at 55 C on the reference, Beta, Delta and Omicron
#     genomes (run_temperature_model: first model + Tm at each site)
#   - for NCBI's designs, where our two versions rank the same product
#
# Run it with:  python3 journal/compare_omicron_spike.py
# =============================================================================

import glob
import os

from compare_ncbi_blast import primer_blast_designs
from first_model_factor_specificity import read_fasta
from primer_design import design_primers, find_target, possible_designs
from second_model_annealing_temp import QPCR_CONDITIONS, assay_temperature_check, run_temperature_model

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPIKE = (21497, 25309)   # 1-based, inclusive, in OL672836.1


def genomes():
    found = {"Wuhan ref": read_fasta(glob.glob(os.path.join(PROJECT, "data", "reference", "*"))[0])[1]}
    for path in sorted(glob.glob(os.path.join(PROJECT, "data", "variants", "*"))):
        found[os.path.basename(path).split("_")[0]] = read_fasta(path)[1]
    return found


def report(label, forward, reverse, probe, start, end, omicron, all_genomes, extra=""):
    check = assay_temperature_check(forward, reverse, probe, omicron[start - 1:end],
                                    conditions=QPCR_CONDITIONS)
    tms = "/".join(f"{check[k]['tm']:.1f}" for k in ("forward", "reverse", "probe"))
    rules = "pass" if check["ok"] else "; ".join(p.split(":")[0].split(" (")[0] for p in check["problems"])
    status = ", ".join(f"{name} {run_temperature_model(forward, reverse, g, probe=probe, max_mismatches=2)['assay_at_ta']}"
                       for name, g in all_genomes.items())
    print(f"  {label}  {start}-{end} ({end - start + 1} bp){extra}")
    print(f"      F {forward}  R {reverse}  P {probe}")
    print(f"      qPCR Tm F/R/P {tms} C; Version 2 rules: {rules}")
    print(f"      at 55 C: {status}")


if __name__ == "__main__":
    all_genomes = genomes()
    omicron = all_genomes["Omicron"]
    spike = omicron[SPIKE[0] - 1:SPIKE[1]]
    assert spike.startswith("ATG") and spike.endswith("TAA") and len(spike) % 3 == 0
    start, end, _ = find_target(omicron, spike)

    for title, annealing in (("VERSION 1 (GC only)", None), ("VERSION 2 (temperature)", 55.0)):
        result = design_primers(omicron, spike, 2, 60, 150, annealing_c=annealing)
        print(f"\n{title}: {result['possible']} possible designs")
        for n, d in enumerate(result["designs"], 1):
            report(f"design {n}", d["forward"]["seq"], d["reverse"]["seq"], d["probe"]["seq"],
                   d["product_start"] + 1, d["product_end"], omicron, all_genomes)

    ranked = {label: {(d["product_start"], d["product_end"]): i + 1
                      for i, d in reversed(list(enumerate(possible_designs(omicron, start, end,
                                                                           annealing_c=ta))))}
              for label, ta in (("v1", None), ("v2", 55.0))}
    print("\nNCBI PRIMER-BLAST")
    path = os.path.join(PROJECT, "data", "ncbi_blast", "primer_blast_omicron_spike.txt")
    for n, d in enumerate(primer_blast_designs(path), 1):
        f, r, p = d["forward"], d["reverse"], d["probe"]
        product = (start + f[1] - 1, start + r[1])
        ranks = ", ".join(f"{label} rank {ranked[label].get(product, '-')}" for label in ("v1", "v2"))
        report(f"design {n}", f[0], r[0], p[0], product[0] + 1, product[1], omicron, all_genomes,
               f"; in our lists: {ranks}")
