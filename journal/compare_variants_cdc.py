# =============================================================================
# COMPARE: CDC N1/N2/N3 ON ALPHA, BETA AND OMICRON
#          our models vs NCBI BLAST, and the three website versions
# =============================================================================
# The published primers are the CDC's N1, N2 and N3 assays (data/primers/).
# Genomes: the reference NC_045512.2 and one genome each of Alpha
# (OL461702.1, checked for all nine Alpha signature changes and both spike
# deletions), Beta (MZ068155) and Omicron BA.1 (OL672836.1).
#
# PART 1: where each CDC oligo binds on each genome
#   ours  run_temperature_model (first model + Tm): mismatches, where they
#         are counted from the 3' end, the Tm change (dTm), and the assay's
#         status at 55 C
#   NCBI  two-sequence BLAST of the nine oligos against each genome
#         (data/ncbi_blast/cdc_vs_<genome>.tsv), each hit extended to the
#         whole oligo and its mismatches counted (IUPAC letters such as the
#         N3 probe's Y count as matches when the genome base fits)
#
# PART 2: the three website versions on each genome's N gene (found by
#   the closest match, <= 3 mismatches, to the reference N gene's first and
#   last 30 bases)
#   Version 1  design_primers(annealing_c=None)            GC only
#   Version 2  design_primers(annealing_c=55, structure_checks=False)
#   Version 3  design_primers(annealing_c=55)              + dimers, hairpins
#   For each design: whether it overlaps a CDC product, and on which of the
#   four genomes it gives the product at 55 C. For each CDC assay:
#   whether it passes Version 2's and Version 3's rules on that genome.
#
# Run it with:  python3 journal/compare_variants_cdc.py
# =============================================================================

import csv
import glob
import os

from compare_ncbi_blast import read_hits
from first_model_factor_specificity import base_matches, read_fasta, reverse_complement
from primer_design import design_primers
from second_model_annealing_temp import (
    CONDITIONS, QPCR_CONDITIONS, assay_temperature_check, run_temperature_model,
)

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAVED = os.path.join(PROJECT, "data", "ncbi_blast")
N_GENE = (28274, 29533)          # reference, 1-based inclusive
GENOMES = ("reference", "Alpha", "Beta", "Omicron")


def load_genomes():
    found = {"reference": read_fasta(glob.glob(os.path.join(PROJECT, "data", "reference", "*"))[0])[1]}
    for name in GENOMES[1:]:
        found[name] = read_fasta(glob.glob(os.path.join(PROJECT, "data", "variants", f"{name}_*"))[0])[1]
    return found


def cdc_assays():
    rows = csv.DictReader(open(os.path.join(PROJECT, "data", "primers", "CDC primer sequences.tsv")),
                          delimiter="\t")
    assays = {}
    for r in rows:
        if r["assay"] != "RP":
            assays.setdefault(r["assay"], {})[r["oligo"]] = r["sequence_5to3"]
    return assays


def closest(genome, anchor, near, window=400, most=3):
    """Start of the stretch near `near` that differs least from `anchor` (<= most mismatches)."""
    best = None
    for i in range(max(0, near - window), min(len(genome) - len(anchor), near + window) + 1):
        mismatches = sum(a != b for a, b in zip(genome[i:i + len(anchor)], anchor))
        if best is None or mismatches < best[0]:
            best = (mismatches, i)
    if best[0] > most:
        raise ValueError("N gene end not found")
    return best[1]


def n_gene(genome, reference):
    """0-based start and end of the N gene, found from its first and last 30 bases
    (variants such as Alpha D3L change a few of them, so up to 3 may differ)."""
    gene = reference[N_GENE[0] - 1:N_GENE[1]]
    start = closest(genome, gene[:30], N_GENE[0] - 1)
    end = closest(genome, gene[-30:], N_GENE[1] - 30) + 30
    return start, end


# =============================================================================
# PART 1: OUR SITES vs NCBI BLAST
# =============================================================================

def blast_best_site(hits, name, oligo, genome):
    """The hit with the fewest whole-oligo mismatches: (start 1-based, strand,
    whole-oligo mismatches, mismatches BLAST shows, aligned length)."""
    best = None
    for h in (h for h in hits if h["query"] == name):
        if h["s_start"] < h["s_end"]:
            strand, start = "+", h["s_start"] - (h["q_start"] - 1)
            pattern = oligo
        else:
            strand, start = "-", h["s_end"] - (len(oligo) - h["q_end"])
            pattern = reverse_complement(oligo)
        site = genome[start - 1:start - 1 + len(oligo)]
        if len(site) < len(oligo):
            continue
        mismatches = sum(not base_matches(g, p) for g, p in zip(site, pattern))
        candidate = (start, strand, mismatches, h["mismatches"], h["aligned"])
        if best is None or mismatches < best[2]:
            best = candidate
    return best


def part1(genomes, assays):
    print("PART 1: CDC oligos on each genome: our model vs NCBI BLAST")
    print("  (mismatches over the whole oligo; 'from 3'' = position counted from the 3' end)\n")
    for genome_name in GENOMES:
        genome = genomes[genome_name]
        hits = read_hits(os.path.join(SAVED, f"cdc_vs_{genome_name}.tsv"))
        print(f"{genome_name.upper()}")
        for assay, oligos in assays.items():
            result = run_temperature_model(oligos["forward"], oligos["reverse"], genome,
                                           probe=oligos["probe"], max_mismatches=2,
                                           conditions=CONDITIONS)
            print(f"  {assay}: assay at 55 C = {result['assay_at_ta']}")
            for role, seq in oligos.items():
                sites = [s for s in result["sites"] if s["primer"] == role]
                ours = min(sites, key=lambda s: s["mismatches"]) if sites else None
                blast = blast_best_site(hits, f"CDC_{assay}_{role}", seq, genome)
                if ours:
                    t = ours["temperature"]
                    tm = "" if t["tm"] is None else f", dTm {t['dTm']:+.1f} C" if t["dTm"] is not None else ""
                    ours_text = (f"{ours['start'] + 1}{ours['strand']} {ours['mismatches']} mm"
                                 + (f" (from 3': {ours['mismatch_from_3prime']})" if ours["mismatches"] else "")
                                 + tm)
                else:
                    ours_text = "no site"
                if blast:
                    shown = "" if blast[2] == blast[3] and blast[4] == len(seq) else \
                        f" [BLAST shows {blast[3]} mm over {blast[4]}/{len(seq)}]"
                    blast_text = f"{blast[0]}{blast[1]} {blast[2]} mm{shown}"
                else:
                    blast_text = "no hit"
                same = ours and blast and ours["start"] + 1 == blast[0] and ours["mismatches"] == blast[2]
                print(f"    {role:<8} ours {ours_text:<44} BLAST {blast_text:<40} {'same' if same else 'DIFFERENT'}")
        print()


# =============================================================================
# PART 2: THE THREE WEBSITE VERSIONS vs THE CDC DESIGNS
# =============================================================================

VERSIONS = {
    "Version 1": {"annealing_c": None},
    "Version 2": {"annealing_c": 55.0, "structure_checks": False},
    "Version 3": {"annealing_c": 55.0},
}


def cdc_products(genome, assays):
    """Where each CDC product sits on this genome (0-based start, end), from the first model."""
    products = {}
    for assay, oligos in assays.items():
        result = run_temperature_model(oligos["forward"], oligos["reverse"], genome,
                                       probe=oligos["probe"], max_mismatches=2)
        if result["products"]:
            p = result["products"][0]
            products[assay] = (p["start"], p["end"])
    return products


def part2(genomes, assays):
    print("PART 2: the three website versions on each genome's N gene, vs CDC N1/N2/N3")
    for genome_name in GENOMES:
        genome = genomes[genome_name]
        start, end = n_gene(genome, genomes["reference"])
        products = cdc_products(genome, assays)
        print(f"\n{genome_name.upper()}: N gene {start + 1}-{end}; CDC products "
              + ", ".join(f"{a} {s + 1}-{e}" for a, (s, e) in products.items()))
        for assay, oligos in assays.items():
            s, e = products[assay]
            amplicon = genome[s:e]
            plain = assay_temperature_check(oligos["forward"], oligos["reverse"], oligos["probe"],
                                            amplicon, conditions=QPCR_CONDITIONS, dimers=False)
            full = assay_temperature_check(oligos["forward"], oligos["reverse"], oligos["probe"],
                                           amplicon, conditions=QPCR_CONDITIONS)
            print(f"  CDC {assay} vs our rules: Version 2 {'pass' if plain['ok'] else 'FAIL ' + str(plain['problems'])}; "
                  f"Version 3 {'pass' if full['ok'] else 'FAIL ' + str(full['problems'])}")
        for version, options in VERSIONS.items():
            result = design_primers(genome, genome[start:end], 2, 60, 150, **options)
            print(f"  {version}: {result['possible']} possible designs")
            for n, d in enumerate(result["designs"], 1):
                overlaps = [a for a, (s, e) in products.items()
                            if d["product_start"] < e and s < d["product_end"]]
                works = [other for other in GENOMES
                         if run_temperature_model(d["forward"]["seq"], d["reverse"]["seq"], genomes[other],
                                                  probe=d["probe"]["seq"], max_mismatches=2)["assay_at_ta"]
                         == "product"]
                print(f"    design {n} {d['product_start'] + 1}-{d['product_end']} ({d['product_length']} bp)"
                      f"  F {d['forward']['seq']}  R {d['reverse']['seq']}"
                      f"  overlaps CDC: {', '.join(overlaps) if overlaps else 'none'}"
                      f"  works on: {'all four' if len(works) == len(GENOMES) else ', '.join(works)}")


if __name__ == "__main__":
    genomes = load_genomes()
    assays = cdc_assays()
    part1(genomes, assays)
    part2(genomes, assays)
