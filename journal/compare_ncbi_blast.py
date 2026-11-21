# =============================================================================
# COMPARE WITH NCBI BLAST: the website's five test cases
# =============================================================================
# Do NCBI BLAST and our first model agree on where each designed primer and
# probe binds?
#
# For every oligo the website designs in test_cases_web_specificity.txt
# (cases 1-4; case 5 stops at the GC check, so it has no oligos), plus the
# top-ranked candidate the designer had to REJECT in cases 3 and 4:
#
#   1. OUR SITES: every site on the case's genome, on either strand, with
#      0-3 mismatches (a plain scan of every position, no shortcuts).
#   2. NCBI "ALIGN TWO SEQUENCES" (BlastAlign.cgi): the case's oligos
#      against the case's genome. blastn, word size 7 (the smallest NCBI
#      allows, as for short queries), expect 1000, no low-complexity filter.
#   3. NCBI core_nt: every oligo against NCBI's whole nucleotide collection,
#      to see whether the made-up sequences match anything real.
#      blastn, word size 11, expect 30000 (Primer-BLAST's value), top 50 hits.
#
# Asking NCBI takes minutes, so the answers are saved in data/ncbi_blast/ and
# this script reads them from there. Run with --ask-ncbi to search again
# (it waits 15 s between submissions and 60 s between status checks, as
# NCBI's usage rules ask).
#
# Run it with:  python3 journal/compare_ncbi_blast.py [--ask-ncbi]
# =============================================================================

import json
import os
import re
import subprocess
import sys
import time

from first_model_factor_specificity import reverse_complement
from primer_design import clean_sequence, design_primers, find_target, possible_designs

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SAVED = os.path.join(PROJECT, "data", "ncbi_blast")
BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
BLAST_TWO = "https://blast.ncbi.nlm.nih.gov/BlastAlign.cgi"   # two-sequence searches


# =============================================================================
# SECTION 1: THE OLIGOS AND OUR OWN SITES
# =============================================================================

def test_cases():
    with open(os.path.join(PROJECT, "website", "example_data.js")) as handle:
        line = next(l for l in handle if l.startswith("window.TEST_EXAMPLES = "))
    return json.loads(line[len("window.TEST_EXAMPLES = "):].rstrip().rstrip(";"))


def our_sites(oligo, genome, most=3):
    """Every site with <= most mismatches: (1-based start, strand, mismatches)."""
    sites = []
    for strand, pattern in (("+", oligo), ("-", reverse_complement(oligo))):
        for i in range(len(genome) - len(pattern) + 1):
            mismatches = sum(a != b for a, b in zip(genome[i:i + len(pattern)], pattern))
            if mismatches <= most:
                sites.append((i + 1, strand, mismatches))
    return sorted(sites)


def collect_oligos(tests):
    """One row per oligo (per case): its sequence, role, and which designs use it."""
    rows = []
    for number, test in enumerate(tests, 1):
        genome, target = clean_sequence(test["genome"]), clean_sequence(test["target"])
        uses = {}
        for allowance in (0, 1, 2):
            try:
                result = design_primers(genome, target, allowance)
            except ValueError:
                continue            # case 5: GC too high, nothing designed
            picked = [(f"allowance {allowance} design {k}", d)
                      for k, d in enumerate(result["designs"], 1)]
            if not picked and result["possible"]:
                start, end, _ = find_target(genome, target)
                picked = [(f"allowance {allowance} top candidate (rejected)",
                           possible_designs(genome, start, end)[0])]
            for label, d in picked:
                for role in ("forward", "reverse", "probe"):
                    uses.setdefault((role, d[role]["seq"]), []).append(label)
        for (role, seq), labels in uses.items():
            rows.append({"case": number, "role": role, "seq": seq, "uses": labels,
                         "ours": our_sites(seq, genome),
                         "name": f"c{number}_{role[0].upper()}{len(rows)}"})
    return rows


# =============================================================================
# SECTION 2: ASKING NCBI (only with --ask-ncbi)
# -----------------------------------------------------------------------------
# curl is used rather than Python's urllib, which rejects some networks'
# certificates. No email or personal details are sent.
# =============================================================================

def ncbi(url, params):
    args = ["curl", "-sS", "--max-time", "120", url]
    for key, value in dict(params, TOOL="compare_ncbi_blast").items():
        args += ["--data-urlencode", f"{key}={value}"]
    return subprocess.run(args, capture_output=True, text=True, check=True).stdout


def submit(url, params):
    page = ncbi(url, dict(params, CMD="Put"))
    rid = re.search(r"RID = (\S+)", page)
    if not rid:
        sys.exit("NCBI did not accept the search:\n" + page[:2000])
    time.sleep(15)
    return rid.group(1)


def fetch(rid):
    while True:
        info = ncbi(BLAST, {"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"})
        status = re.search(r"Status=(\w+)", info).group(1)
        if status == "READY":
            return ncbi(BLAST, {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "Text",
                                "ALIGNMENT_VIEW": "Tabular", "HITLIST_SIZE": "50"})
        if status != "WAITING":
            sys.exit(f"NCBI search {rid} ended with status {status}")
        time.sleep(60)


def ask_ncbi(tests, rows):
    os.makedirs(SAVED, exist_ok=True)
    rids = {}
    for number in sorted({r["case"] for r in rows}):
        query = "".join(f">{r['name']}\n{r['seq']}\n" for r in rows if r["case"] == number)
        rids[f"case{number}_genome"] = submit(BLAST_TWO, {
            "PROGRAM": "blastn", "QUERY": query,
            "SUBJECTS": f">case{number}_genome\n" + clean_sequence(tests[number - 1]["genome"]),
            "WORD_SIZE": "7", "EXPECT": "1000", "FILTER": "F", "DUST": "no"})
    unique = {}
    for r in rows:
        unique.setdefault(r["seq"], r["name"])
    rids["core_nt"] = submit(BLAST, {
        "PROGRAM": "blastn", "DATABASE": "core_nt",
        "QUERY": "".join(f">{name}\n{seq}\n" for seq, name in unique.items()),
        "WORD_SIZE": "11", "EXPECT": "30000", "FILTER": "F", "HITLIST_SIZE": "50"})
    for name, rid in rids.items():
        print(f"waiting for NCBI search {rid} ({name}) ...", flush=True)
        with open(os.path.join(SAVED, f"{name}.tsv"), "w") as handle:
            handle.write(f"# NCBI RID {rid}, {time.strftime('%Y-%m-%d')}\n" + fetch(rid))


# =============================================================================
# SECTION 3: READING NCBI'S TABLES
# -----------------------------------------------------------------------------
# BLAST reports a LOCAL alignment: it may leave off the oligo's ends, and a
# base it leaves off is not counted as a mismatch. For the case genomes we
# have the sequence, so the whole oligo is lined up at the hit's site and
# its mismatches counted exactly. For core_nt we do not, so the best we can
# say is "at most  mismatches + gaps + bases left off".
# =============================================================================

def read_hits(path):
    hits = []
    seen = set()
    with open(path) as handle:
        for line in handle:
            line = re.sub(r"<[^>]+>", "", line).strip()
            if not line or line.startswith("#") or "\t" not in line:
                continue
            f = line.split("\t")
            if tuple(f) in seen:
                continue            # the same query sent twice
            seen.add(tuple(f))
            hits.append({"query": f[0], "subject": f[1], "identity": float(f[2]),
                         "aligned": int(f[3]), "mismatches": int(f[4]), "gaps": int(f[5]),
                         "q_start": int(f[6]), "q_end": int(f[7]),
                         "s_start": int(f[8]), "s_end": int(f[9]), "evalue": float(f[10])})
    return hits


def whole_oligo_mismatches(hit, length):
    """Upper bound: every base BLAST left off counted as a mismatch."""
    return hit["mismatches"] + hit["gaps"] + (hit["q_start"] - 1) + (length - hit["q_end"])


def whole_oligo_site(hit, oligo, genome):
    """The hit's site extended to the whole oligo: (1-based start, strand, mismatches)."""
    length = len(oligo)
    if hit["s_start"] < hit["s_end"]:
        strand, start = "+", hit["s_start"] - (hit["q_start"] - 1)
        pattern = oligo
    else:
        strand, start = "-", hit["s_end"] - (length - hit["q_end"])
        pattern = reverse_complement(oligo)
    if hit["gaps"] or start < 1 or start - 1 + length > len(genome):
        return start, strand, whole_oligo_mismatches(hit, length)
    site = genome[start - 1:start - 1 + length]
    return start, strand, sum(a != b for a, b in zip(site, pattern))


# =============================================================================
# SECTION 4: THE COMPARISON
# =============================================================================

def compare_with_genomes(rows, tests):
    print("PART 1: each case's oligos against the case's own genome")
    print("(a 'site' = the whole oligo lines up with <= 3 mismatches; BLAST's "
          "short partial hits of 7-16 bases are counted separately)\n")
    for number in sorted({r["case"] for r in rows}):
        hits = read_hits(os.path.join(SAVED, f"case{number}_genome.tsv"))
        genome = clean_sequence(tests[number - 1]["genome"])
        print(f"CASE {number}")
        for r in (r for r in rows if r["case"] == number):
            length = len(r["seq"])
            blast_sites, partial = [], 0
            for h in (h for h in hits if h["query"] == r["name"]):
                start, strand, total = whole_oligo_site(h, r["seq"], genome)
                if total <= 3:
                    blast_sites.append((start, strand, total, h["mismatches"], h["aligned"]))
                else:
                    partial += 1
            agree = sorted(s[:3] for s in blast_sites) == r["ours"]
            print(f"  {r['role']:<8}{r['seq']:<27}{'; '.join(r['uses'])}")
            print(f"    ours : {[f'{p}{s} {m}mm' for p, s, m in r['ours']]}")
            print(f"    BLAST: {[f'{p}{s} {m}mm (BLAST shows {shown}mm over {al}/{length} bases)' if al < length or shown != m else f'{p}{s} {m}mm' for p, s, m, shown, al in sorted(blast_sites)]}"
                  f"{f' + {partial} short partial hit(s)' if partial else ''}")
            if not agree:
                missed = [s for s in r["ours"] if s not in [b[:3] for b in blast_sites]]
                print(f"    -> BLAST MISSED {missed}")
        print()


def compare_with_core_nt(rows):
    path = os.path.join(SAVED, "core_nt.tsv")
    if not os.path.exists(path):
        print("PART 2: no saved core_nt search yet (run with --ask-ncbi)")
        return
    hits = read_hits(path)
    print("PART 2: every oligo against NCBI core_nt (all of GenBank's main nucleotides)")
    print("  hits grouped by AT MOST how many mismatches over the whole oligo (bases")
    print("  BLAST left off count as mismatches); NCBI returns the top 50 subjects only\n")
    lengths = {r["name"]: len(r["seq"]) for r in rows}
    for name, length in lengths.items():
        mine = [h for h in hits if h["query"] == name]
        counts = {}
        for h in mine:
            m = whole_oligo_mismatches(h, length)
            counts[m] = counts.get(m, 0) + 1
        best = min(mine, key=lambda h: whole_oligo_mismatches(h, length), default=None)
        if not mine:
            continue
        print(f"  {name:<8}{length} nt  hits: {len(mine):>3}  "
              f"by <=mismatches: {dict(sorted(counts.items()))}  "
              f"best: {best['subject']} (<={whole_oligo_mismatches(best, length)} mm, E={best['evalue']:g})")


# =============================================================================
# SECTION 5: THE N GENE AND THE CDC PRIMERS (Primer-BLAST and BLAST)
# -----------------------------------------------------------------------------
# Saved from NCBI on 2026-11-19 (searched by hand from this script's helpers;
# the job keys are at the top of each file):
#   primer_blast_n_gene_design.txt  Primer-BLAST designing on the N gene
#       (28274-29533) with our settings: primers 18-24 nt, probe 18-30 nt,
#       40-60 % GC, product 60-150 bp; its own Tm defaults (57-63 C, opt 60,
#       pair within 3 C); specificity vs the human reference genome.
#   primer_blast_<pair>.txt  our 3 N gene designs and CDC N1-N3 given as
#       fixed pairs: Primer-BLAST's Tm, and specificity vs human.
#   n_gene_core_nt_no_sarscov2.json  BLAST of the same oligos against
#       core_nt with SARS-CoV-2 (txid2697049) excluded.
# Primer-BLAST's default conditions are 50 mM Na+, 1.5 mM Mg2+, 0.6 mM dNTP,
# 50 nM oligo; ours (CONDITIONS) are qPCR's 3 mM Mg2+, 0.8 mM dNTP, 250 nM.
# =============================================================================

PRIMER_BLAST_CONDITIONS = {"monovalent_mM": 50.0, "divalent_mM": 1.5, "dntp_mM": 0.6, "oligo_nM": 50.0}
PAIRS = ("ours1", "ours2", "ours3", "CDC_N1", "CDC_N2", "CDC_N3")


def primer_blast_designs(path):
    """Each 'Primer pair' block: {role: (seq, start, end, tm, gc)} + product length."""
    text = open(path).read()
    designs = []
    for block in re.split(r"\nPrimer pair \d+\n", text)[1:]:
        design = {}
        for role, label in (("forward", "Forward primer"), ("reverse", "Reverse primer"),
                            ("probe", "Internal oligo")):
            m = re.search(label + r"\t([ACGT]+)\t\w+\t(\d+)\t(\d+)\t(\d+)\t([\d.]+)\t([\d.]+)", block)
            if m:
                design[role] = (m.group(1), int(m.group(3)), int(m.group(4)),
                                float(m.group(5)), float(m.group(6)))
        design["product"] = int(re.search(r"Product length\t(\d+)", block).group(1))
        designs.append(design)
    return designs


def compare_n_gene():
    from first_model_factor_specificity import read_fasta
    from second_model_annealing_temp import CONDITIONS, perfect_tm
    import glob
    if not os.path.exists(os.path.join(SAVED, "primer_blast_n_gene_design.txt")):
        print("PART 3: no saved Primer-BLAST results")
        return
    print("\nPART 3: THE N GENE: OUR DESIGNER vs NCBI PRIMER-BLAST\n")

    print("3a. Primer Tm: Primer-BLAST vs our second model (perfect_tm)")
    print(f"  {'pair':<8}{'oligo':<26}{'Primer-BLAST':>13}{'ours, its cond.':>16}{'ours, qPCR':>12}")
    for pair in PAIRS:
        d = primer_blast_designs(os.path.join(SAVED, f"primer_blast_{pair}.txt"))[0]
        tms = []
        for role in ("forward", "reverse"):
            seq, *_, tm, _ = d[role]
            tms.append(tm)
            print(f"  {pair:<8}{seq:<26}{tm:>13.2f}{perfect_tm(seq, PRIMER_BLAST_CONDITIONS):>16.2f}"
                  f"{perfect_tm(seq, CONDITIONS):>12.2f}")
        print(f"  {'':<8}{'difference between the primers':<34}{abs(tms[0] - tms[1]):>5.2f}"
              f"{'  <- over Primer-BLAST limit of 3 C' if abs(tms[0] - tms[1]) > 3 else ''}"
              f"{'  <- below its 57 C minimum' if min(tms) < 57 else ''}")

    print("\n3b. Primer-BLAST's 10 N gene designs, and where our designer ranks the same product")
    reference = read_fasta(glob.glob(os.path.join(PROJECT, "data", "reference", "*"))[0])[1]
    gene_start = 28274 - 1
    start, end, _ = find_target(reference, reference[gene_start:29533])
    ours = possible_designs(reference, start, end)
    rank = {(d["product_start"], d["product_end"]): i for i, d in reversed(list(enumerate(ours)))}
    for n, d in enumerate(primer_blast_designs(os.path.join(SAVED, "primer_blast_n_gene_design.txt")), 1):
        f, r, p = d["forward"], d["reverse"], d.get("probe")
        key = (gene_start + f[1] - 1, gene_start + r[1])
        i = rank.get(key)
        where = (f"our rank {i + 1:>6} of {len(ours)} (penalty {ours[i]['match_penalty']})"
                 if i is not None else "not in our list")
        print(f"  {n:>2}  F {f[0]} {f[3]:.1f}C  R {r[0]} {r[3]:.1f}C  "
              f"P {p[0] if p else '-'} {p[3] if p else '':}  {d['product']} bp  "
              f"genome {key[0] + 1}-{key[1]}  {where}")

    path = os.path.join(SAVED, "n_gene_core_nt_no_sarscov2.json")
    print("\n3c. BLAST vs core_nt, SARS-CoV-2 excluded: top organisms with <= 2 mismatches")
    for report in json.load(open(path))["BlastOutput2"]:
        search = report["report"]["results"]["search"]
        length, counts = search["query_len"], {}
        for hit in search.get("hits", []):
            name = hit["description"][0].get("sciname", "?")
            if any(length - hsp["identity"] <= 2 for hsp in hit["hsps"]):
                counts[name] = counts.get(name, 0) + 1
        top = sorted(counts.items(), key=lambda kv: -kv[1])[:3]
        print(f"  {search['query_title']:<16}{sum(counts.values()):>3} hits  {top}")


if __name__ == "__main__":
    tests = test_cases()
    rows = collect_oligos(tests)
    if "--ask-ncbi" in sys.argv:
        ask_ncbi(tests, rows)
    compare_with_genomes(rows, tests)
    compare_with_core_nt(rows)
    compare_n_gene()
