# =============================================================================
# CHECK: does the website's JavaScript designer agree with primer_design.py?
# =============================================================================
# There are three websites, each with a JavaScript translation of the designer:
#   website/              Version 1, GC only: design_primers(annealing_c=None)
#   website_temperature/  Version 2, temperature: annealing_c=55,
#                         structure_checks=False
#   website_structures/   Version 3, temperature + primer dimers + hairpins:
#                         annealing_c=55, structure_checks=True
# This script runs the Python and each JavaScript version on the same inputs
# and checks they return the same designs (sequences, positions, product,
# penalty, and for the temperature site each oligo's Tm to 1e-9 C), the
# same counts and the same warnings or error.
#
# Inputs: the fake genomes from test_primer_design.py, plus the real
# reference genome (N gene, and the CDC N1 region).
#
# JavaScript is run with macOS's built-in osascript, so no Node is needed.
#
# Run it with:  python3 journal/check_website_port.py
# =============================================================================

import glob
import json
import os
import subprocess
import tempfile

from first_model_factor_specificity import read_fasta, reverse_complement
from primer_design import design_primers
from second_model_annealing_temp import QPCR_CONDITIONS
from test_primer_design import GENOME, SHIFT, TARGET, basic_genome, random_bases

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SITES = {   # name: (JavaScript file, annealing temperature, dimer and hairpin checks)
    "Version 1: website (GC only)": (os.path.join(PROJECT, "website", "primer_design.js"), None, False),
    "Version 2: website_temperature": (os.path.join(PROJECT, "website_temperature", "primer_design.js"), 55.0, False),
    "Version 3: website_structures": (os.path.join(PROJECT, "website_structures", "primer_design.js"), 55.0, True),
}


# =============================================================================
# SECTION 1: THE CASES (genome, target, allowance, product range)
# =============================================================================

def cases():
    reference = read_fasta(glob.glob(os.path.join(PROJECT, "data", "reference", "*"))[0])[1]
    n_gene = reference[28274 - 1:29533]
    near_copy = "".join(SHIFT[b] if i % 12 == 0 else b for i, b in enumerate(TARGET))
    with_n = TARGET[:180] + "N" * 30 + TARGET[210:]
    island = TARGET[:100] + random_bases(40, seed=6, gc=1.0) + TARGET[140:]
    medium = random_bases(400, seed=11, gc=0.65)
    low = random_bases(400, seed=12, gc=0.33)
    return [
        ("basic", GENOME, TARGET, 2, (60, 150)),
        ("reverse-complement target", GENOME, reverse_complement(TARGET), 2, (60, 150)),
        ("N block", basic_genome(with_n), with_n, 2, (60, 150)),
        ("GC island", basic_genome(island), island, 2, (60, 150)),
        ("two copies", basic_genome(TARGET) + TARGET, TARGET, 2, (60, 150)),
        ("near copy, allowance 2", basic_genome(TARGET) + near_copy, TARGET, 2, (60, 150)),
        ("near copy, allowance 0", basic_genome(TARGET) + near_copy, TARGET, 0, (60, 150)),
        ("medium GC target", basic_genome(medium), medium, 2, (60, 150)),
        ("low GC target", basic_genome(low), low, 2, (60, 150)),
        ("no 40-60 % stretch", basic_genome("ATTA" * 50), "ATTA" * 50, 2, (60, 150)),
        ("product 54-80", GENOME, TARGET, 2, (54, 80)),
        ("product 250-350", GENOME, TARGET, 2, (250, 350)),
        ("allowance 1", GENOME, TARGET, 1, (60, 150)),
        ("error: GC too high", GENOME, random_bases(400, seed=9, gc=0.8), 2, (60, 150)),
        ("error: not in genome", GENOME, random_bases(400, seed=14, gc=0.45), 2, (60, 150)),
        ("error: shortest product 53", GENOME, TARGET, 2, (53, 150)),
        ("N gene, reference", reference, n_gene, 2, (60, 150)),
        ("N gene, reference, 54-80 bp", reference, n_gene, 2, (54, 80)),
        ("CDC N1 region, reference", reference, reference[28236:28408], 2, (54, 80)),
        ("error: backwards range", GENOME, TARGET, 2, (150, 60)),
    ] + web_cases()


def web_cases():
    """The five cases in test_cases_web_specificity.txt, read from the website's
    example_data.js, each at allowance 0, 1 and 2."""
    with open(os.path.join(PROJECT, "website", "example_data.js")) as handle:
        line = next(l for l in handle if l.startswith("window.TEST_EXAMPLES = "))
    tests = json.loads(line[len("window.TEST_EXAMPLES = "):].rstrip().rstrip(";"))
    return [(f"{t['label'][:6]}, allowance {a}", t["genome"], t["target"], a, (60, 150))
            for t in tests for a in (0, 1, 2)]


# =============================================================================
# SECTION 2: RUN BOTH AND COMPARE
# =============================================================================

def summary(result):
    """The parts of a result both versions must agree on."""
    keep = ("forward", "reverse", "probe")
    tms = [None if d.get("tm") is None else
           {k: round(v, 9) for k, v in d["tm"].items()} for d in result["designs"]]
    dimers = [None if not d.get("dimers") else
              [[x["pair"], round(x["any"]["dG"], 6), x["any"]["top"], x["any"]["bottom"],
                None if x["three_prime"] is None else round(x["three_prime"]["dG"], 6)]
               for x in d["dimers"]] for d in result["designs"]]
    hairpins = [None if not d.get("hairpins") else
                {k: None if f is None else [round(f["dG"], 6), f["stem_5"], f["loop"], f["three_prime"]]
                 for k, f in d["hairpins"].items()} for d in result["designs"]]
    return {
        "tms": tms,
        "dimers": dimers,
        "hairpins": hairpins,
        "dimer_rejected": result.get("dimer_rejected", 0),
        "hairpin_rejected": result.get("hairpin_rejected", 0),
        "warnings": result["warnings"],
        "target": result["target"],
        "possible": result["possible"],
        "checked": result["checked"],
        "designs": [{**{o: [d[o]["seq"], d[o]["start"], d[o]["end"]] for o in keep},
                     "product": [d["product_start"], d["product_end"], d["product_length"]],
                     "penalty": d["match_penalty"]} for d in result["designs"]],
    }


def run_python(genome, target, allowance, product, annealing, conditions=None, structures=True):
    try:
        return summary(design_primers(genome, target, allowance, *product, annealing_c=annealing,
                                      conditions=conditions or QPCR_CONDITIONS,
                                      structure_checks=structures))
    except ValueError as error:
        return {"error": str(error)}


def run_javascript(js_file, all_cases, structures=True):
    inputs = [{"genome": c[1], "target": c[2], "allowance": c[3], "product": list(c[4]),
               "annealing": c[5], "conditions": c[6] if len(c) > 6 else None,
               "structures": structures} for c in all_cases]
    runner = """
    function round9(x) { return Math.round(x * 1e9) / 1e9; }
    function round6(x) { return Math.round(x * 1e6) / 1e6; }
    var out = INPUTS.map(function (c) {
      try {
        var r = c.annealing === "gc-site"
          ? PrimerDesign.designPrimers(c.genome, c.target, c.allowance, c.product[0], c.product[1])
          : PrimerDesign.designPrimers(c.genome, c.target, c.allowance, c.product[0], c.product[1], 3, c.annealing,
                                       c.conditions || PrimerDesign.QPCR_CONDITIONS, c.structures);
        return {
          dimers: r.designs.map(function (d) {
            if (!d.dimers) return null;
            return d.dimers.map(function (x) {
              return [x.pair, round6(x.any.dG), x.any.top, x.any.bottom,
                      x.three_prime === null ? null : round6(x.three_prime.dG)];
            });
          }),
          dimer_rejected: r.dimer_rejected,
          hairpin_rejected: r.hairpin_rejected,
          hairpins: r.designs.map(function (d) {
            if (!d.hairpins) return null;
            var h = {};
            for (var k in d.hairpins) {
              var f = d.hairpins[k];
              h[k] = f === null ? null : [round6(f.dG), f.stem_5, f.loop, f.three_prime];
            }
            return h;
          }),
          tms: r.designs.map(function (d) {
            if (!d.tm) return null;
            var t = {};
            for (var k in d.tm) t[k] = round9(d.tm[k]);
            return t;
          }),
          warnings: r.warnings, target: r.target, possible: r.possible, checked: r.checked,
          designs: r.designs.map(function (d) {
            return {
              forward: [d.forward.seq, d.forward.start, d.forward.end],
              reverse: [d.reverse.seq, d.reverse.start, d.reverse.end],
              probe: [d.probe.seq, d.probe.start, d.probe.end],
              product: [d.product_start, d.product_end, d.product_length],
              penalty: d.match_penalty,
            };
          }),
        };
      } catch (e) { return { error: e.message }; }
    });
    JSON.stringify(out);
    """
    with open(js_file) as handle:
        script = handle.read() + "\nvar INPUTS = " + json.dumps(inputs) + ";\n" + runner
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False) as handle:
        handle.write(script)
        path = handle.name
    try:
        output = subprocess.run(["osascript", "-l", "JavaScript", path],
                                capture_output=True, text=True, check=True).stdout
    finally:
        os.remove(path)
    return json.loads(output)


def site_cases(site_annealing):
    """Each case with its annealing temperature: the GC-only site uses None
    (its designPrimers has no annealing argument: marked "gc-site"); the
    temperature site runs every case at 55 C, plus a few at 60 C and without."""
    if site_annealing is None:
        return [(name, g, t, a, p, "gc-site") for name, g, t, a, p in cases()]
    base = [(name, g, t, a, p, site_annealing) for name, g, t, a, p in cases()]
    extra = [(f"{name}, {label}", g, t, a, p, ta) for name, g, t, a, p in cases()
             if name in ("basic", "N gene, reference", "medium GC target", "low GC target")
             for label, ta in (("annealing 60 C", 60.0), ("no temperature", None))]
    errors = [("error: annealing 40 C", GENOME, TARGET, 2, (60, 150), 40.0),
              ("error: annealing 80 C", GENOME, TARGET, 2, (60, 150), 80.0)]
    mixes = [("5 mM Mg2+, 250 nM primer", dict(QPCR_CONDITIONS, divalent_mM=5.0, primer_nM=250.0)),
             ("100 mM K+, 250 nM probe", dict(QPCR_CONDITIONS, monovalent_mM=100.0, probe_nM=250.0)),
             ("error: 20 mM Mg2+", dict(QPCR_CONDITIONS, divalent_mM=20.0))]
    conditions = [(f"{name}, {label}", g, t, a, p, site_annealing, mix) for name, g, t, a, p in cases()
                  if name in ("basic", "N gene, reference") for label, mix in mixes]
    return base + extra + errors + conditions


if __name__ == "__main__":
    total_failures = 0
    for site, (js_file, site_annealing, structures) in SITES.items():
        print(f"\n=== {site} ===")
        all_cases = site_cases(site_annealing)
        javascript = run_javascript(js_file, all_cases, structures)
        failures = 0
        for case, js in zip(all_cases, javascript):
            name, genome, target, allowance, product, annealing = case[:6]
            py = run_python(genome, target, allowance, product,
                            None if annealing == "gc-site" else annealing,
                            case[6] if len(case) > 6 else None, structures)
            if annealing == "gc-site" and "error" not in js:
                js = {**js, "tms": [None] * len(js.get("designs", [])),
                      "dimers": [None] * len(js.get("designs", [])), "dimer_rejected": 0,
                      "hairpins": [None] * len(js.get("designs", [])), "hairpin_rejected": 0}
            same = py == js
            failures += not same
            if "error" in py:
                detail = "error: " + py["error"][:50]
            else:
                detail = f"{len(py['designs'])} designs, {py['possible']} possible"
            print(f"  {'PASS' if same else 'FAIL'}: {name:40} {detail}")
            if not same:
                print("    Python:    ", json.dumps(py)[:400])
                print("    JavaScript:", json.dumps(js)[:400])
        print(f"{len(all_cases) - failures} of {len(all_cases)} cases give the same result in both.")
        total_failures += failures
    print("\nALL SAME" if not total_failures else f"\n{total_failures} DIFFERENCES")
