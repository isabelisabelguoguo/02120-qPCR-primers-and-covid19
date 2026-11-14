# =============================================================================
# CHECK: does the website's JavaScript designer agree with primer_design.py?
# =============================================================================
# The website (website/index.html) runs website/primer_design.js, a
# JavaScript translation of primer_design.py and the first-model functions it
# uses. This script runs BOTH on the same inputs and checks they return the
# same designs (sequences, positions, product, penalty), the same counts and
# the same warnings or error.
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
from test_primer_design import GENOME, SHIFT, TARGET, basic_genome, random_bases

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JS_FILE = os.path.join(PROJECT, "website", "primer_design.js")


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
    return {
        "warnings": result["warnings"],
        "target": result["target"],
        "possible": result["possible"],
        "checked": result["checked"],
        "designs": [{**{o: [d[o]["seq"], d[o]["start"], d[o]["end"]] for o in keep},
                     "product": [d["product_start"], d["product_end"], d["product_length"]],
                     "penalty": d["match_penalty"]} for d in result["designs"]],
    }


def run_python(genome, target, allowance, product):
    try:
        return summary(design_primers(genome, target, allowance, *product))
    except ValueError as error:
        return {"error": str(error)}


def run_javascript(all_cases):
    inputs = [{"genome": g, "target": t, "allowance": a, "product": list(p)} for _, g, t, a, p in all_cases]
    runner = """
    var out = INPUTS.map(function (c) {
      try {
        var r = PrimerDesign.designPrimers(c.genome, c.target, c.allowance, c.product[0], c.product[1]);
        return {
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
    with open(JS_FILE) as handle:
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


if __name__ == "__main__":
    all_cases = cases()
    javascript = run_javascript(all_cases)
    failures = 0
    for (name, genome, target, allowance, product), js in zip(all_cases, javascript):
        py = run_python(genome, target, allowance, product)
        same = py == js
        failures += not same
        if "error" in py:
            detail = "error: " + py["error"][:50]
        else:
            detail = f"{len(py['designs'])} designs, {py['possible']} possible"
        print(f"  {'PASS' if same else 'FAIL'}: {name:28} {detail}")
        if not same:
            print("    Python:    ", json.dumps(py)[:300])
            print("    JavaScript:", json.dumps(js)[:300])
    print(f"\n{len(all_cases) - failures} of {len(all_cases)} cases give the same result in both.")
