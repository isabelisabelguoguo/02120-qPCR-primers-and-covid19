// =============================================================================
// PRIMER DESIGN FOR THE WEBSITE: JavaScript version of journal/primer_design.py
// =============================================================================
// 02-120 project: Do early SARS-CoV-2 qPCR primers still match later variants?
//
// TEMPERATURE VERSION: website_temperature/ (website/ keeps the GC-only
// version 3 designer).
//
// A browser cannot run Python, so this file translates, function by function,
//   - journal/primer_design.py (version 4: with temperature), all of it,
//   - the perfect-match Tm, the qPCR temperature rules and the primer-dimer
//     scan of journal/second_model_annealing_temp.py (version 3, Sections
//     1-4, 9, 10), and
//   - the parts of journal/first_model_factor_specificity.py it uses:
//     the k-mer index, find_binding_sites, find_indel_sites,
//     pair_into_products, probe_inside, product_status and their helpers.
// Near-miss search and printing are not needed by the designer, so they are
// left out. Each JavaScript function names the Python function it copies.
//
// The Python files stay the reference. journal/check_website_port.py runs
// both on the same inputs and checks they give the same designs.
//
// Positions are 0-based inside, end exclusive, exactly as in Python.
// Everything is exported on one object, PrimerDesign.
// =============================================================================

(function (root) {
  "use strict";

  // ===========================================================================
  // FIRST MODEL, SECTION 1: DNA HELPERS
  // ===========================================================================

  const COMPLEMENT = {
    A: "T", T: "A", C: "G", G: "C",
    R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K",
    B: "V", V: "B", D: "H", H: "D", N: "N",
  };

  // which genome bases each primer letter accepts
  const IUPAC = {
    A: "A", C: "C", G: "G", T: "T",
    R: "AG", Y: "CT", S: "CG", W: "AT", K: "GT", M: "AC",
    B: "CGT", D: "AGT", H: "ACT", V: "ACG", N: "ACGT",
  };

  const KNOWN_BASES = "ACGT";


  // ===========================================================================
  // SECOND MODEL, SECTIONS 1-4 AND 9: PERFECT-MATCH Tm AND THE qPCR RULES
  // ---------------------------------------------------------------------------
  // Only what the designer needs: the Tm of an oligo on its perfectly
  // matching site, and the product Tm. Values from
  // data/thermodynamics/santalucia_1998_nn.tsv (SantaLucia 1998), 1 M NaCl,
  // dH kcal/mol, dS cal/(K*mol).
  // ===========================================================================

  const MATCHED = {
    "AA/TT": [-7.9, -22.2], "AT/TA": [-7.2, -20.4], "TA/AT": [-7.2, -21.3],
    "CA/GT": [-8.5, -22.7], "GT/CA": [-8.4, -22.4], "CT/GA": [-7.8, -21.0],
    "GA/CT": [-8.2, -22.2], "CG/GC": [-10.6, -27.2], "GC/CG": [-9.8, -24.4],
    "GG/CC": [-8.0, -19.9],
    init_w_term_GC: [0.1, -2.8], init_w_term_AT: [2.3, 4.1], symmetry_correction: [0.0, -1.4],
  };
  const R = 1.9872;   // gas constant, cal/(K*mol)

  // data/thermodynamics/allawi_santalucia_internal_mismatch_nn.tsv (51 stacks)
  const INTERNAL_MISMATCH = {"AG/TT": [1.0, 0.9], "AT/TG": [-2.5, -8.3], "CG/GT": [-4.1, -11.7], "CT/GG": [-2.8, -8.0], "GG/CT": [3.3, 10.4], "GG/TT": [5.8, 16.3], "GT/CG": [-4.4, -12.3], "GT/TG": [4.1, 9.5], "TG/AT": [-0.1, -1.7], "TG/GT": [-1.4, -6.2], "TT/AG": [-1.3, -5.3], "AA/TG": [-0.6, -2.3], "AG/TA": [-0.7, -2.3], "CA/GG": [-0.7, -2.3], "CG/GA": [-4.0, -13.2], "GA/CG": [-0.6, -1.0], "GG/CA": [0.5, 3.2], "TA/AG": [0.7, 0.7], "TG/AA": [3.0, 7.4], "AC/TT": [0.7, 0.2], "AT/TC": [-1.2, -6.2], "CC/GT": [-0.8, -4.5], "CT/GC": [-1.5, -6.1], "GC/CT": [2.3, 5.4], "GT/CC": [5.2, 13.5], "TC/AT": [1.2, 0.7], "TT/AC": [1.0, 0.7], "AA/TC": [2.3, 4.6], "AC/TA": [5.3, 14.6], "CA/GC": [1.9, 3.7], "CC/GA": [0.6, -0.6], "GA/CC": [5.2, 14.2], "GC/CA": [-0.7, -3.8], "TA/AC": [3.4, 8.0], "TC/AA": [7.6, 20.2], "AA/TA": [1.2, 1.7], "CA/GA": [-0.9, -4.2], "GA/CA": [-2.9, -9.8], "TA/AA": [4.7, 12.9], "AC/TC": [0.0, -4.4], "CC/GC": [-1.5, -7.2], "GC/CC": [3.6, 8.9], "TC/AC": [6.1, 16.4], "AG/TG": [-3.1, -9.5], "CG/GG": [-4.9, -15.3], "GG/CG": [-6.0, -15.8], "TG/AG": [1.6, 3.6], "AT/TT": [-2.7, -10.8], "CT/GT": [-5.0, -15.8], "GT/CT": [-2.2, -8.4], "TT/AT": [0.2, -1.5]};

  // the second model's user settings
  const ANNEALING_TEMP_C = 55.0;
  const CONDITIONS = { monovalent_mM: 50.0, divalent_mM: 3.0, dntp_mM: 0.8, oligo_nM: 250.0 };
  const DENATURATION_C = 95.0;
  const PRIMER_NM = 500.0;
  const PROBE_NM = 125.0;
  const PRIMER_TM_ABOVE_ANNEALING = [6.0, 13.0];
  const PRIMER_TM_MAX_DIFFERENCE = 4.0;
  const PROBE_ABOVE_PRIMERS_C = 2.0;
  const PRODUCT_BELOW_DENATURATION_C = 5.0;
  const QPCR_CONDITIONS = Object.assign({}, CONDITIONS, { primer_nM: PRIMER_NM, probe_nM: PROBE_NM });
  const CONDITION_LIMITS = {
    monovalent_mM: [10.0, 200.0], divalent_mM: [0.0, 10.0], dntp_mM: [0.0, 2.0],
    primer_nM: [50.0, 2000.0], probe_nM: [50.0, 1000.0],
  };
  const DIMER_TEMP_C = 37.0;
  const DIMER_MAX_DG = -9000.0;
  const DIMER_3END_MAX_DG = -5000.0;

  // stack_value
  function stackValue(top, bottom, table = MATCHED) {
    const key = `${top}/${bottom}`;
    if (key in table) return table[key];
    return table[key.split("").reverse().join("")];
  }

  // equivalent_sodium_molar
  function equivalentSodiumMolar(conditions) {
    const freeMg = conditions.divalent_mM - conditions.dntp_mM;
    const sodiumMM = conditions.monovalent_mM + (freeMg > 0 ? 120 * Math.sqrt(freeMg) : 0);
    return sodiumMM / 1000;
  }

  // conditions_for
  function conditionsFor(role, conditions = CONDITIONS) {
    const amount = role === "probe"
      ? ("probe_nM" in conditions ? conditions.probe_nM : PROBE_NM)
      : ("primer_nM" in conditions ? conditions.primer_nM : PRIMER_NM);
    return Object.assign({}, conditions, { oligo_nM: amount });
  }

  // check_conditions
  function checkConditions(conditions) {
    for (const name of Object.keys(CONDITION_LIMITS)) {
      const [low, high] = CONDITION_LIMITS[name];
      const value = conditions[name];
      if (typeof value !== "number" || value !== value || !(low <= value && value <= high)) {
        throw new Error(`${name.replace("_", " ")} must be ${low}-${high}`);
      }
    }
  }

  // salt_corrected_entropy
  function saltCorrectedEntropy(ds, length, conditions) {
    return ds + 0.368 * (length - 1) * Math.log(equivalentSodiumMolar(conditions));
  }

  // duplex_energy, for a duplex whose two END pairs are Watson-Crick (all the
  // dimer scan needs; the terminal-mismatch branch of the Python never runs
  // there). Returns [dH, dS] or null for two mismatches side by side.
  function duplexEnergy(top, bottom) {
    const n = top.length;
    const paired = [];
    for (let i = 0; i < n; i++) paired.push(bottom[i] === COMPLEMENT[top[i]]);
    let dh = 0.0, ds = 0.0;
    for (let i = 0; i < n - 1; i++) {
      if (!paired[i] && !paired[i + 1]) return null;
      const table = paired[i] && paired[i + 1] ? MATCHED : INTERNAL_MISMATCH;
      const value = stackValue(top.slice(i, i + 2), bottom.slice(i, i + 2), table);
      dh += value[0];
      ds += value[1];
    }
    for (const end of [top[0], top[n - 1]]) {
      const value = "GC".includes(end) ? MATCHED.init_w_term_GC : MATCHED.init_w_term_AT;
      dh += value[0];
      ds += value[1];
    }
    if (paired.every((p) => p) && top === reverseComplement(top)) {
      dh += MATCHED.symmetry_correction[0];
      ds += MATCHED.symmetry_correction[1];
    }
    return [dh, ds];
  }

  // dimer_scan
  function dimerScan(a, b, conditions, tempC = DIMER_TEMP_C, extendable = [true, true]) {
    const below = b.split("").reverse().join("");
    let bestAny = null, bestEnd = null;
    const logSalt = Math.log(equivalentSodiumMolar(conditions));
    const kelvin = tempC + 273.15;
    for (let offset = -(below.length - 1); offset < a.length; offset++) {
      const first = Math.max(0, offset), last = Math.min(a.length, below.length + offset);
      const paired = [];
      for (let i = first; i < last; i++) paired.push(below[i - offset] === COMPLEMENT[a[i]]);
      for (let start = first; start < last; start++) {
        if (!paired[start - first]) continue;
        let dh = 0.0, ds = 0.0;
        for (let k = start; k < last - 1; k++) {
          const here = paired[k - first], next = paired[k + 1 - first];
          if (!here && !next) break;
          const table = here && next ? MATCHED : INTERNAL_MISMATCH;
          const value = stackValue(a.slice(k, k + 2), below.slice(k - offset, k - offset + 2), table);
          dh += value[0];
          ds += value[1];
          if (!next) continue;
          const end = k + 1;
          const top = a.slice(start, end + 1), bottom = below.slice(start - offset, end - offset + 1);
          let regionDh = dh, regionDs = ds;
          for (const base of [top[0], top[top.length - 1]]) {
            const init = "GC".includes(base) ? MATCHED.init_w_term_GC : MATCHED.init_w_term_AT;
            regionDh += init[0];
            regionDs += init[1];
          }
          let allPaired = true;
          for (let p = start; p <= end; p++) if (!paired[p - first]) { allPaired = false; break; }
          if (allPaired && top === reverseComplement(top)) {
            regionDh += MATCHED.symmetry_correction[0];
            regionDs += MATCHED.symmetry_correction[1];
          }
          const length = end - start + 1;
          const dg = 1000 * regionDh - kelvin * (regionDs + 0.368 * (length - 1) * logSalt);
          const aEnd = extendable[0] && end === a.length - 1;
          const bEnd = extendable[1] && start - offset === 0;
          if (bestAny !== null && dg >= bestAny.dG && !(aEnd || bEnd)) continue;
          const region = { dG: dg, dH: regionDh, dS: regionDs, top, bottom, a_start: start,
                           b_start: b.length - 1 - (end - offset), length };
          if (bestAny === null || dg < bestAny.dG) bestAny = region;
          if ((aEnd || bEnd) && (bestEnd === null || dg < bestEnd.dG)) bestEnd = region;
        }
      }
    }
    return { any: bestAny, three_prime: bestEnd };
  }

  // dimer_check
  function dimerCheck(oligos, conditions, cache = null) {
    const names = Object.keys(oligos);
    const dimers = [], problems = [];
    for (let x = 0; x < names.length; x++) {
      for (let y = x; y < names.length; y++) {
        const first = names[x], second = names[y];
        const extendable = [first !== "probe", second !== "probe"];
        const key = `${oligos[first]}|${oligos[second]}|${extendable}`;
        let scan;
        if (cache !== null && cache.has(key)) {
          scan = cache.get(key);
        } else {
          scan = dimerScan(oligos[first], oligos[second], conditions, DIMER_TEMP_C, extendable);
          if (cache !== null) cache.set(key, scan);
        }
        const pair = first === second ? first : `${first}/${second}`;
        dimers.push({ pair, any: scan.any, three_prime: scan.three_prime });
        if (scan.any && scan.any.dG < DIMER_MAX_DG) problems.push(`${pair} dimer`);
        if (scan.three_prime && scan.three_prime.dG < DIMER_3END_MAX_DG) problems.push(`${pair} 3' dimer`);
      }
    }
    return [dimers, problems];
  }

  // perfect_tm (duplex_energy for a perfect duplex, then melting_temperature)
  function perfectTm(oligo, conditions) {
    for (const base of oligo) if (!KNOWN_BASES.includes(base)) return null;
    const n = oligo.length;
    let bottom = "";
    for (const base of oligo) bottom += COMPLEMENT[base];
    let dh = 0.0, ds = 0.0;
    for (let i = 0; i < n - 1; i++) {
      const value = stackValue(oligo.slice(i, i + 2), bottom.slice(i, i + 2));
      dh += value[0];
      ds += value[1];
    }
    for (const end of [oligo[0], oligo[n - 1]]) {
      const value = "GC".includes(end) ? MATCHED.init_w_term_GC : MATCHED.init_w_term_AT;
      dh += value[0];
      ds += value[1];
    }
    const selfComplementary = oligo === reverseComplement(oligo);
    if (selfComplementary) {
      dh += MATCHED.symmetry_correction[0];
      ds += MATCHED.symmetry_correction[1];
    }
    // melting_temperature
    ds = ds + 0.368 * (n - 1) * Math.log(equivalentSodiumMolar(conditions));
    let concentration = conditions.oligo_nM * 1e-9;
    if (!selfComplementary) concentration /= 4;
    return 1000 * dh / (ds + R * Math.log(concentration)) - 273.15;
  }

  // product_tm
  function productTm(product, conditions = QPCR_CONDITIONS) {
    let gc = 0;
    for (const base of product) if (base === "G" || base === "C") gc++;
    const gcPercent = 100 * gc / product.length;
    return 81.5 + 16.6 * Math.log10(equivalentSodiumMolar(conditions)) + 0.41 * gcPercent - 600 / product.length;
  }

  function isKnown(base) {
    return base === "A" || base === "C" || base === "G" || base === "T";
  }

  // reverse_complement
  function reverseComplement(seq) {
    let out = "";
    for (let i = seq.length - 1; i >= 0; i--) out += COMPLEMENT[seq[i]];
    return out;
  }

  // base_matches
  function baseMatches(genomeBase, primerBase) {
    return IUPAC[primerBase].includes(genomeBase);
  }

  // genome_vs_primer
  function genomeVsPrimer(genomeBase, primerBase) {
    const couldBe = IUPAC[genomeBase];
    if (couldBe === undefined) return "unknown";
    const accepted = IUPAC[primerBase];
    let all = true, any = false;
    for (const base of couldBe) {
      if (accepted.includes(base)) any = true;
      else all = false;
    }
    if (all) return "match";
    if (!any) return "mismatch";
    return "unknown";
  }

  // expand_degenerate
  function expandDegenerate(seq) {
    let options = [""];
    for (const base of seq) {
      const next = [];
      for (const option of options) for (const choice of IUPAC[base]) next.push(option + choice);
      options = next;
    }
    return options;
  }

  // ===========================================================================
  // FIRST MODEL, SECTION 2: K-MER INDEX
  // ===========================================================================

  // build_kmer_index
  function buildKmerIndex(target, k) {
    const index = new Map();
    for (let position = 0; position + k <= target.length; position++) {
      const kmer = target.slice(position, position + k);
      if (!index.has(kmer)) index.set(kmer, []);
      index.get(kmer).push(position);
    }
    return index;
  }

  // seed_hits
  function seedHits(index, seed) {
    const hits = [];
    for (const plain of expandDegenerate(seed)) {
      const found = index.get(plain);
      if (found) hits.push(...found);
    }
    return hits;
  }

  // ===========================================================================
  // FIRST MODEL, SECTION 3 + 3b: COMPARING A SITE, SEVERITY, 3' RISK
  // ===========================================================================

  // compare_site
  function compareSite(target, start, pattern) {
    const mismatches = [], unknowns = [];
    for (let offset = 0; offset < pattern.length; offset++) {
      const position = start + offset;
      if (position < 0 || position >= target.length) {
        unknowns.push(offset);
        continue;
      }
      const verdict = genomeVsPrimer(target[position], pattern[offset]);
      if (verdict === "mismatch") mismatches.push(offset);
      else if (verdict === "unknown") unknowns.push(offset);
    }
    return [mismatches, unknowns];
  }

  // distances_from_3prime
  function distancesFrom3prime(offsets, length, strand) {
    const out = strand === "+" ? offsets.map((o) => length - o) : offsets.map((o) => o + 1);
    return out.sort((a, b) => a - b);
  }

  const SEVERE_PAIRS = new Set(["AA", "GA", "AG", "CC"]);   // primer + template

  // mismatch_severity
  function mismatchSeverity(primerBase, genomeBase) {
    const template = COMPLEMENT[genomeBase];
    for (const p of IUPAC[primerBase]) {
      for (const t of IUPAC[template]) if (SEVERE_PAIRS.has(p + t)) return "severe";
    }
    return "mild";
  }

  // describe_changes
  function describeChanges(target, pattern, pairs, strand) {
    const changes = pairs.map(([offset, position]) => {
      let primerBase = pattern[offset], genomeBase = target[position];
      if (strand === "-") {
        primerBase = COMPLEMENT[primerBase];
        genomeBase = COMPLEMENT[genomeBase];
      }
      return {
        from_3prime: distancesFrom3prime([offset], pattern.length, strand)[0],
        primer: primerBase,
        genome: genomeBase,
        severity: mismatchSeverity(primerBase, genomeBase),
      };
    });
    return changes.sort((a, b) => a.from_3prime - b.from_3prime);
  }

  const THREE_PRIME_WINDOW = 3;
  const RISK_ORDER = [null, "possible", "low", "high"];

  // three_prime_risk
  function threePrimeRisk(changes, unknownFrom3prime) {
    const near = changes.filter((c) => c.from_3prime <= THREE_PRIME_WINDOW);
    if (near.some((c) => c.severity === "severe")) return "high";
    if (near.length) return "low";
    if (unknownFrom3prime.some((d) => d <= THREE_PRIME_WINDOW)) return "possible";
    return null;
  }

  // worst_risk
  function worstRisk(...risks) {
    return risks.reduce((a, b) => (RISK_ORDER.indexOf(b) > RISK_ORDER.indexOf(a) ? b : a));
  }

  // ===========================================================================
  // FIRST MODEL, SECTION 4: FINDING BINDING SITES ON BOTH STRANDS
  // ===========================================================================

  const MIN_KNOWN_BASES = 10;

  // Positions of N/ambiguous bases, worked out once per target (the Python
  // loops over the whole target every time; the answer is the same).
  const unknownCache = { target: null, positions: [] };
  function unknownPositions(target) {
    if (unknownCache.target !== target) {
      const positions = [];
      for (let i = 0; i < target.length; i++) if (!isKnown(target[i])) positions.push(i);
      unknownCache.target = target;
      unknownCache.positions = positions;
    }
    return unknownCache.positions;
  }

  // candidate_starts
  function candidateStarts(pattern, target, index, k, maxMismatches) {
    const pieces = maxMismatches + 1;
    const pieceLength = Math.floor(pattern.length / pieces);
    const candidates = new Set();
    for (let piece = 0; piece < pieces; piece++) {
      const offset = piece * pieceLength;
      for (const hit of seedHits(index, pattern.slice(offset, offset + k))) candidates.add(hit - offset);
    }
    for (const position of unknownPositions(target)) {
      for (let s = position - pattern.length + 1; s <= position; s++) candidates.add(s);
    }
    for (let s = -(pattern.length - MIN_KNOWN_BASES); s < 0; s++) candidates.add(s);
    for (let s = target.length - pattern.length + 1; s < target.length - MIN_KNOWN_BASES + 1; s++) {
      candidates.add(s);
    }
    return candidates;
  }

  // find_binding_sites
  function findBindingSites(primerName, primer, target, index, k, maxMismatches) {
    if (Math.floor(primer.length / (maxMismatches + 1)) < k) {
      throw new Error(`${primerName}: ${primer.length} bases is too short to seed with k=${k} ` +
                      `while allowing ${maxMismatches} mismatches; lower k or the tolerance.`);
    }
    const sites = [];
    for (const [strand, pattern] of [["+", primer], ["-", reverseComplement(primer)]]) {
      const starts = [...candidateStarts(pattern, target, index, k, maxMismatches)].sort((a, b) => a - b);
      for (const start of starts) {
        const [mismatches, unknowns] = compareSite(target, start, pattern);
        const known = pattern.length - unknowns.length;
        if (mismatches.length > maxMismatches || known < MIN_KNOWN_BASES) continue;
        let outcome, reason;
        if (!unknowns.length) [outcome, reason] = ["binds", ""];
        else if (mismatches.length + unknowns.length <= maxMismatches) {
          [outcome, reason] = ["binds", "binds even if every unknown base is a mismatch"];
        } else if (start < 0) [outcome, reason] = ["not assessable", "runs off the start of the sequence"];
        else if (start + pattern.length > target.length) {
          [outcome, reason] = ["not assessable", "runs off the end of the sequence"];
        } else [outcome, reason] = ["not assessable", "overlaps N/ambiguous bases"];
        const changes = describeChanges(target, pattern, mismatches.map((o) => [o, start + o]), strand);
        const unknownFrom3prime = distancesFrom3prime(unknowns, pattern.length, strand);
        sites.push({
          primer: primerName, strand, start, end: start + pattern.length, outcome, reason,
          mismatches: mismatches.length, changes, indel: null,
          three_prime_risk: threePrimeRisk(changes, unknownFrom3prime),
        });
      }
    }
    return sites;
  }

  // ===========================================================================
  // FIRST MODEL, SECTION 5 + 5b: SMALL INSERTION/DELETION IN THE 5' HALF
  // ===========================================================================

  const MAX_BINDING_INDEL = 1;

  // gapped_pairs
  function gappedPairs(length, leftStart, gapAt, gapLength, kind) {
    const pairs = [];
    if (kind === "none") {
      for (let o = 0; o < length; o++) pairs.push([o, leftStart + o]);
      return [pairs, leftStart + length];
    }
    for (let o = 0; o < gapAt; o++) pairs.push([o, leftStart + o]);
    if (kind === "deletion") {
      for (let o = gapAt + gapLength; o < length; o++) pairs.push([o, leftStart + o - gapLength]);
      return [pairs, leftStart + length - gapLength];
    }
    for (let o = gapAt; o < length; o++) pairs.push([o, leftStart + o + gapLength]);
    return [pairs, leftStart + length + gapLength];
  }

  // placements
  function placements(anchor, length, maxIndel) {
    const out = [[anchor, 0, 0, "none"]];
    for (let gapLength = 1; gapLength <= maxIndel; gapLength++) {
      for (let gapAt = 1; gapAt < length; gapAt++) {
        for (const leftStart of [anchor, anchor - gapLength]) out.push([leftStart, gapAt, gapLength, "insertion"]);
        if (gapAt + gapLength < length) {
          for (const leftStart of [anchor, anchor + gapLength]) out.push([leftStart, gapAt, gapLength, "deletion"]);
        }
      }
    }
    return out;
  }

  // describe_indel
  function describeIndel(kind, gapAt, gapLength, length, strand) {
    const offsets = [];
    if (kind === "deletion") for (let o = gapAt; o < gapAt + gapLength; o++) offsets.push(o);
    else offsets.push(gapAt - 1, gapAt);
    const span = distancesFrom3prime(offsets, length, strand);
    return { kind, length: gapLength, from_3prime: [span[0], span[span.length - 1]] };
  }

  // anchors_for
  function anchorsFor(pattern, index, k) {
    const anchors = new Set();
    for (let offset = 0; offset + k <= pattern.length; offset++) {
      for (const hit of seedHits(index, pattern.slice(offset, offset + k))) anchors.add(hit - offset);
    }
    return anchors;
  }

  // find_indel_sites
  function findIndelSites(primerName, primer, target, index, k, maxMismatches, existing) {
    if (maxMismatches < 1) return [];
    const found = [];
    for (const [strand, pattern] of [["+", primer], ["-", reverseComplement(primer)]]) {
      const length = pattern.length;
      const candidates = [];
      for (const anchor of anchorsFor(pattern, index, k)) {
        for (const [leftStart, gapAt, gapLength, kind] of placements(anchor, length, MAX_BINDING_INDEL)) {
          if (kind === "none") continue;
          const indel = describeIndel(kind, gapAt, gapLength, length, strand);
          if (indel.from_3prime[0] <= Math.floor(length / 2)) continue;
          const [pairs, end] = gappedPairs(length, leftStart, gapAt, gapLength, kind);
          if (leftStart < 0 || end > target.length) continue;
          let allKnown = true;
          for (let p = leftStart; p < end; p++) if (!isKnown(target[p])) { allKnown = false; break; }
          if (!allKnown) continue;
          const mismatched = pairs.filter(([o, p]) => !baseMatches(target[p], pattern[o]));
          const events = mismatched.length + 1;
          if (events <= maxMismatches) candidates.push([events, leftStart, end, indel, mismatched]);
        }
      }
      const taken = existing.filter((s) => s.primer === primerName && s.strand === strand)
                            .map((s) => [s.start, s.end]);
      candidates.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
      for (const [events, leftStart, end, indel, mismatched] of candidates) {
        if (taken.some(([ts, te]) => leftStart < te && ts < end)) continue;
        taken.push([leftStart, end]);
        const changes = describeChanges(target, pattern, mismatched, strand);
        found.push({
          primer: primerName, strand, start: leftStart, end, outcome: "binds", reason: "",
          mismatches: events, changes, indel, three_prime_risk: threePrimeRisk(changes, []),
        });
      }
    }
    return found;
  }

  // ===========================================================================
  // FIRST MODEL, SECTIONS 6-8: PRODUCTS, STATUS, run_model (products only)
  // ===========================================================================

  // pair_into_products
  function pairIntoProducts(sites, minLength, maxLength) {
    const usable = sites.filter((s) => s.primer === "forward" || s.primer === "reverse");
    const plus = usable.filter((s) => s.strand === "+");
    const minus = usable.filter((s) => s.strand === "-");
    const products = [];
    for (const left of plus) {
      for (const right of minus) {
        if (right.end <= left.start) continue;
        const length = right.end - left.start;
        if (minLength <= length && length <= maxLength && right.start >= left.start) {
          products.push({
            left, right, start: left.start, end: right.end, length,
            three_prime_risk: worstRisk(left.three_prime_risk, right.three_prime_risk),
          });
        }
      }
    }
    products.sort((a, b) => a.start - b.start || a.length - b.length);
    return products;
  }

  // oligo_outcome
  function oligoOutcome(name, sites) {
    const outcomes = new Set(sites.filter((s) => s.primer === name).map((s) => s.outcome));
    if (outcomes.has("binds")) return "binds";
    if (outcomes.has("not assessable")) return "not assessable";
    return "mismatched";
  }

  // probe_inside
  function probeInside(product, sites) {
    const inside = sites.filter((s) => s.primer === "probe" && s.start >= product.start && s.end <= product.end);
    return oligoOutcome("probe", inside);
  }

  // product_status
  function productStatus(product, probeGiven) {
    if (probeGiven && product.probe === "mismatched") return "no signal";
    const parts = [product.left.outcome, product.right.outcome];
    if (probeGiven) parts.push(product.probe);
    if (parts.includes("not assessable")) return "not assessable";
    if (product.three_prime_risk === "high") return "at risk";
    return "product";
  }

  // run_model (sites and products; the designer needs nothing else)
  function runModel(forward, reverse, target, probe, maxMismatches, k, minLength, maxLength) {
    const index = buildKmerIndex(target, k);
    const oligos = [["forward", forward], ["reverse", reverse]];
    if (probe) oligos.push(["probe", probe]);
    let sites = [];
    for (const [name, seq] of oligos) {
      const plain = findBindingSites(name, seq, target, index, k, maxMismatches);
      sites = sites.concat(plain, findIndelSites(name, seq, target, index, k, maxMismatches, plain));
    }
    const products = pairIntoProducts(sites, minLength, maxLength);
    for (const product of products) {
      if (probe) product.probe = probeInside(product, sites);
      product.status = productStatus(product, Boolean(probe));
    }
    return { sites, products };
  }

  // ===========================================================================
  // PRIMER DESIGN, SECTION 1: SETTINGS (same values as primer_design.py)
  // ===========================================================================

  const K = 6;
  const OLIGO_GC = [0.40, 0.60];
  const IDEAL_GC = 0.50;
  const TARGET_GC_LOW_BELOW = 0.40;
  const GC_MEDIUM_ABOVE = 0.60;
  const GC_TOO_HIGH_ABOVE = 0.70;
  const PRIMER_LENGTH = [18, 24];
  const PROBE_LENGTH = [18, 30];
  const PRODUCT_RANGE = [60, 150];
  const MIN_PRODUCT = 54;
  const QPCR_PRODUCT = [60, 150];
  const MAX_ALLOWED_MISMATCHES = 2;
  const DESIGNS_WANTED = 3;
  // temperature (version 4): each oligo's length is aimed at annealing + these
  const PRIMER_TM_TARGET_ABOVE = 8.0;
  const PROBE_TM_TARGET_ABOVE = 15.0;
  const ANNEALING_RANGE = [45.0, 72.0];

  // ===========================================================================
  // PRIMER DESIGN, SECTION 2: READING THE INPUT
  // ===========================================================================

  // clean_sequence
  function cleanSequence(text) {
    const lines = text.split(/\r?\n/).filter((line) => !line.startsWith(">"));
    return lines.join("").replace(/\s+/g, "").toUpperCase();
  }

  // Website only: the Python crashes on a letter that is not DNA; here the
  // user gets a message instead.
  function checkLetters(name, seq) {
    const bad = [...new Set(seq)].filter((c) => !(c in IUPAC)).sort();
    if (bad.length) {
      throw new Error(`the ${name} contains letters that are not DNA: ${bad.join(" ")}`);
    }
  }

  // find_target
  function findTarget(genome, target) {
    let start = genome.indexOf(target);
    if (start >= 0) return [start, start + target.length, "+"];
    start = genome.indexOf(reverseComplement(target));
    if (start >= 0) return [start, start + target.length, "-"];
    throw new Error("the target gene is not in the genome (checked both strands)");
  }

  // ===========================================================================
  // PRIMER DESIGN, SECTION 2b: GC CONTENT
  // ===========================================================================

  // gc_fraction
  function gcFraction(seq) {
    let known = 0, gc = 0;
    for (const base of seq) {
      if (!isKnown(base)) continue;
      known++;
      if (base === "G" || base === "C") gc++;
    }
    return known ? gc / known : 0.0;
  }

  // gc_level
  function gcLevel(seq) {
    const gc = gcFraction(seq);
    if (gc > GC_TOO_HIGH_ABOVE) return "too high";
    if (gc > GC_MEDIUM_ABOVE) return "medium";
    if (gc < TARGET_GC_LOW_BELOW) return "low";
    return "ok";
  }

  // gc_percent
  function gcPercent(seq) {
    return `${(100 * gcFraction(seq)).toFixed(1)} %`;
  }

  // check_target_gc
  function checkTargetGc(target) {
    const level = gcLevel(target);
    if (level === "too high") {
      throw new Error(
        `the target's GC content is ${gcPercent(target)}, above ` +
        `${Math.round(100 * GC_TOO_HIGH_ABOVE)} %: primers here would bind too strongly, ` +
        "also at partly matching sites, and the GC-rich template folds on itself " +
        "and copies poorly. Choose a less GC-rich target.");
    }
    if (level === "medium") {
      return `target GC content is ${gcPercent(target)} (medium, above ` +
             `${Math.round(100 * GC_MEDIUM_ABOVE)} %): primers may bind more strongly than ` +
             "intended; check the designs before ordering";
    }
    if (level === "low") {
      return `target GC content is ${gcPercent(target)} (low, below ` +
             `${Math.round(100 * TARGET_GC_LOW_BELOW)} %): fewer stretches reach 40-60 % GC, ` +
             "so there may be fewer designs";
    }
    return null;
  }

  // ===========================================================================
  // PRIMER DESIGN, SECTION 3: EVERY POSSIBLE DESIGN IN THE TARGET
  // ===========================================================================

  // known
  function known(seq) {
    for (const base of seq) if (!isKnown(base)) return false;
    return true;
  }

  // usable
  function usable(seq) {
    const gc = gcFraction(seq);
    return known(seq) && OLIGO_GC[0] <= gc && gc <= OLIGO_GC[1];
  }

  // target_tm
  function targetTm(role, annealingC) {
    return annealingC + (role === "probe" ? PROBE_TM_TARGET_ABOVE : PRIMER_TM_TARGET_ABOVE);
  }

  // choose_length: returns [start, end, Tm or null]
  function chooseLength(genome, anchor, lengths, regionStart, regionEnd, grow, role, annealingC, conditions) {
    let best = null;
    for (let length = lengths[0]; length <= lengths[1]; length++) {
      const [start, end] = grow === "right" ? [anchor, anchor + length] : [anchor - length, anchor];
      if (start < regionStart || end > regionEnd) break;
      const window = genome.slice(start, end);
      if (!usable(window)) continue;
      let tm = null;
      let distance;
      if (annealingC === null) {
        distance = Math.abs(gcFraction(window) - IDEAL_GC);
      } else {
        tm = perfectTm(window, conditionsFor(role, conditions));
        if (role !== "probe" && !(annealingC + PRIMER_TM_ABOVE_ANNEALING[0] <= tm &&
                                  tm <= annealingC + PRIMER_TM_ABOVE_ANNEALING[1])) continue;
        distance = Math.abs(tm - targetTm(role, annealingC));
      }
      if (best === null || distance < best[0]) best = [distance, start, end, tm];
    }
    return best === null ? null : [best[1], best[2], best[3]];
  }

  // all_windows (a Map keeps the anchors in the order they were added, like a Python dict)
  function allWindows(genome, anchors, lengths, regionStart, regionEnd, grow, role, annealingC, conditions) {
    const windows = new Map();
    for (const anchor of anchors) {
      const window = chooseLength(genome, anchor, lengths, regionStart, regionEnd, grow, role, annealingC, conditions);
      if (window) windows.set(anchor, window);
    }
    return windows;
  }

  // centred_probe (with temperature: only probes with Tm >= minTm)
  function centredProbe(probes, gapStart, gapEnd, minTm = null) {
    const middle = (gapStart + gapEnd) / 2;
    let best = null;
    for (let start = gapStart; start <= gapEnd - PROBE_LENGTH[0]; start++) {
      const window = probes.get(start);
      if (window === undefined || window[1] > gapEnd) continue;
      if (minTm !== null && window[2] < minTm) continue;
      const distance = Math.abs((window[0] + window[1]) / 2 - middle);
      if (best === null || distance < best[0]) best = [distance, window];
    }
    return best === null ? null : best[1];
  }

  function range(from, to) {   // from .. to - 1, like Python's range
    const out = [];
    for (let i = from; i < to; i++) out.push(i);
    return out;
  }

  // match_penalty (round(x, 2) as in Python); shorter products rank higher
  function matchPenalty(forward, reverseWindow, probe, productLength, minProduct) {
    let gcPoints = 0;
    for (const seq of [forward, reverseWindow, probe]) gcPoints += Math.abs(100 * gcFraction(seq) - 100 * IDEAL_GC);
    return Math.round((gcPoints + (productLength - minProduct) / 10) * 100) / 100;
  }

  // temperature_penalty
  function temperaturePenalty(forwardTm, reverseTm, probeTm, productLength, minProduct, annealingC) {
    const tmPoints = Math.abs(forwardTm - targetTm("forward", annealingC)) +
                     Math.abs(reverseTm - targetTm("reverse", annealingC)) +
                     Math.abs(probeTm - targetTm("probe", annealingC));
    return pyRound2(tmPoints + (productLength - minProduct) / 10);
  }

  // Python's round(x, 2) rounds the exact binary value half-to-even; for the
  // Tm-based penalty (not a multiple of 0.05 like the GC one) the plain
  // Math.round is the same except at exact halves, which do not occur.
  function pyRound2(x) {
    return Math.round(x * 100) / 100;
  }

  // possible_designs
  function possibleDesigns(genome, regionStart, regionEnd, minProduct, maxProduct, annealingC = null,
                           conditions = QPCR_CONDITIONS) {
    const [shortest, longest] = [minProduct, maxProduct];
    const temperature = annealingC !== null;
    const forwards = allWindows(genome, range(regionStart, regionEnd), PRIMER_LENGTH,
                                regionStart, regionEnd, "right", "forward", annealingC, conditions);
    const probes = allWindows(genome, range(regionStart, regionEnd), PROBE_LENGTH,
                              regionStart, regionEnd, "right", "probe", annealingC, conditions);
    const reverses = allWindows(genome, range(regionStart + 1, regionEnd + 1), PRIMER_LENGTH,
                                regionStart, regionEnd, "left", "reverse", annealingC, conditions);
    const designs = [];
    for (const [forwardStart, [, forwardEnd, forwardTm]] of forwards) {
      for (let productEnd = forwardStart + shortest; productEnd <= forwardStart + longest; productEnd++) {
        const reverseWin = reverses.get(productEnd);
        if (reverseWin === undefined) continue;
        const [reverseStart, , reverseTm] = reverseWin;
        if (reverseStart - forwardEnd < PROBE_LENGTH[0]) continue;
        let minProbeTm = null;
        if (temperature) {
          if (Math.abs(forwardTm - reverseTm) > PRIMER_TM_MAX_DIFFERENCE) continue;
          minProbeTm = Math.max(forwardTm, reverseTm) + PROBE_ABOVE_PRIMERS_C;
        }
        const probeWindow = centredProbe(probes, forwardEnd, reverseStart, minProbeTm);
        if (probeWindow === null) continue;
        const forward = genome.slice(forwardStart, forwardEnd);
        const reverseWindow = genome.slice(reverseStart, productEnd);
        const probe = genome.slice(probeWindow[0], probeWindow[1]);
        let ampliconTm = null;
        if (temperature) {
          ampliconTm = productTm(genome.slice(forwardStart, productEnd), conditions);
          if (ampliconTm > DENATURATION_C - PRODUCT_BELOW_DENATURATION_C) continue;
        }
        const penalty = temperature
          ? temperaturePenalty(forwardTm, reverseTm, probeWindow[2], productEnd - forwardStart, shortest, annealingC)
          : matchPenalty(forward, reverseWindow, probe, productEnd - forwardStart, minProduct);
        designs.push({
          forward: { seq: forward, start: forwardStart, end: forwardEnd },
          reverse: { seq: reverseComplement(reverseWindow), start: reverseStart, end: productEnd },
          probe: { seq: probe, start: probeWindow[0], end: probeWindow[1] },
          product_start: forwardStart,
          product_end: productEnd,
          product_length: productEnd - forwardStart,
          match_penalty: penalty,
          tm: temperature ? { forward: forwardTm, reverse: reverseTm, probe: probeWindow[2], product: ampliconTm } : null,
        });
      }
    }
    designs.sort((a, b) => a.match_penalty - b.match_penalty || a.product_start - b.product_start);
    return designs;
  }

  // ===========================================================================
  // PRIMER DESIGN, SECTION 4: SPECIFICITY CHECK WITH THE FIRST MODEL
  // ===========================================================================

  // binds_once
  function bindsOnce(oligo, genome, index, maxMismatches, cache) {
    if (!cache.has(oligo)) {
      cache.set(oligo, findBindingSites("oligo", oligo, genome, index, K, maxMismatches).length === 1);
    }
    return cache.get(oligo);
  }

  // specific
  function specific(design, genome, index, maxMismatches, cache) {
    const oligos = [design.forward.seq, design.reverse.seq, design.probe.seq];
    if (!oligos.every((o) => bindsOnce(o, genome, index, maxMismatches, cache))) return false;
    const result = runModel(oligos[0], oligos[1], genome, oligos[2], maxMismatches, K,
                            20, Math.max(1000, design.product_length));
    const products = result.products;
    return products.length === 1 && products[0].status === "product" &&
           products[0].start === design.product_start && products[0].end === design.product_end;
  }

  // ===========================================================================
  // PRIMER DESIGN, SECTION 5: THE WHOLE DESIGN IN ONE CALL
  // ===========================================================================

  // design_primers
  function designPrimers(genomeText, targetText, maxMismatches = 2, minProduct = PRODUCT_RANGE[0],
                         maxProduct = PRODUCT_RANGE[1], designsWanted = DESIGNS_WANTED,
                         annealingC = ANNEALING_TEMP_C, conditions = QPCR_CONDITIONS) {
    const genome = cleanSequence(genomeText);
    const target = cleanSequence(targetText);
    if (!genome) throw new Error("the genome box is empty");
    if (!target) throw new Error("the target gene box is empty");
    checkLetters("genome", genome);
    checkLetters("target gene", target);
    if (!(Number.isInteger(maxMismatches) && 0 <= maxMismatches && maxMismatches <= MAX_ALLOWED_MISMATCHES)) {
      throw new Error(`mismatch allowance must be 0 to ${MAX_ALLOWED_MISMATCHES}`);
    }
    if (!Number.isInteger(minProduct) || !Number.isInteger(maxProduct)) {
      throw new Error("the product lengths must be whole numbers");
    }
    if (minProduct < MIN_PRODUCT) {
      throw new Error(`the shortest product must be at least ${MIN_PRODUCT} bp ` +
                      "(room for an 18-base forward primer, probe and reverse primer)");
    }
    if (maxProduct < minProduct) {
      throw new Error(`the longest product (${maxProduct} bp) is shorter than the ` +
                      `shortest (${minProduct} bp)`);
    }
    if (annealingC !== null && !(typeof annealingC === "number" && isFinite(annealingC) &&
        ANNEALING_RANGE[0] <= annealingC && annealingC <= ANNEALING_RANGE[1])) {
      throw new Error(`the annealing temperature must be ${ANNEALING_RANGE[0]}-${ANNEALING_RANGE[1]} C ` +
                      "(qPCR usually anneals at 55-65 C)");
    }
    if (annealingC !== null) checkConditions(conditions);
    if (target.length < minProduct) {
      throw new Error(`the target (${target.length} bases) is shorter than the ` +
                      `shortest product asked for (${minProduct} bp)`);
    }

    // -- the GC check comes first: a too GC-rich target stops here
    const warnings = [];
    const targetWarning = checkTargetGc(target);
    if (targetWarning) warnings.push(targetWarning);
    if (minProduct < QPCR_PRODUCT[0] || maxProduct > QPCR_PRODUCT[1]) {
      warnings.push(`the product range ${minProduct}-${maxProduct} bp goes outside ` +
                    `${QPCR_PRODUCT[0]}-${QPCR_PRODUCT[1]} bp, where qPCR works best`);
    }

    const [regionStart, regionEnd, strand] = findTarget(genome, target);
    const designs = possibleDesigns(genome, regionStart, regionEnd, minProduct, maxProduct, annealingC, conditions);

    const index = buildKmerIndex(genome, K);
    const cache = new Map();
    const chosen = [];
    let checked = 0, dimerRejected = 0;
    const dimerCache = new Map();
    for (const design of designs) {
      if (chosen.length === designsWanted) break;
      if (chosen.some((d) => design.product_start < d.product_end && d.product_start < design.product_end)) continue;
      if (annealingC !== null) {
        const [dimers, problems] = dimerCheck(
          { forward: design.forward.seq, reverse: design.reverse.seq, probe: design.probe.seq },
          conditionsFor("forward", conditions), dimerCache);
        if (problems.length) { dimerRejected++; continue; }
        design.dimers = dimers;
      }
      checked++;
      if (specific(design, genome, index, maxMismatches, cache)) chosen.push(design);
    }

    if (!designs.length) {
      const [shortest, longest] = [minProduct, maxProduct];
      const temperatureRules = annealingC === null ? "" :
        ` and whose temperatures suit annealing at ${fmtG(annealingC)} C ` +
        `(primers ${fmtG(annealingC + PRIMER_TM_ABOVE_ANNEALING[0])}-` +
        `${fmtG(annealingC + PRIMER_TM_ABOVE_ANNEALING[1])} C and within ` +
        `${fmtG(PRIMER_TM_MAX_DIFFERENCE)} C of each other, probe ` +
        `${fmtG(PROBE_ABOVE_PRIMERS_C)} C above them)`;
      warnings.push("no design: the target has no forward, reverse and probe stretches with " +
                    `${Math.round(100 * OLIGO_GC[0])}-${Math.round(100 * OLIGO_GC[1])} % GC` +
                    `${temperatureRules} that fit together into a ${shortest}-${longest} bp product`);
    }

    return {
      warnings,
      target: { start: regionStart, end: regionEnd, strand },
      genome_length: genome.length,
      product_range: [minProduct, maxProduct],
      annealing_c: annealingC,
      possible: designs.length,
      checked,
      dimer_rejected: dimerRejected,
      conditions: annealingC === null ? null : Object.assign({}, conditions),
      designs: chosen,
    };
  }

  // Python's f"{x:g}" for the numbers used here (whole or one decimal)
  function fmtG(x) {
    return String(Math.round(x * 1e6) / 1e6);
  }

  root.PrimerDesign = {
    designPrimers, gcFraction, gcPercent, reverseComplement, cleanSequence,
    perfectTm, productTm, conditionsFor, targetTm, dimerScan, QPCR_CONDITIONS, CONDITION_LIMITS,
    DIMER_MAX_DG, DIMER_3END_MAX_DG,
    MAX_ALLOWED_MISMATCHES, PRODUCT_RANGE, MIN_PRODUCT, QPCR_PRODUCT, PRIMER_LENGTH, PROBE_LENGTH,
    ANNEALING_TEMP_C, ANNEALING_RANGE, PRIMER_TM_ABOVE_ANNEALING, PRIMER_TM_MAX_DIFFERENCE,
    PROBE_ABOVE_PRIMERS_C, DENATURATION_C, PRODUCT_BELOW_DENATURATION_C, PRIMER_NM, PROBE_NM,
  };
})(typeof window !== "undefined" ? window : this);
