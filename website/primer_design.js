// =============================================================================
// PRIMER DESIGN FOR THE WEBSITE: JavaScript version of journal/primer_design.py
// =============================================================================
// 02-120 project: Do early SARS-CoV-2 qPCR primers still match later variants?
//
// A browser cannot run Python, so this file translates, function by function,
//   - journal/primer_design.py (version 2), all of it, and
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

  // choose_length
  function chooseLength(genome, anchor, lengths, regionStart, regionEnd, grow) {
    let best = null;
    for (let length = lengths[0]; length <= lengths[1]; length++) {
      const [start, end] = grow === "right" ? [anchor, anchor + length] : [anchor - length, anchor];
      if (start < regionStart || end > regionEnd) break;
      const window = genome.slice(start, end);
      if (!usable(window)) continue;
      const distance = Math.abs(gcFraction(window) - IDEAL_GC);
      if (best === null || distance < best[0]) best = [distance, start, end];
    }
    return best === null ? null : [best[1], best[2]];
  }

  // all_windows (a Map keeps the anchors in the order they were added, like a Python dict)
  function allWindows(genome, anchors, lengths, regionStart, regionEnd, grow) {
    const windows = new Map();
    for (const anchor of anchors) {
      const window = chooseLength(genome, anchor, lengths, regionStart, regionEnd, grow);
      if (window) windows.set(anchor, window);
    }
    return windows;
  }

  // centred_probe
  function centredProbe(probes, gapStart, gapEnd) {
    const middle = (gapStart + gapEnd) / 2;
    let best = null;
    for (let start = gapStart; start <= gapEnd - PROBE_LENGTH[0]; start++) {
      const window = probes.get(start);
      if (window === undefined || window[1] > gapEnd) continue;
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

  // possible_designs
  function possibleDesigns(genome, regionStart, regionEnd, minProduct, maxProduct) {
    const [shortest, longest] = [minProduct, maxProduct];
    const forwards = allWindows(genome, range(regionStart, regionEnd), PRIMER_LENGTH,
                                regionStart, regionEnd, "right");
    const probes = allWindows(genome, range(regionStart, regionEnd), PROBE_LENGTH,
                              regionStart, regionEnd, "right");
    const reverses = allWindows(genome, range(regionStart + 1, regionEnd + 1), PRIMER_LENGTH,
                                regionStart, regionEnd, "left");
    const designs = [];
    for (const [forwardStart, [, forwardEnd]] of forwards) {
      for (let productEnd = forwardStart + shortest; productEnd <= forwardStart + longest; productEnd++) {
        const reverseWin = reverses.get(productEnd);
        if (reverseWin === undefined) continue;
        const reverseStart = reverseWin[0];
        if (reverseStart - forwardEnd < PROBE_LENGTH[0]) continue;
        const probeWindow = centredProbe(probes, forwardEnd, reverseStart);
        if (probeWindow === null) continue;
        const forward = genome.slice(forwardStart, forwardEnd);
        const reverseWindow = genome.slice(reverseStart, productEnd);
        const probe = genome.slice(probeWindow[0], probeWindow[1]);
        designs.push({
          forward: { seq: forward, start: forwardStart, end: forwardEnd },
          reverse: { seq: reverseComplement(reverseWindow), start: reverseStart, end: productEnd },
          probe: { seq: probe, start: probeWindow[0], end: probeWindow[1] },
          product_start: forwardStart,
          product_end: productEnd,
          product_length: productEnd - forwardStart,
          match_penalty: matchPenalty(forward, reverseWindow, probe, productEnd - forwardStart, minProduct),
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
                         maxProduct = PRODUCT_RANGE[1], designsWanted = DESIGNS_WANTED) {
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
    const designs = possibleDesigns(genome, regionStart, regionEnd, minProduct, maxProduct);

    const index = buildKmerIndex(genome, K);
    const cache = new Map();
    const chosen = [];
    let checked = 0;
    for (const design of designs) {
      if (chosen.length === designsWanted) break;
      if (chosen.some((d) => design.product_start < d.product_end && d.product_start < design.product_end)) continue;
      checked++;
      if (specific(design, genome, index, maxMismatches, cache)) chosen.push(design);
    }

    if (!designs.length) {
      const [shortest, longest] = [minProduct, maxProduct];
      warnings.push("no design: the target has no forward, reverse and probe stretches with " +
                    `${Math.round(100 * OLIGO_GC[0])}-${Math.round(100 * OLIGO_GC[1])} % GC ` +
                    `that fit together into a ${shortest}-${longest} bp product`);
    }

    return {
      warnings,
      target: { start: regionStart, end: regionEnd, strand },
      genome_length: genome.length,
      product_range: [minProduct, maxProduct],
      possible: designs.length,
      checked,
      designs: chosen,
    };
  }

  root.PrimerDesign = {
    designPrimers, gcFraction, gcPercent, reverseComplement, cleanSequence,
    MAX_ALLOWED_MISMATCHES, PRODUCT_RANGE, MIN_PRODUCT, QPCR_PRODUCT, PRIMER_LENGTH, PROBE_LENGTH,
  };
})(typeof window !== "undefined" ? window : this);
