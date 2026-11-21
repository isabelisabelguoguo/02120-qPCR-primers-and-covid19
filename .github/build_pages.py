# =============================================================================
# BUILD THE PUBLIC WEBSITE (GitHub Pages)
# =============================================================================
# The two websites are written for claude.ai, which wraps each page in a
# document skeleton (doctype, charset, viewport, a small reset) when it is
# published. GitHub Pages serves files as they are, so this script makes a
# copy of each site with that skeleton added, in _site/:
#
#   _site/index.html                 home page: links to both versions
#   _site/website/                   Version 1 (GC content only)
#   _site/website_temperature/       Version 2 (with melting temperature)
#
# The pages link to each other with their claude.ai addresses (those only
# work inside the organization); here the links point to the public copies.
#
# The source files in website/ and website_temperature/ are not changed.
# .github/workflows/pages.yml runs this on every push to main and publishes
# _site/. To try it locally:  python3 .github/build_pages.py
# =============================================================================

import os
import shutil

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(PROJECT, "_site")
SITES = {
    "website": "https://claude.ai/artifact/XKd9SSGhpFHH1WXYemP3Bh",
    "website_temperature": "https://claude.ai/artifact/JRi6KEZtGsBqdkP5jUspZ7",
}

# The same skeleton claude.ai adds: charset, viewport (with the side insets
# the pages' gutters expect), and its small reset.
HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<style>:root{color-scheme:light;box-sizing:border-box;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}body{margin:0;padding:0;font:14px -apple-system,BlinkMacSystemFont,sans-serif}img{max-width:100%}[hidden]{display:none!important}</style>
</head>
<body>
"""
TAIL = "\n</body>\n</html>\n"

HOME = """<title>qPCR Primer Designer</title>
<meta name="description" content="Design qPCR primers and probes for a target gene, checked for specificity across the whole genome. 02-120 project on SARS-CoV-2.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Sans+Condensed:wght@600&display=swap">
<style>
:root {
  --bg: #f6f7f5; --surface: #ffffff; --ink: #1b211e; --muted: #5c6862; --line: #d9dfdb;
  --accent: #1f6f5c; --heat: #c2410c;
  --font-body: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --font-head: "IBM Plex Sans Condensed", "IBM Plex Sans", system-ui, sans-serif;
  --font-mono: "IBM Plex Mono", ui-monospace, "SF Mono", Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root { --bg: #121614; --surface: #1a201d; --ink: #e4e9e6; --muted: #9aa7a0; --line: #2d3532;
          --accent: #4fb396; --heat: #f08a4b; color-scheme: dark; }
}
* { box-sizing: border-box; }
body { background: var(--bg); color: var(--ink); font-family: var(--font-body); font-size: 15px; line-height: 1.5; }
.wrap { max-width: 860px; margin: 0 auto; display: grid; gap: 24px;
        padding-left: max(16px, env(safe-area-inset-left, 0px));
        padding-right: max(16px, env(safe-area-inset-right, 0px)); padding-block: 40px 64px; }
.eyebrow { font-family: var(--font-mono); font-size: 12px; letter-spacing: 0.06em; text-transform: uppercase; color: var(--muted); }
h1 { font-family: var(--font-head); font-weight: 600; font-size: 34px; line-height: 1.15; margin: 4px 0 0; text-wrap: balance; }
.lede { margin: 0; color: var(--muted); max-width: 65ch; }
.versions { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px; }
@media (max-width: 600px) { .versions { grid-template-columns: 1fr; } }
a.version { display: grid; gap: 8px; align-content: start; background: var(--surface); border: 1px solid var(--line);
            border-radius: 8px; padding: 20px; color: var(--ink); text-decoration: none; min-width: 0; }
a.version:hover, a.version:focus-visible { border-color: var(--accent); }
a.version:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
.tag { font-family: var(--font-mono); font-size: 12px; color: var(--muted); }
.version h2 { font-family: var(--font-head); font-size: 22px; font-weight: 600; margin: 0; }
.version p { margin: 0; color: var(--muted); font-size: 14px; }
.go { color: var(--accent); font-weight: 600; font-size: 14px; }
.v2 .tag { color: var(--heat); }
footer { font-size: 13px; color: var(--muted); }
footer a { color: var(--accent); }
</style>

<div class="wrap">
  <header>
    <div class="eyebrow">02-120 Programming for Scientists · qPCR primers and SARS-CoV-2</div>
    <h1>qPCR Primer Designer</h1>
    <p class="lede">Paste a genome and a target gene. The designer suggests a forward primer, a reverse
      primer and a probe, and checks with our specificity model that each design makes only one
      product in the whole genome. Both versions run entirely in your browser.</p>
  </header>

  <div class="versions">
    <a class="version v1" href="website/">
      <span class="tag">Version 1 · GC content only</span>
      <h2>qPCR Primer Designer</h2>
      <p>Chooses each primer's length for 40–60 % GC and ranks the shortest products first.
        Does not consider melting temperature.</p>
      <span class="go">Open Version 1 →</span>
    </a>
    <a class="version v2" href="website_temperature/">
      <span class="tag">Version 2 · with melting temperature</span>
      <h2>Temperature-Aware Primer Designer</h2>
      <p>Chooses lengths so primers hold at the qPCR annealing step, keeps the probe above them,
        checks primer dimers, and lets you set the salt and oligo amounts.</p>
      <span class="go">Open Version 2 →</span>
    </a>
  </div>

  <footer>Version 1 led to Version 2: compared with NCBI Primer-BLAST, Version 1 picked primers that melt
    too close to the annealing temperature. Code, tests and the build journal:
    <a href="https://github.com/isabelisabelguoguo/02120-qPCR-primers-and-covid19">GitHub repository</a>.</footer>
</div>
"""


def build():
    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    os.makedirs(OUT)
    for folder in SITES:
        source, target = os.path.join(PROJECT, folder), os.path.join(OUT, folder)
        shutil.copytree(source, target)
        page_path = os.path.join(target, "index.html")
        with open(page_path, encoding="utf-8") as handle:
            page = handle.read()
        for other, claude_url in SITES.items():
            page = page.replace(claude_url, f"../{other}/")
        page = page.replace('target="_blank" rel="noopener" ', "")   # same site now: same tab
        with open(page_path, "w", encoding="utf-8") as handle:
            handle.write(HEAD + page + TAIL)
    with open(os.path.join(OUT, "index.html"), "w", encoding="utf-8") as handle:
        handle.write(HEAD + HOME + TAIL)
    open(os.path.join(OUT, ".nojekyll"), "w").close()   # serve files as they are
    print(f"built {OUT}: index.html, " + ", ".join(f"{f}/" for f in SITES))


if __name__ == "__main__":
    build()
