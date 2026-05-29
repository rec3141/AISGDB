#!/usr/bin/env python3
"""GLANSIS species list → NCBI Datasets API → tiered genome inventory.

Reads the USGS NAS Darwin Core archive for the GLANSIS subset
(https://nas.er.usgs.gov/ipt/archive.do?r=nas_glansis), extracts the
distinct `scientificName` values, queries NCBI Datasets v2 + Entrez for
each, and classifies the species into one of five tiers:

    Tier A  Chromosome-scale RefSeq or GenBank reference assembly
    Tier B  Scaffold-, contig-, or other draft-level assembly
    Tier C  No assembly, but transcriptome / EST records exist
            (TSA or EST GenBank divisions — often a precursor to WGS)
    Tier D  No assembly or transcriptome — markers only (COI/18S/mtDNA)
    Tier E  No usable public sequence

Writes:
    glansis_inventory.tsv   one row per species with tier + accession + length
    glansis_inventory.html  rendered table for github.io
    glansis_inventory.json  machine-readable counts + per-tier species lists

A local cache at <cache_path> stops re-queries on subsequent runs; pass
--refresh to invalidate it.
"""

import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

NCBI_GENOME_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxon}/dataset_report"
)
NCBI_ESEARCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    "?db=nuccore&term={term}&retmode=json&retmax=0"
)
# Entrez filters for assembled transcriptomes (tsa_master) and EST records
# (gbdiv_est). Both flag "expressed-region-only" sequence resources that pre-
# date a full WGS assembly for the same taxon — useful as Tier-C references
# for shotgun mapping but biased toward coding regions.
TRANSCRIPTOME_FILTER = "(gbdiv_est[PROP] OR tsa_master[PROP])"

USER_AGENT = "danaSeq-GLANSIS-inventory/1.0 (https://github.com/rec3141/danaSeq)"


def fetch_json(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8", "replace"))


def species_from_dwc(occurrence_path):
    """Return the distinct scientificName values from the GLANSIS DwC archive."""
    names = set()
    with open(occurrence_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = header.index("scientificName")
        for line in f:
            row = line.rstrip("\n").split("\t")
            if len(row) <= idx:
                continue
            n = row[idx].strip()
            if n:
                names.add(n)
    return sorted(names)


def classify_assemblies(reports):
    """Pick the most informative assembly across all reports for one species."""
    if not reports:
        return None
    # Rank by assembly level. RefSeq beats GenBank within a level.
    level_rank = {
        "Chromosome": 4, "Complete Genome": 5,
        "Scaffold": 2, "Contig": 1, None: 0,
    }
    def key(r):
        info = r.get("assembly_info", {}) or {}
        lvl = level_rank.get(info.get("assembly_level"), 0)
        is_refseq = (r.get("accession") or "").startswith("GCF_")
        return (lvl, int(is_refseq))
    best = max(reports, key=key)
    info = best.get("assembly_info", {}) or {}
    stats = best.get("assembly_stats", {}) or {}
    org = best.get("organism", {}) or {}
    return {
        "accession": best.get("accession"),
        "assembly_name": info.get("assembly_name"),
        "assembly_level": info.get("assembly_level"),
        "total_length": stats.get("total_sequence_length"),
        "organism_name": org.get("organism_name"),
        "refseq_category": info.get("refseq_category"),
        "n_reports": len(reports),
    }


def query_ncbi_genome(taxon):
    url = NCBI_GENOME_URL.format(taxon=urllib.parse.quote(taxon))
    try:
        d = fetch_json(url)
    except Exception as e:
        return None, str(e)
    return d.get("reports") or [], None


def query_entrez_count(term):
    url = NCBI_ESEARCH_URL.format(term=urllib.parse.quote(term))
    try:
        d = fetch_json(url)
        return int(d["esearchresult"].get("count", 0))
    except Exception:
        return 0


def query_nuccore_count(taxon):
    return query_entrez_count(f'{taxon}[Organism]')


def query_transcriptome_count(taxon):
    return query_entrez_count(f'{taxon}[Organism] AND {TRANSCRIPTOME_FILTER}')


def tier_for(genome_info, nuccore_count, transcriptome_count):
    """Map (best assembly, transcriptome count, total nuccore count) → tier."""
    if genome_info:
        lvl = (genome_info.get("assembly_level") or "").lower()
        if "chromosome" in lvl or "complete" in lvl:
            return "A", "Chromosome-scale WGS"
        return "B", f"{genome_info.get('assembly_level') or 'Draft'} WGS"
    if transcriptome_count and transcriptome_count > 0:
        return "C", f"Transcriptome / EST ({transcriptome_count} records)"
    if nuccore_count and nuccore_count > 0:
        if nuccore_count >= 25:
            return "D", f"Marker pool ({nuccore_count} records)"
        return "D", f"Markers only ({nuccore_count} records)"
    return "E", "No public sequence"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--occurrence", default="/tmp/glansis/occurrence.txt",
                    help="Path to extracted GLANSIS DwC occurrence.txt")
    ap.add_argument("--out-tsv",  default="/matika/projects/project_zebra/docs/glansis_inventory.tsv")
    ap.add_argument("--out-html", default="/matika/projects/project_zebra/docs/glansis_inventory.html")
    ap.add_argument("--out-json", default="/matika/projects/project_zebra/docs/glansis_inventory.json")
    ap.add_argument("--cache",    default="/tmp/glansis_ncbi_cache.json")
    ap.add_argument("--refresh",  action="store_true",
                    help="Invalidate the cache and re-query NCBI.")
    ap.add_argument("--delay-ms", type=int, default=120,
                    help="Sleep between NCBI calls (default 120 ms).")
    ap.add_argument("--limit",    type=int, default=0,
                    help="Stop after N species (debug).")
    args = ap.parse_args()

    species = species_from_dwc(args.occurrence)
    if args.limit:
        species = species[: args.limit]
    print(f"[inventory] GLANSIS species: {len(species)}", file=sys.stderr)

    cache = {}
    if os.path.exists(args.cache) and not args.refresh:
        try:
            with open(args.cache) as f:
                cache = json.load(f)
            print(f"[inventory] cache hits available for {len(cache)} species",
                  file=sys.stderr)
        except Exception:
            pass

    results = []
    for i, name in enumerate(species, 1):
        entry = cache.get(name) or {}
        need_genome = "genome" not in entry
        # Fill in any missing fields. Transcriptome + nuccore are only meaningful
        # when there's no genome, so skip them in that case to save API calls.
        if need_genome:
            reports, err = query_ncbi_genome(name)
            time.sleep(args.delay_ms / 1000.0)
            entry["genome"] = classify_assemblies(reports) if reports else None
            entry["err"] = err
        if not entry.get("genome"):
            if "transcriptome_count" not in entry:
                entry["transcriptome_count"] = query_transcriptome_count(name)
                time.sleep(args.delay_ms / 1000.0)
            if "nuccore_count" not in entry:
                entry["nuccore_count"] = query_nuccore_count(name)
                time.sleep(args.delay_ms / 1000.0)
        else:
            entry.setdefault("transcriptome_count", 0)
            entry.setdefault("nuccore_count", 0)
        cache[name] = entry
        if i % 25 == 0:
            with open(args.cache, "w") as f:
                json.dump(cache, f)
            print(f"[inventory] {i}/{len(species)}", file=sys.stderr)
        tier, label = tier_for(
            entry.get("genome"),
            entry.get("nuccore_count") or 0,
            entry.get("transcriptome_count") or 0,
        )
        results.append({
            "scientific_name": name,
            "tier": tier,
            "tier_label": label,
            **(entry.get("genome") or {}),
            "transcriptome_count": entry.get("transcriptome_count") or 0,
            "nuccore_count": entry.get("nuccore_count") or 0,
        })

    with open(args.cache, "w") as f:
        json.dump(cache, f)

    # Sort: WGS first, then transcriptome, then markers, then nothing.
    # Within tier, sort alphabetically.
    tier_order = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}
    results.sort(key=lambda r: (tier_order[r["tier"]], r["scientific_name"]))

    # Counts
    counts = defaultdict(int)
    for r in results:
        counts[r["tier"]] += 1

    # ---- TSV ----
    cols = ["tier", "scientific_name", "accession", "assembly_level",
            "assembly_name", "total_length", "refseq_category",
            "transcriptome_count", "nuccore_count"]
    with open(args.out_tsv, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in results:
            f.write("\t".join(str(r.get(c, "") or "") for c in cols) + "\n")

    # ---- JSON (machine-readable summary) ----
    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_species": len(results),
        "tier_counts": {k: counts[k] for k in "ABCDE"},
        "species": results,
    }
    with open(args.out_json, "w") as f:
        json.dump(summary, f, indent=2)

    # ---- HTML ----
    write_html(args.out_html, results, counts)

    print(
        f"[inventory] done. Tier A={counts['A']}  B={counts['B']}  "
        f"C={counts['C']}  D={counts['D']}  E={counts['E']}",
        file=sys.stderr,
    )


def fmt_bp(n):
    if not n or n == "":
        return ""
    try:
        n = int(n)
    except (TypeError, ValueError):
        return str(n)
    if n >= 1e9: return f"{n/1e9:.2f} Gb"
    if n >= 1e6: return f"{n/1e6:.0f} Mb"
    if n >= 1e3: return f"{n/1e3:.0f} kb"
    return str(n)


def write_html(path, results, counts):
    tier_meta = {
        "A": ("Chromosome-scale WGS",                "#047857"),
        "B": ("Draft (scaffold / contig) WGS",       "#1d4ed8"),
        "C": ("Transcriptome / EST",                 "#7c3aed"),
        "D": ("Marker / mitochondrion only",         "#b45309"),
        "E": ("No usable public sequence",           "#be123c"),
    }
    n_total = len(results)
    rows_html = []
    for r in results:
        tier = r["tier"]
        label, color = tier_meta[tier]
        acc = r.get("accession") or ""
        # Link to the BioProject search rather than the canonical
        # /datasets/genome/<acc>/ page: as of 2026-05, the NCBI Datasets
        # genome landing pages return errors (likely an NIH funding /
        # service-outage issue). BioProject search is a stable fallback
        # that resolves the underlying study. Version suffix (".1") is
        # stripped — the search matches the base accession.
        acc_base = acc.split(".", 1)[0]
        acc_html = (
            f'<a href="https://www.ncbi.nlm.nih.gov/bioproject/?term={acc_base}">{acc}</a>'
            if acc else ""
        )
        tsc = r.get("transcriptome_count") or ""
        nuc = r.get("nuccore_count") or ""
        rows_html.append(f"""
        <tr>
          <td><span class="badge" style="background:{color}1a;color:{color}">{tier}</span></td>
          <td><em>{r["scientific_name"]}</em></td>
          <td>{r.get("assembly_level") or ""}</td>
          <td class="mono">{acc_html}</td>
          <td class="mono right">{fmt_bp(r.get("total_length"))}</td>
          <td class="mono right">{tsc}</td>
          <td class="mono right">{nuc}</td>
        </tr>
        """)

    tier_summary_html = " ".join(
        f'<span class="tier-pill" style="background:{tier_meta[t][1]}1a;color:{tier_meta[t][1]}">'
        f'<b>{t}</b> {tier_meta[t][0]} &middot; <b>{counts[t]}</b> '
        f'({counts[t]/n_total*100:.0f}%)</span>'
        for t in "ABCDE"
    )

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GLANSIS species → NCBI assembly inventory</title>
<style>
  body {{
    font: 15px/1.5 -apple-system, BlinkMacSystemFont, "Inter", "Segoe UI", Roboto, sans-serif;
    color: #1f2937; max-width: 56rem; margin: 2rem auto; padding: 0 1rem;
  }}
  h1 {{ font-size: 1.5rem; margin-bottom: 0.25rem; }}
  .meta {{ color: #475569; font-size: 0.85rem; margin-bottom: 1rem; }}
  .tier-summary {{ margin: 1rem 0 1.5rem; display: flex; flex-wrap: wrap; gap: 0.5rem; }}
  .tier-pill {{
    display: inline-block; padding: 0.4rem 0.8rem; border-radius: 6px;
    font-size: 0.9rem;
  }}
  .filter-bar {{ margin: 1rem 0; }}
  .filter-bar input {{
    padding: 0.4rem 0.7rem; font-size: 0.9rem; border: 1px solid #cbd5e1;
    border-radius: 4px; width: 100%; max-width: 24rem;
  }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; }}
  th, td {{
    text-align: left; padding: 0.4rem 0.6rem; border-bottom: 1px solid #e2e8f0;
    vertical-align: top;
  }}
  th {{ font-weight: 600; color: #1e40af; background: #f8fafc; position: sticky; top: 0; }}
  .mono {{ font-family: "JetBrains Mono", ui-monospace, monospace; font-size: 0.85em; }}
  .right {{ text-align: right; }}
  .badge {{
    display: inline-block; padding: 0.05em 0.55em; border-radius: 9999px;
    font-weight: 600; font-size: 0.85em;
  }}
  a {{ color: #1e40af; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  footer {{ margin-top: 2rem; color: #475569; font-size: 0.8rem; }}
</style>
</head><body>
<h1>GLANSIS species &rarr; NCBI assembly inventory</h1>
<div class="meta">
  Generated {time.strftime("%Y-%m-%d")} &middot;
  {n_total} species from the
  <a href="https://nas.er.usgs.gov/ipt/resource?r=nas_glansis">USGS NAS GLANSIS Darwin Core archive</a>
  &middot; per-species lookup via the
  <a href="https://api.ncbi.nlm.nih.gov/datasets/">NCBI Datasets v2 API</a>
  and <a href="https://eutils.ncbi.nlm.nih.gov/">Entrez E-utilities</a>.
</div>

<div class="tier-summary">{tier_summary_html}</div>

<div class="filter-bar">
  <input id="filter" placeholder="Filter species (case-insensitive substring)" oninput="filterRows(event)">
</div>

<table>
  <thead>
    <tr>
      <th>Tier</th><th>Species</th><th>Assembly level</th>
      <th>Accession</th><th class="right">Length</th>
      <th class="right">EST/TSA</th><th class="right">nuccore</th>
    </tr>
  </thead>
  <tbody id="rows">
    {''.join(rows_html)}
  </tbody>
</table>

<footer>
  Built by <a href="https://github.com/rec3141/danaSeq/blob/main/docs/glansis_inventory.py">glansis_inventory.py</a>
  &middot; methodology described in the
  <a href="glansis-ncbi-panel-proposal.html">GLANSIS &rarr; NCBI panel proposal</a>.
</footer>

<script>
  function filterRows(e) {{
    const q = e.target.value.toLowerCase();
    const rows = document.querySelectorAll('#rows tr');
    for (const r of rows) {{
      r.style.display = r.textContent.toLowerCase().includes(q) ? '' : 'none';
    }}
  }}
</script>

</body></html>
"""
    with open(path, "w") as f:
        f.write(html)


if __name__ == "__main__":
    main()
