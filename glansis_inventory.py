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
import csv
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

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

# GLANSIS Tier-2 risk-assessment CSV (UTF-16 LE, tab-delimited). Powers the
# Group, Common Name, and Risk columns. The upstream URL bumps versions
# destructively (older _v* files 404 immediately on a bump) and has no DOI /
# NCEI archive, so we vendor a gzipped copy in this repo and fall back to it
# if the live URL is gone.
GLANSIS_RA_URL = "https://www.glerl.noaa.gov/glansis/data/RALevel2_v6.txt"

# GLANSIS's Group field has uncontrolled singular/plural and " - " vs "-"
# variants. Collapse to the most-common spelling.
GROUP_ALIASES = {
    "Plant":                 "Plants",
    "Nematode":              "Nematodes",
    "Rotifer":               "Rotifers",
    "Mollusks-Bivalves":     "Mollusk-Bivalve",
    "Mollusks-Gastropods":   "Mollusk-Gastropoda",
    "Mollusk - Slugs":       "Mollusk-Gastropoda",
    "Crustaceans-Amphipods": "Crustacean-Amphipod",
    "Crustaceans-Copepod":   "Crustacean-Copepod",
    "Crustaceans-Cladocerans": "Crustacean-Cladoceran",
    "Crustaceans-Crayfish":  "Crustacean-Crayfish",
    "Crustaceans-Mysids":    "Crustacean-Mysid",
    "Crustaceans-Tanaids":   "Crustacean-Tanaid",
    "Crustaceans-Crab":      "Crustacean-Crab",
    "Reptile - Turtle":      "Reptile",
    "Annelids-Polychaetes":  "Annelids",
    "Cnidarian":             "Coelenterates",
    "Platyhelminthes":       "Flatworm",
    "Myxozoa-Malacosporea":  "Myxosporean",
}

# Severity order for the Risk column — Invasive first because it's a
# confirmed-fact label, not a predicted score.
RISK_PRIORITY = ["Invasive", "High", "Watchlist", "Moderate", "Low"]

# GLANSIS regulations table — one row per (species, jurisdiction). UTF-8 TSV.
# Same URL-stability caveat as RALevel2_v6.txt: vendor a snapshot.
GLANSIS_REGS_URL = "https://www.glerl.noaa.gov/glansis/data/invasiveRegs.txt"

# Most-restrictive regulation level wins for the displayed cell.
REG_PRIORITY = ["Prohibited", "Restricted", "Other"]

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


def _normalize_group(g):
    g = (g or "").strip()
    return GROUP_ALIASES.get(g, g)


_RE_QUANT = re.compile(r"\b(High|Moderate|Medium|Low)\b", re.I)


def _extract_verdict(overall):
    """Map a single 'Overall' string to a severity label, or None."""
    t = (overall or "").strip()
    if not t:
        return None
    lt = t.lower()
    # GLANSIS categorical tags (whole-cell values).
    if lt == "invasive":
        return "Invasive"
    if lt.startswith("watchlist") or lt == "recommend watchlist":
        return "Watchlist"
    # Quantitative ladder used by USFWS ERSS and similar tools.
    m = _RE_QUANT.search(t)
    if m:
        v = m.group(1).title()
        return "Moderate" if v == "Medium" else v
    # Substring-level invasive detection — catches longer phrases like
    # "Invasive with Benefits" without colliding with "Non-Invasive".
    if "invasive" in lt and "non" not in lt:
        return "Invasive"
    return None


def load_ra_file(path):
    """Read GLANSIS RALevel2 CSV (UTF-16 LE TSV, or .gz of same).

    Returns a dict keyed by (genus_lower, species_lower) →
        {group, common_name, risk, n_ra, organizations}
    """
    if not path or not os.path.exists(path):
        return {}
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as gz:
            raw = gz.read()
    else:
        with open(path, "rb") as f:
            raw = f.read()
    # The upstream file is UTF-16 LE with BOM.
    text = raw.decode("utf-16")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    by_species = defaultdict(list)
    for row in reader:
        g = (row.get("Genus") or "").strip()
        s = (row.get("Species") or "").strip()
        if not (g and s):
            continue
        by_species[(g.lower(), s.lower())].append(row)

    enrichment = {}
    for key, rows in by_species.items():
        groups = Counter(_normalize_group(r.get("Group")) for r in rows
                         if (r.get("Group") or "").strip())
        commons = Counter((r.get("Common Names") or "").strip() for r in rows
                          if (r.get("Common Names") or "").strip())
        orgs = Counter((r.get("Organization") or "").strip() for r in rows
                       if (r.get("Organization") or "").strip())
        verdicts = [v for v in (_extract_verdict(r.get("Overall")) for r in rows) if v]
        risk = next((v for v in RISK_PRIORITY if v in verdicts), None)
        enrichment[key] = {
            "group":        groups.most_common(1)[0][0] if groups else "",
            "common_name":  commons.most_common(1)[0][0] if commons else "",
            "risk":         risk or "",
            "n_ra":         len(rows),
            "organizations": sorted(orgs),
        }
    return enrichment


def _enrichment_key(scientific_name):
    """Best-effort (genus, species) from a DwC scientificName."""
    parts = scientific_name.split()
    if len(parts) < 2:
        return None
    return (parts[0].lower(), parts[1].lower())


def load_regs(path):
    """Read GLANSIS invasiveRegs.txt (UTF-8 TSV or .gz of same).

    Returns a dict keyed by (genus_lower, species_lower) →
        {top_level, n_jurisdictions, jurisdictions_by_level}
    The dataset has family-level rows (empty Genus/Species) too — those are
    skipped here, which costs some recall on multi-species families.
    """
    if not path or not os.path.exists(path):
        return {}
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            text = fh.read()
    else:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    by_species = defaultdict(list)
    for row in reader:
        g = (row.get("Genus") or "").strip()
        s = (row.get("Species") or "").strip()
        if not (g and s):
            continue
        by_species[(g.lower(), s.lower())].append(row)

    regs = {}
    for key, rows in by_species.items():
        by_level = defaultdict(list)
        for r in rows:
            level = (r.get("RegulationLevel") or "").strip() or "Other"
            jur = (r.get("Jurisdiction") or "").strip()
            if jur:
                by_level[level].append(jur)
        top = next((lvl for lvl in REG_PRIORITY if lvl in by_level),
                   next(iter(by_level), ""))
        regs[key] = {
            "top_level":         top,
            "n_jurisdictions":   sum(len(v) for v in by_level.values()),
            "by_level":          {lvl: sorted(jur) for lvl, jur in by_level.items()},
        }
    return regs


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
    ap.add_argument("--ra",
                    default=os.path.join(os.path.dirname(__file__) or ".",
                                         "data", "RALevel2_v6.txt.gz"),
                    help="GLANSIS Tier-2 risk-assessment CSV (UTF-16 TSV or .gz). "
                         f"Live source: {GLANSIS_RA_URL}")
    ap.add_argument("--regs",
                    default=os.path.join(os.path.dirname(__file__) or ".",
                                         "data", "invasiveRegs.txt.gz"),
                    help="GLANSIS regulatory-listings TSV (UTF-8 or .gz). "
                         f"Live source: {GLANSIS_REGS_URL}")
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

    ra_enrich = load_ra_file(args.ra)
    print(f"[inventory] RA enrichment loaded for {len(ra_enrich)} (genus, species) pairs",
          file=sys.stderr)
    regs_enrich = load_regs(args.regs)
    print(f"[inventory] regs loaded for {len(regs_enrich)} (genus, species) pairs",
          file=sys.stderr)

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
        ek = _enrichment_key(name) or ("", "")
        ra = ra_enrich.get(ek, {})
        rg = regs_enrich.get(ek, {})
        results.append({
            "scientific_name": name,
            "tier": tier,
            "tier_label": label,
            **(entry.get("genome") or {}),
            "transcriptome_count": entry.get("transcriptome_count") or 0,
            "nuccore_count": entry.get("nuccore_count") or 0,
            "group":       ra.get("group", ""),
            "common_name": ra.get("common_name", ""),
            "risk":        ra.get("risk", ""),
            "n_ra":        ra.get("n_ra", 0),
            "reg_level":   rg.get("top_level", ""),
            "reg_n_juris": rg.get("n_jurisdictions", 0),
            "reg_by_level": rg.get("by_level", {}),
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
    cols = ["tier", "scientific_name", "common_name", "group", "risk", "n_ra",
            "reg_level", "reg_n_juris",
            "accession", "assembly_level",
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
    risk_color = {
        "Invasive":  "#be123c",
        "High":      "#dc2626",
        "Watchlist": "#ea580c",
        "Moderate":  "#ca8a04",
        "Low":       "#16a34a",
    }
    reg_color = {
        "Prohibited": "#be123c",
        "Restricted": "#ea580c",
        "Other":      "#64748b",
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
        grp = r.get("group") or ""
        cname = r.get("common_name") or ""
        risk = r.get("risk") or ""
        n_ra = r.get("n_ra") or 0
        rc = risk_color.get(risk, "#94a3b8")
        risk_cell = (
            f'<span class="badge" style="background:{rc}1a;color:{rc}">{risk}</span>'
            f' <span class="muted">({n_ra})</span>' if risk else
            (f'<span class="muted">({n_ra})</span>' if n_ra else "")
        )
        species_cell = f'<em>{r["scientific_name"]}</em>'
        if cname:
            species_cell += f'<br><span class="muted small">{cname}</span>'
        reg_level = r.get("reg_level") or ""
        reg_n = r.get("reg_n_juris") or 0
        reg_by_level = r.get("reg_by_level") or {}
        reg_tooltip = " | ".join(
            f"{lvl}: {', '.join(juris)}"
            for lvl in REG_PRIORITY if lvl in reg_by_level
            for juris in [reg_by_level[lvl]]
        ) or ""
        rcc = reg_color.get(reg_level, "#94a3b8")
        regs_cell = (
            f'<span class="badge" style="background:{rcc}1a;color:{rcc}" '
            f'title="{reg_tooltip}">{reg_level}</span>'
            f' <span class="muted">×{reg_n}</span>' if reg_level else ""
        )

        # Per-column sort keys. Numeric where it makes sense; pre-ranked
        # for categorical columns with a natural severity order. Empty
        # strings sort to the end regardless of direction (handled in JS).
        tier_rank = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}[tier]
        risk_rank = RISK_PRIORITY.index(risk) if risk in RISK_PRIORITY else 99
        reg_rank  = REG_PRIORITY.index(reg_level) if reg_level in REG_PRIORITY else 99
        # Composite: severity major, count minor (more jurisdictions = more
        # severe within a tier). Encoded so ascending sort = "worst first".
        reg_sort  = reg_rank * 10000 - (reg_n if reg_level else 0)
        asm_rank  = {"Complete Genome": 0, "Chromosome": 1,
                     "Scaffold": 2, "Contig": 3}.get(
                         r.get("assembly_level") or "", 99)
        len_sort  = r.get("total_length") or ""
        tsc_sort  = tsc if isinstance(tsc, int) else (tsc or "")
        nuc_sort  = nuc if isinstance(nuc, int) else (nuc or "")
        species_sort = r["scientific_name"].lower()
        group_sort   = (grp or "").lower()
        acc_sort     = (acc or "").lower()

        rows_html.append(f"""
        <tr data-group="{grp}" data-risk="{risk}" data-reg="{reg_level}">
          <td data-sort="{tier_rank}"><span class="badge" style="background:{color}1a;color:{color}">{tier}</span></td>
          <td data-sort="{species_sort}">{species_cell}</td>
          <td class="small" data-sort="{group_sort}">{grp}</td>
          <td data-sort="{risk_rank}">{risk_cell}</td>
          <td data-sort="{reg_sort}">{regs_cell}</td>
          <td data-sort="{asm_rank}">{r.get("assembly_level") or ""}</td>
          <td class="mono" data-sort="{acc_sort}">{acc_html}</td>
          <td class="mono right" data-sort="{len_sort}">{fmt_bp(r.get("total_length"))}</td>
          <td class="mono right" data-sort="{tsc_sort}">{tsc}</td>
          <td class="mono right" data-sort="{nuc_sort}">{nuc}</td>
        </tr>
        """)

    tier_summary_html = " ".join(
        f'<span class="tier-pill" style="background:{tier_meta[t][1]}1a;color:{tier_meta[t][1]}">'
        f'<b>{t}</b> {tier_meta[t][0]} &middot; <b>{counts[t]}</b> '
        f'({counts[t]/n_total*100:.0f}%)</span>'
        for t in "ABCDE"
    )

    # Group + risk dropdown option lists, ordered by frequency.
    group_counts = Counter(r.get("group","") for r in results if r.get("group"))
    risk_counts  = Counter(r.get("risk","")  for r in results if r.get("risk"))
    group_options = "".join(
        f'<option value="{g}">{g} ({n})</option>'
        for g, n in group_counts.most_common()
    )
    risk_options = "".join(
        f'<option value="{r}">{r} ({n})</option>'
        for r, n in sorted(risk_counts.items(),
                           key=lambda x: RISK_PRIORITY.index(x[0])
                                         if x[0] in RISK_PRIORITY else 99)
    )
    n_with_ra = sum(1 for r in results if r.get("n_ra"))
    n_with_reg = sum(1 for r in results if r.get("reg_n_juris"))
    reg_options = "".join(
        f'<option value="{lvl}">{lvl}</option>' for lvl in REG_PRIORITY
        if any(r.get("reg_level") == lvl for r in results)
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
  .filter-bar {{
    margin: 1rem 0; display: flex; flex-wrap: wrap; gap: 0.5rem; align-items: center;
  }}
  .filter-bar input, .filter-bar select {{
    padding: 0.4rem 0.7rem; font-size: 0.9rem; border: 1px solid #cbd5e1;
    border-radius: 4px;
  }}
  .filter-bar input {{ flex: 1; min-width: 18rem; max-width: 28rem; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; }}
  th, td {{
    text-align: left; padding: 0.4rem 0.6rem; border-bottom: 1px solid #e2e8f0;
    vertical-align: top;
  }}
  th {{
    font-weight: 600; color: #1e40af; background: #f8fafc; position: sticky; top: 0;
    cursor: pointer; user-select: none;
  }}
  th:hover {{ background: #eef2ff; }}
  th.sorted-asc::after  {{ content: " \\25B2"; opacity: 0.6; font-size: 0.8em; }}
  th.sorted-desc::after {{ content: " \\25BC"; opacity: 0.6; font-size: 0.8em; }}
  .mono {{ font-family: "JetBrains Mono", ui-monospace, monospace; font-size: 0.85em; }}
  .right {{ text-align: right; }}
  .small {{ font-size: 0.85em; }}
  .muted {{ color: #64748b; }}
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
  <input id="filter" placeholder="Filter species (substring)" oninput="applyFilters()">
  <select id="groupFilter" onchange="applyFilters()">
    <option value="">All groups</option>
    {group_options}
  </select>
  <select id="riskFilter" onchange="applyFilters()">
    <option value="">All risk</option>
    {risk_options}
  </select>
  <select id="regFilter" onchange="applyFilters()">
    <option value="">All regulation</option>
    {reg_options}
  </select>
</div>
<div class="meta">
  Group, Common Name, and Risk are joined from the
  <a href="https://www.glerl.noaa.gov/glansis/raT2Explorer.html">GLANSIS Tier&nbsp;2 Risk-Assessment Clearinghouse</a>
  (<code>RALevel2_v6.txt</code>, vendored at <code>data/</code>) &mdash;
  {n_with_ra} of {n_total} species have ≥1 assessment.
  Risk reduces all assessments per species to a single label
  (priority: Invasive &gt; High &gt; Watchlist &gt; Moderate &gt; Low); count
  in parentheses is the number of source assessments.
  Regulated is joined from
  <a href="https://www.glerl.noaa.gov/glansis/raT2Explorer.html"><code>invasiveRegs.txt</code></a>
  &mdash; {n_with_reg} of {n_total} species have at least one regulatory listing
  in a Great Lakes jurisdiction. Hover the badge for the per-jurisdiction
  breakdown; cell shows the most-restrictive level and total jurisdictions.
</div>

<table>
  <thead>
    <tr>
      <th onclick="sortBy(0)">Tier</th>
      <th onclick="sortBy(1)">Species</th>
      <th onclick="sortBy(2)">Group</th>
      <th onclick="sortBy(3)">Risk</th>
      <th onclick="sortBy(4)">Regulated</th>
      <th onclick="sortBy(5)">Assembly level</th>
      <th onclick="sortBy(6)">Accession</th>
      <th class="right" onclick="sortBy(7)">Length</th>
      <th class="right" onclick="sortBy(8)">EST/TSA</th>
      <th class="right" onclick="sortBy(9)">nuccore</th>
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
  function applyFilters() {{
    const q   = document.getElementById('filter').value.toLowerCase();
    const grp = document.getElementById('groupFilter').value;
    const rsk = document.getElementById('riskFilter').value;
    const reg = document.getElementById('regFilter').value;
    const rows = document.querySelectorAll('#rows tr');
    for (const r of rows) {{
      let show = true;
      if (q   && !r.textContent.toLowerCase().includes(q)) show = false;
      if (grp && r.dataset.group !== grp) show = false;
      if (rsk && r.dataset.risk  !== rsk) show = false;
      if (reg && r.dataset.reg   !== reg) show = false;
      r.style.display = show ? '' : 'none';
    }}
  }}

  // Click-to-sort. First click ascending, second descending.
  // Empty cells always sink to the bottom regardless of direction.
  let sortState = {{ col: null, dir: 1 }};
  function sortBy(col) {{
    if (sortState.col === col) sortState.dir = -sortState.dir;
    else                       {{ sortState.col = col; sortState.dir = 1; }}
    const tbody = document.getElementById('rows');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    rows.sort((a, b) => {{
      const av = a.children[col].dataset.sort ?? '';
      const bv = b.children[col].dataset.sort ?? '';
      if (av === '' && bv !== '') return  1;
      if (bv === '' && av !== '') return -1;
      const an = parseFloat(av), bn = parseFloat(bv);
      if (!isNaN(an) && !isNaN(bn)) return (an - bn) * sortState.dir;
      return av.localeCompare(bv) * sortState.dir;
    }});
    rows.forEach(r => tbody.appendChild(r));
    document.querySelectorAll('th').forEach((th, i) => {{
      th.classList.toggle('sorted-asc',  i === col && sortState.dir > 0);
      th.classList.toggle('sorted-desc', i === col && sortState.dir < 0);
    }});
  }}
</script>

</body></html>
"""
    with open(path, "w") as f:
        f.write(html)


if __name__ == "__main__":
    main()
