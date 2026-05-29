# AISGDB

**Aquatic Invasive Species Genomic Database** — a working space for a curated
whole-genome reference panel for shotgun-metagenomic AIS surveillance.

Browse the rendered docs at **<https://rec3141.github.io/AISGDB/>**.

## Contents

| File | What it is |
|---|---|
| [`index.html`](index.html) | Landing page for github.io |
| [`ais-shotgun-metagenomics-landscape.html`](ais-shotgun-metagenomics-landscape.html) | Literature landscape doc — where shotgun metagenomic AIS surveillance stands in 2026 |
| [`glansis-ncbi-panel-proposal.html`](glansis-ncbi-panel-proposal.html) | Funder-facing pitch — HTML version |
| [`glansis-ncbi-panel-proposal.tex`](glansis-ncbi-panel-proposal.tex) | Same content as a LaTeX article (compiles with `pdflatex`) |
| [`glansis_inventory.html`](glansis_inventory.html) | Filterable per-species inventory — all 370 GLANSIS species tiered |
| [`glansis_inventory.tsv`](glansis_inventory.tsv) | Machine-readable inventory |
| [`glansis_inventory.json`](glansis_inventory.json) | Same as TSV plus per-tier summary |
| [`glansis_inventory.py`](glansis_inventory.py) | Reproducible script — re-runnable as new NCBI assemblies are deposited |
| [`data/RALevel2_v6.txt.gz`](data/RALevel2_v6.txt.gz) | Vendored snapshot of the GLANSIS Tier-2 risk-assessment CSV (UTF-16 TSV, gzipped). Source: `https://www.glerl.noaa.gov/glansis/data/RALevel2_v6.txt`. Mirrored because the upstream URL bumps versions destructively (older `_v*` files 404 immediately) and the dataset has no DOI / NCEI archive. |

## First-pass results

Of the 370 species in the NOAA GLANSIS database, NCBI sequence status as of
2026-05-29:

| Tier | Definition | Species | Share |
|---|---|---|---|
| A | Chromosome-scale RefSeq or GenBank reference assembly | 95 | 26% |
| B | Scaffold- or contig-level draft assembly | 30 | 8% |
| C | Transcriptome / EST records (TSA or EST GenBank divisions) | 13 | 4% |
| D | Markers only (COI / 18S / mtDNA records) | 92 | 25% |
| E | No usable public sequence | 140 | 38% |

Tier A alone is a 95-species panel for the Great Lakes basin — roughly 20×
larger than what's currently deployed in companion pipelines.

Tier C species — those with transcriptome or EST data but no assembly —
are interesting precisely because they have the resources needed for a WGS
project in progress. They're the highest-yield candidates for an
opportunistic "next WGS" prioritisation list.

Per-species rows are also enriched from the **GLANSIS Tier-2 Risk-Assessment
Clearinghouse** (`RALevel2_v6.txt`) with three additional columns:

| Column | Source | Coverage |
|---|---|---|
| `common_name` | Modal `Common Names` across all assessments for the species | 197/370 (53%) |
| `group`       | Modal `Group` (functional group: Fishes, Plants, Mollusk-Bivalve, …) with sing/plural variants collapsed | 252/370 (68%) |
| `risk`        | Most-severe verdict across all `Overall` strings, by priority *Invasive > High > Watchlist > Moderate > Low* | 193/370 (52%) |

The HTML view exposes Group and Risk as dropdown filters, mirroring the
GLANSIS Tier-2 Explorer UI.

## Re-running the inventory

```bash
# pull the GLANSIS Darwin Core archive
curl -sS -L -o /tmp/glansis.zip \
    "https://nas.er.usgs.gov/ipt/archive.do?r=nas_glansis"
mkdir -p /tmp/glansis && cd /tmp/glansis && unzip -o /tmp/glansis.zip

# query NCBI (uses a local cache; pass --refresh to invalidate)
python3 glansis_inventory.py \
    --occurrence /tmp/glansis/occurrence.txt
```

Outputs `glansis_inventory.{tsv,html,json}`. Safe to re-run on a cron — new
chromosome-scale assemblies move species from Tier E/D/C/B → A as they're
deposited, and the inventory delta is the news.

## Related projects

- **[danaSeq](https://github.com/rec3141/danaSeq)** — the nanopore pipeline
  whose `nanopore_live/--mapping_refs` module consumes an AISGDB-shaped
  reference directory. See [`docs/mapping-references.md`](https://github.com/rec3141/danaSeq/blob/main/docs/mapping-references.md)
  for the reference schema and pipeline integration.
- **[microscape.app](https://microscape.app)** — the SPA frontend; the AIS
  view at `/<run>/#ais` surfaces per-sample identity histograms,
  genome-position distributions, and a live HQ-identity cutoff slider.

## License

CC BY 4.0 for documents; MIT for code. See individual files for headers.
