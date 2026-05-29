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
| [`data/`](data/) | Vendored snapshots of the four data files that back the GLANSIS Risk-Assessment and Regulations Explorers. See table below. |

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

Per-species rows are enriched from two GLANSIS sources, all joined on
case-folded `(Genus, Species)`:

| Column | Source | Coverage |
|---|---|---|
| `common_name`  | Modal `Common Names` across all assessments for the species (`RALevel2_v6.txt`) | 197/370 (53%) |
| `group`        | Modal `Group` (Fishes / Plants / Mollusk-Bivalve / …), with sing/plural variants collapsed (`RALevel2_v6.txt`) | 252/370 (68%) |
| `risk`         | Most-severe verdict across all `Overall` strings, by priority *Invasive > High > Watchlist > Moderate > Low* (`RALevel2_v6.txt`) | 193/370 (52%) |
| `reg_level`    | Most-restrictive regulation level across all jurisdictions, by priority *Prohibited > Restricted > Other* (`invasiveRegs.txt` + Manitoba MR 173/2015 Sched. A) | 134/370 (36%) |
| `reg_n_juris`  | Total number of regulatory listings — e.g. round goby (*Neogobius melanostomus*) is Prohibited in 9 jurisdictions and Restricted in 3 (`invasiveRegs.txt` + Manitoba) | same |

The HTML view exposes Group, Risk, and Regulation level as dropdown
filters, mirroring the GLANSIS Tier-2 Explorer UI; the regulation cell
hover-tooltip lists the per-jurisdiction breakdown.

### Vendored data files

The GLANSIS Explorer pages load these four TSVs at runtime via PapaParse.
None has a DOI, NCEI listing, or data.noaa.gov entry — and the file names
bump versions destructively (older `_v*` 404 the moment a new one ships) —
so we mirror them here. All four are tab-delimited; encoding noted below.

| File | Source URL | Encoding | What it is |
|---|---|---|---|
| [`data/RALevel2_v6.txt.gz`](data/RALevel2_v6.txt.gz) | <https://www.glerl.noaa.gov/glansis/data/RALevel2_v6.txt> | UTF-16 LE | Tier-2 risk-assessment table — 4,997 rows × 20 cols. Powers `group`, `common_name`, `risk`. |
| [`data/invasiveRegs.txt.gz`](data/invasiveRegs.txt.gz) | <https://www.glerl.noaa.gov/glansis/data/invasiveRegs.txt> | UTF-8 | Regulatory listings — 757 rows × 19 cols across 13 jurisdictions (WI, MN, OH, NY, IN, PA, QC, IL, US, MI, ON, Canada, …). Powers `reg_level`, `reg_n_juris`. |
| [`data/RA_Content_v3_forLevel2.txt.gz`](data/RA_Content_v3_forLevel2.txt.gz) | <https://www.glerl.noaa.gov/glansis/data/RA_Content_v3_forLevel2.txt> | UTF-8 | Methodology lookup — 18 rows describing each risk-assessment method (USFWS ERSS, Canadian AqWRA, …). Not joined into the inventory; archived for future use as per-assessment tooltip text. |
| [`data/IllinoisWhiteList.txt.gz`](data/IllinoisWhiteList.txt.gz) | <https://www.glerl.noaa.gov/glansis/data/IllinoisWhiteList.txt> | CP-1252 | Illinois aquaculture-approved list — 367 species. Not joined: it's the *inverse* signal (permitted in IL), almost no overlap with the GLANSIS DwC list. |
| [`data/manitoba_schedule_a.tsv`](data/manitoba_schedule_a.tsv) | extracted from PDF (see below) | UTF-8 | Manitoba *Aquatic Invasive Species Regulation* (MR 173/2015) Schedule A — 128 prohibited species in 4 groups (Fish 81, Invertebrates 24, Plants 21, Algae 2). 125 carry a clean binomial; 62 overlap with the GLANSIS DwC list and are folded into `reg_level`/`reg_n_juris` as one additional Prohibited jurisdiction. The 63 MB-only species (Russian/Beluga sturgeons, killer shrimp, yabby, etc.) are documented in the TSV but not yet expanded as inventory rows. |
| [`data/manitoba_schedule_a.parse.py`](data/manitoba_schedule_a.parse.py) | (this repo) | — | Reproducible parser — converts the `pdftotext -layout` dump back into a clean TSV. Re-run when the regulation is amended. |
| [`data/manitoba_AIS_regulation_173-2015.pdf.gz`](data/manitoba_AIS_regulation_173-2015.pdf.gz) | <https://web2.gov.mb.ca/laws/regs/current/_pdf-regs.php?reg=173/2015> | PDF | Pinned snapshot of the source regulation (current 25 May 2023). Manitoba doesn't ship a structured table — CanLII blocks scrapers, and the official site only offers HTML and PDF. The PDF parses more reliably than the HTML. |

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
