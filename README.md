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
| [`glansis_inventory.html`](glansis_inventory.html) | Filterable + sortable per-species inventory — all 370 GLANSIS species tiered |
| [`glansis_inventory.{tsv,json}`](glansis_inventory.tsv) | Machine-readable variants of the GLANSIS inventory |
| [`glansis_inventory.py`](glansis_inventory.py) | Reproducible inventory builder (also exports helpers used by the Canadian inventory) |
| [`canadian_ais_inventory.html`](canadian_ais_inventory.html) | Sibling inventory page for Canada's federal SOR/2015-121 Schedule (Parts 2 + 3) — 101 species |
| [`canadian_ais_inventory.{tsv,json}`](canadian_ais_inventory.tsv) | Machine-readable variants of the Canadian inventory |
| [`canadian_ais_inventory.py`](canadian_ais_inventory.py) | Thin entry point that imports `glansis_inventory` helpers and runs the federal species list |
| [`data/`](data/) | Vendored snapshots of every external data file the inventories depend on. See table below. |

## First-pass results

NCBI sequence status as of 2026-05-29, across two species lists:

**GLANSIS** (NOAA Great Lakes watchlist, 370 species):

| Tier | Definition | Species | Share |
|---|---|---|---|
| A | Chromosome-scale RefSeq or GenBank reference assembly | 95 | 26% |
| B | Scaffold- or contig-level draft assembly | 30 | 8% |
| C | Transcriptome / EST records (TSA or EST GenBank divisions) | 13 | 4% |
| D | Markers only (COI / 18S / mtDNA records) | 92 | 25% |
| E | No usable public sequence | 140 | 38% |

**Canada SOR/2015-121** (federal AISR Schedule Parts 2 + 3, 101 species):

| Tier | Species | Share |
|---|---|---|
| A | 47 | 47% |
| B | 14 | 14% |
| C |  3 |  3% |
| D | 32 | 32% |
| E |  5 |  5% |

The federal list is more sequence-rich — half its species already have
chromosome-scale assemblies because it's biased toward economically- and
ecologically-important fish that have been the targets of dedicated WGS
projects (the Asian carps, the prairie sturgeons, the Canadian sport
fish in Part 3). The GLANSIS list is dragged down by ~25% obscure plant
markers-only and ~38% no-sequence species — a feature of a broader
basin-wide watchlist.

Tier A alone is a 95-species GLANSIS panel + 47 federal-list species
(with substantial overlap on zebra/quagga/Asian carps/round goby etc.).

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
| `reg_level`    | Most-restrictive regulation level across all jurisdictions, by priority *Prohibited > Restricted > Other* (`invasiveRegs.txt` + federal SOR/2015-121 Schedule Parts 2 + 3) | 122/370 (33%) |
| `reg_n_juris`  | Total number of regulatory listings — federal Part 2 zones are decomposed per-province (BC/AB/SK/MB, ON+MB, etc.) and de-duplicated against GLANSIS-listed jurisdictions. E.g. zebra mussel: 10 Prohibited (Canada-wide + BC/AB/SK/MB province-specific + GLANSIS US states) + 2 Restricted | same |

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
| [`data/canada_sor_2015_121_schedule.tsv`](data/canada_sor_2015_121_schedule.tsv) | parsed from PDF (see below) | UTF-8 | Federal *Aquatic Invasive Species Regulations* (SOR/2015-121) Schedule — 89 Part-2 prohibited species + 14 Part-3 species-at-risk = 103 entries. Each row carries the species' Condition (e.g. "Dead and eviscerated") and the four per-prohibition zones (Importation / Possession / Transport / Release — e.g. zebra mussel possession prohibited in BC/AB/SK/MB; bighead carp prohibited Canada-wide). Drives the Canadian inventory page and is merged into the GLANSIS Regulated column as the authoritative Canadian source. |
| [`data/canada_sor_2015_121.parse.py`](data/canada_sor_2015_121.parse.py) | (this repo) | — | Reproducible parser — uses `pdftohtml -xml` and item-number anchor binding to handle the schedule's wrapping cells and split-glyph quirks. Re-run when the regulation is amended. |
| [`data/canada_sor_2015_121.pdf.gz`](data/canada_sor_2015_121.pdf.gz) | <https://laws-lois.justice.gc.ca/PDF/SOR-2015-121.pdf> | PDF | Pinned snapshot of the source regulation (current April 2026, last amended June 2021). |
| [`data/manitoba_schedule_a.tsv`](data/manitoba_schedule_a.tsv) | parsed from PDF (see below) | UTF-8 | Manitoba MR 173/2015 Schedule A — kept as a historical reference and cross-check. **Not used as inventory input**: the federal SOR/2015-121 covers the same species (and more) with richer per-province scope. Reading the two side-by-side made the relationship explicit: MR 173/2015 is the provincial-enforcement mirror of the federal regulation. |
| [`data/manitoba_schedule_a.parse.py`](data/manitoba_schedule_a.parse.py) | (this repo) | — | Reproducible parser for the MB schedule — kept for the same reason as the data file (historical reference / cross-check against federal). |
| [`data/manitoba_AIS_regulation_173-2015.pdf.gz`](data/manitoba_AIS_regulation_173-2015.pdf.gz) | <https://web2.gov.mb.ca/laws/regs/current/_pdf-regs.php?reg=173/2015> | PDF | Pinned snapshot of the Manitoba regulation (current 25 May 2023). |

## Re-running the inventories

```bash
# pull the GLANSIS Darwin Core archive (drives the GLANSIS page)
curl -sS -L -o /tmp/glansis.zip \
    "https://nas.er.usgs.gov/ipt/archive.do?r=nas_glansis"
mkdir -p /tmp/glansis && cd /tmp/glansis && unzip -o /tmp/glansis.zip

# query NCBI for the GLANSIS list (uses a local cache; pass --refresh to invalidate)
python3 glansis_inventory.py --occurrence /tmp/glansis/occurrence.txt

# query NCBI for the federal SOR/2015-121 list (shares the same cache —
# overlapping species aren't re-queried)
python3 canadian_ais_inventory.py
```

Outputs `glansis_inventory.{tsv,html,json}` and `canadian_ais_inventory.{tsv,html,json}`.
Safe to re-run on a cron — new chromosome-scale assemblies move species
from Tier E/D/C/B → A as they're deposited, and the inventory delta is the news.

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
