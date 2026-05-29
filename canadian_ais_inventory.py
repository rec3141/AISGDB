#!/usr/bin/env python3
"""Canadian SOR/2015-121 → NCBI Datasets → tiered genome inventory.

Sibling of glansis_inventory.py — same NCBI classification logic and
enrichment sources (Tier-2 RA, invasiveRegs.txt + the federal regulation
itself), but the species list comes from the federal Aquatic Invasive
Species Regulations Schedule (Parts 2 + 3) rather than the GLANSIS DwC
archive. 89 Part-2 prohibited species + 14 Part-3 species-at-risk = 103
species, tiered into:

    Tier A  Chromosome-scale RefSeq or GenBank reference assembly
    Tier B  Scaffold-, contig-, or other draft-level assembly
    Tier C  No assembly, but transcriptome / EST records exist
    Tier D  No assembly or transcriptome — markers only
    Tier E  No usable public sequence

Writes:
    canadian_ais_inventory.tsv   one row per species with tier + accession
    canadian_ais_inventory.html  rendered table for github.io
    canadian_ais_inventory.json  machine-readable counts + species list

Uses the same NCBI cache file as glansis_inventory by default — many
species overlap and we don't want to re-query them.
"""

import argparse
import os
import sys
import time

import glansis_inventory as gi


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(__file__) or "."
    ap.add_argument("--canada-schedule",
                    default=os.path.join(here, "data",
                                         "canada_sor_2015_121_schedule.tsv"),
                    help="Federal SOR/2015-121 schedule (Parts 2+3, parsed TSV). "
                         f"Source: {gi.CANADA_AIS_REGULATION_URL}")
    ap.add_argument("--ra",
                    default=os.path.join(here, "data", "RALevel2_v6.txt.gz"),
                    help="GLANSIS Tier-2 risk-assessment CSV (UTF-16 TSV or .gz).")
    ap.add_argument("--regs",
                    default=os.path.join(here, "data", "invasiveRegs.txt.gz"),
                    help="GLANSIS regulatory-listings TSV (UTF-8 or .gz).")
    ap.add_argument("--out-tsv",  default=os.path.join(here, "canadian_ais_inventory.tsv"))
    ap.add_argument("--out-html", default=os.path.join(here, "canadian_ais_inventory.html"))
    ap.add_argument("--out-json", default=os.path.join(here, "canadian_ais_inventory.json"))
    ap.add_argument("--cache",    default="/tmp/glansis_ncbi_cache.json",
                    help="Shared NCBI cache; defaults to the GLANSIS cache "
                         "so species overlapping both lists aren't re-queried.")
    ap.add_argument("--refresh",  action="store_true")
    ap.add_argument("--delay-ms", type=int, default=120)
    ap.add_argument("--limit",    type=int, default=0)
    args = ap.parse_args()

    canada_entries = gi.load_canada_schedule(args.canada_schedule)
    species = gi.species_from_canada_schedule(canada_entries)
    if args.limit:
        species = species[: args.limit]
    n_part2 = sum(1 for v in canada_entries.values() if v["part"] == "2")
    n_part3 = sum(1 for v in canada_entries.values() if v["part"] == "3")
    print(f"[canada] SOR/2015-121 species: {len(species)} "
          f"(Part 2: {n_part2}, Part 3: {n_part3})", file=sys.stderr)

    ra_enrich = gi.load_ra_file(args.ra)
    regs_enrich = gi.load_regs(args.regs)
    regs_enrich = gi.merge_canada_into_regs(regs_enrich, canada_entries)

    results, counts = gi.classify_species(
        species, ra_enrich, regs_enrich,
        cache_path=args.cache, delay_ms=args.delay_ms, refresh=args.refresh,
    )

    intro = (
        f'Generated {time.strftime("%Y-%m-%d")} &middot; '
        f'{len(results)} species from the '
        f'<a href="{gi.CANADA_AIS_REGULATION_URL}">'
        f'federal Aquatic Invasive Species Regulations (SOR/2015-121)</a> '
        f'Schedule, Parts 2 ({n_part2} prohibited species) and '
        f'3 ({n_part3} species at risk) &middot; '
        f'per-species lookup via the '
        f'<a href="https://api.ncbi.nlm.nih.gov/datasets/">NCBI Datasets v2 API</a> '
        f'and <a href="https://eutils.ncbi.nlm.nih.gov/">Entrez E-utilities</a>. '
        f'Note: Manitoba\'s MR 173/2015 Schedule A is a provincial-enforcement '
        f'mirror of this list, not an independent source.'
    )

    gi.write_outputs(
        results, counts,
        out_tsv=args.out_tsv, out_json=args.out_json, out_html=args.out_html,
        source_name="Canada SOR/2015-121",
        html_title="Canadian AIS Regulation (SOR/2015-121) &rarr; NCBI assembly inventory",
        html_intro_html=intro,
    )

    print(
        f"[canada] done. Tier A={counts['A']}  B={counts['B']}  "
        f"C={counts['C']}  D={counts['D']}  E={counts['E']}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
