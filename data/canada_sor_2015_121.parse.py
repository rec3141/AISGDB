#!/usr/bin/env python3
"""Parse Schedule Part 2 of SOR/2015-121 from pdftohtml -xml dump.

Layout quirk worth knowing: for each entry, the scientific name's GENUS
sits on the row *above* the item number, and the species on the row
below. So a row-grouping parser misattributes the genus to the previous
entry. We side-step that by treating item numbers as anchors and binding
each text block to the closest anchor on the same page (within ±25 px).

English Part 2 spans pages 19–23; PARTIE 2 starts page 24."""

import re
import sys
import xml.etree.ElementTree as ET

SRC = "/data/scratch/sor_2015_121/out.xml"
PART2_PAGES = range(19, 24)
PART3_PAGES = range(28, 30)
# On the first Part 2 page (page 19), Part 1 / Definitions content sits
# above this y-coordinate. The "PART 2" heading is at top=666; the table
# header begins at top=750. Skip anything above 700 on page 19 to drop
# definition rows whose item numbers ("1", "2") would otherwise be mistaken
# for entry anchors.
PART2_START_TOP_ON_FIRST_PAGE = 700
# Part 3 also starts on its first page below a heading; skip above this.
PART3_START_TOP_ON_FIRST_PAGE = 720

# Part 3 has its own column layout (only 3 columns: Item, Common, Sci) at
# different x-coordinates from Part 2.
PART3_COL_BINS = [
    ("item",        70, 130),
    ("common_name", 130, 270),
    ("scientific",  270, 460),
]
PART3_COL_NAMES = [c[0] for c in PART3_COL_BINS]

# Hand-patched typos. The anchor-based binder picks up stray glyphs at
# ±20 px boundaries — easier to correct two known cases than to chase
# the geometry.
TYPO_FIXES = {
    "New Zealand ud snail":   "New Zealand mud snail",
    "Orconectes rusticusi":   "Orconectes rusticus",
}

# Column bins by `left` coordinate. The Importation-column continuation
# wraps with a slight indent at left ~399, so the bin extends from 375.
COL_BINS = [
    ("item",            70, 100),
    ("common_name",    100, 185),
    ("scientific",     185, 300),
    ("condition",      300, 375),
    ("import_area",    375, 510),
    ("possess_area",   510, 615),
    ("transport_area", 615, 740),
    ("release_area",   740, 900),
]
COL_NAMES = [c[0] for c in COL_BINS]

HEADER_TEXTS = {
    "Item", "Name", "Common", "Scientific Name", "Condition",
    "Area", "Prohibited", "Area — Importation", "Area — Possession",
    "Area — Transportation", "Area — Release",
    "Column 1", "Column 2", "Column 3", "Column 4",
    "Column 5", "Column 6", "Column 7",
    "Species Subject to Prohibitions and Controls",
    "SCHEDULE", "ANNEXE",
    "Aquatic Invasive Species Regulations",
    "Règlement sur les espèces aquatiques envahissantes",
}
PAGE_FOOTER_RE = re.compile(
    r"^(Current to\b|Last amended\b|À jour au\b|Dernière modification\b)"
)


def bin_text(left):
    for name, lo, hi in COL_BINS:
        if lo <= left < hi:
            return name
    return None


def main():
    tree = ET.parse(SRC)

    # Pass 1: collect all text cells from the Part 2 pages, in (page, top, left, text) form.
    cells = []
    for page in tree.getroot().findall("page"):
        pnum = int(page.get("number"))
        if pnum not in PART2_PAGES:
            continue
        for t in page.findall("text"):
            top   = int(t.get("top"))
            left  = int(t.get("left"))
            text  = "".join(t.itertext()).strip()
            if not text or text in HEADER_TEXTS or PAGE_FOOTER_RE.match(text):
                continue
            if pnum == min(PART2_PAGES) and top < PART2_START_TOP_ON_FIRST_PAGE:
                continue   # Part 1 / definitions on page 19
            if top > 1100:
                continue   # page-footer band (page numbers, "Current to" lines)
            cells.append((pnum, top, left, text))

    # Pass 2: find item-number anchors. An anchor is an item-column integer.
    anchors = []   # list of (pnum, top, item_number)
    for pnum, top, left, text in cells:
        if bin_text(left) == "item" and re.fullmatch(r"\d+", text) and len(text) <= 3:
            anchors.append((pnum, top, text))
    anchors.sort()
    print(f"Found {len(anchors)} item-number anchors", file=sys.stderr)

    # Pass 3: bind each cell to the nearest anchor on the same page,
    # within a vertical tolerance. Anchor's top is the row of the item
    # number itself; genus appears ~13 px above; species ~11 px below.
    # Use ±25 px to capture both directions safely.
    by_item = {}     # (pnum, anchor_top) → dict[col]→list of (top, text)
    for pnum, top, anchor_top, _ in [
        (p, t, at, item) for p, t, _, _ in cells
        for ap, at, item in anchors if False  # placeholder; rewritten below
    ]:
        pass

    # Index anchors per page for fast lookup.
    page_anchors = {}
    for ap, at, item in anchors:
        page_anchors.setdefault(ap, []).append((at, item))
    for ap in page_anchors:
        page_anchors[ap].sort()

    def assign(pnum, top):
        """Find the entry-anchor this text belongs to.
        - Cells within 20 px ABOVE the next anchor → that anchor (genus row).
        - Otherwise → the most recent anchor at or above (item row + continuation).
        """
        ap_list = page_anchors.get(pnum, [])
        if not ap_list:
            return None
        prev = next_a = None
        for at, item in ap_list:
            if at <= top:
                prev = (at, item)
            elif next_a is None:
                next_a = (at, item)
                break
        if next_a and (next_a[0] - top) <= 20:
            return (pnum, next_a[0], next_a[1])
        if prev:
            return (pnum, prev[0], prev[1])
        return None

    entries = {}   # (pnum, anchor_top) → entry dict
    for pnum, top, left, text in cells:
        anchor = assign(pnum, top)
        if anchor is None:
            continue
        key = (anchor[0], anchor[1])
        item_num = anchor[2]
        e = entries.setdefault(key, {"item": item_num,
                                     **{c: "" for c in COL_NAMES[1:]}})
        col = bin_text(left)
        if col is None or col == "item":
            continue
        sep = " " if e[col] else ""
        e[col] += sep + text

    # Sort entries by (page, anchor_top).
    ordered = [entries[k] for k in sorted(entries.keys())]

    # Cleanup whitespace.
    for e in ordered:
        for k in e:
            e[k] = re.sub(r"\s+", " ", e[k]).strip()

    # ---------------- Part 3 ----------------
    part3_entries = []
    for page in tree.getroot().findall("page"):
        pnum = int(page.get("number"))
        if pnum not in PART3_PAGES:
            continue
        rows = {}
        for t in page.findall("text"):
            top  = int(t.get("top"))
            left = int(t.get("left"))
            text = "".join(t.itertext()).strip()
            if not text or text in HEADER_TEXTS or PAGE_FOOTER_RE.match(text):
                continue
            if pnum == min(PART3_PAGES) and top < PART3_START_TOP_ON_FIRST_PAGE:
                continue
            if left >= 470:
                continue  # French side on bilingual Part 3 pages
            rows.setdefault(round(top / 5) * 5, []).append((left, text))
        for top in sorted(rows):
            cells = sorted(rows[top])
            entry = {n: "" for n in PART3_COL_NAMES}
            for left, text in cells:
                for name, lo, hi in PART3_COL_BINS:
                    if lo <= left < hi:
                        sep = " " if entry[name] else ""
                        entry[name] += sep + text
                        break
            if entry["item"].strip().isdigit() and entry["scientific"]:
                part3_entries.append({**entry, "part": "3",
                                      "condition": "", "import_area": "",
                                      "possess_area": "", "transport_area": "",
                                      "release_area": ""})

    # Apply typo fixes (Part 2 only — Part 3 was clean).
    for e in ordered:
        e["part"] = "2"
        e["common_name"] = TYPO_FIXES.get(e["common_name"], e["common_name"])
        e["scientific"]  = TYPO_FIXES.get(e["scientific"],  e["scientific"])

    all_entries = ordered + part3_entries
    cols = ["part"] + COL_NAMES
    print("\t".join(cols))
    for e in all_entries:
        print("\t".join(e[c] for c in cols))

    print(f"\nPARSED Part 2: {len(ordered)} entries", file=sys.stderr)
    print(f"PARSED Part 3: {len(part3_entries)} entries", file=sys.stderr)


if __name__ == "__main__":
    main()
