#!/usr/bin/env python3
"""Parse Schedule A of Manitoba's Aquatic Invasive Species Regulation (M.R.
173/2015) from a pdftotext -layout dump and emit a clean TSV.

Strategy: split each line on 2+ consecutive spaces to get cell boundaries
(the PDF uses ≥3 spaces between columns and ≤1 within a cell). Stitches
wrapped continuation lines back onto the preceding entry. Forward-fills
the General Classification, which appears only on the first row of each
group (Fish / Invertebrates / Plants / Algae)."""

import re
import sys

src = "/data/scratch/manitoba_ais/reg.txt"

# Drop these recurring lines without trying to parse them.
SKIP_RE = re.compile(
    r"^\s*("
    r"WATER PROTECTION|PROTECTION DES EAUX|"
    r"Accessed:|Current from|Insert Date|"
    r"Schedule [ABC]|SCHEDULE [ABC]|"
    r"Column \d|Colonne \d|"
    r"Common name|Scientific name|Condition|"
    r"General|Classification|"
    r"M\.R\.|"
    r"\d+\s*$|"
    r"AQUATIC INVASIVE|ANNEXE|"
    r"\(Section \d+\)"
    r")"
)


def is_new_entry(parts):
    """A new entry has 2+ non-empty cells where cell 2 (common name)
    starts with a capital letter or 'Any '. Plus optional cell 1
    (classification)."""
    if len(parts) < 2:
        return False
    # Find the common-name cell: the first capitalized cell that's not a
    # bare classification.
    for i, p in enumerate(parts):
        p = p.strip()
        if not p:
            continue
        if p in ("Fish", "Invertebrates", "Plants", "Algae"):
            continue
        return p[0].isupper() or p.startswith("Any ")
    return False


def main():
    lines = open(src).read().splitlines()
    # Restrict to the English Schedule A region.
    start = None
    end = None
    for i, ln in enumerate(lines):
        if ln.strip() == "SCHEDULE A" and start is None:
            start = i
        if ln.strip() in ("ANNEXE A", "M.R. 52/2017") and start is not None:
            end = i
            break
    region = lines[start:end]
    print(f"Schedule A region: lines {start}..{end} ({end-start} lines)", file=sys.stderr)

    rows = []
    current_class = ""
    current_entry = None

    for ln in region:
        if not ln.strip():
            continue
        if SKIP_RE.match(ln):
            continue

        # Detect a leading classification on this line.
        m = re.match(r"^\s*(Fish|Invertebrates|Plants|Algae)\b", ln)
        if m:
            current_class = m.group(1)
            ln = ln[m.end():]   # strip it so it doesn't pollute the split

        # Split on 2+ spaces.
        parts = [p for p in re.split(r"\s{2,}", ln.strip()) if p]

        if not parts:
            continue

        if is_new_entry(parts):
            # commit previous
            if current_entry:
                rows.append(current_entry)
            # Slot the parts into (common, scientific, condition). The
            # number of cells per line is 1–3 depending on what wrapped.
            common = parts[0] if len(parts) >= 1 else ""
            sci    = parts[1] if len(parts) >= 2 else ""
            cond   = parts[2] if len(parts) >= 3 else ""
            current_entry = {
                "classification":   current_class,
                "common_name":      common,
                "scientific_name":  sci,
                "condition":        cond,
            }
        else:
            # Continuation of the previous entry.
            if not current_entry:
                continue
            # Heuristic: if there's one part, append to whichever cell
            # is least filled. If multiple parts, append in order.
            if len(parts) == 1:
                p = parts[0]
                # Common name continuations are usually short (1-3 words);
                # scientific continuations look Latin (capitalized genus or
                # all-lowercase); conditions are "eviscerated" / "hybrids".
                if p in ("eviscerated", "Dead", "Dead and eviscerated"):
                    current_entry["condition"] += " " + p
                elif re.match(r"^[a-z]", p) or re.match(r"^[A-Z][a-z]+\b", p):
                    # If scientific cell already ends with a comma or genus,
                    # this is its continuation.
                    current_entry["scientific_name"] += " " + p
                else:
                    current_entry["common_name"] += " " + p
            else:
                # Multiple parts — slot left to right by what's missing.
                idx = 0
                for k in ("common_name", "scientific_name", "condition"):
                    if idx >= len(parts):
                        break
                    if not current_entry[k]:
                        current_entry[k] = parts[idx]
                        idx += 1
                    else:
                        current_entry[k] += " " + parts[idx]
                        idx += 1

    if current_entry:
        rows.append(current_entry)

    # Light cleanup.
    for r in rows:
        for k in r:
            r[k] = re.sub(r"\s+", " ", r[k]).strip()

    print("classification\tcommon_name\tscientific_name\tcondition")
    for r in rows:
        print(f"{r['classification']}\t{r['common_name']}\t"
              f"{r['scientific_name']}\t{r['condition']}")

    from collections import Counter
    cc = Counter(r["classification"] for r in rows)
    print(f"\nROWS: {len(rows)}  |  by class: {dict(cc)}", file=sys.stderr)


if __name__ == "__main__":
    main()
