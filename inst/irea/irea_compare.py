"""Compare cached portal/local IREA outputs; no remote request is made.

python inst/irea/irea_compare.py CACHE --local-dir OUTPUT --require-concordance
Optional --analysis and --method select a validation path.
"""
import argparse
import csv
import math
import re
from decimal import Decimal
from pathlib import Path


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def canonical(term):
    aliases = {"ADSF": "Resistin", "AdipoQ": "Adiponectin",
               "CT-1": "Cardiotrophin-1", "FGF-\u03b2": "FGF-basic",
               "IGF-1": "IGF-I", "NP": "Neuropoietin",
               "PRL": "Prolactin", "PSPN": "Persephin"}
    term = aliases.get(term, term).casefold()
    for word, symbol, letter in (("alpha", "\u03b1", "a"), ("beta", "\u03b2", "b"),
                                 ("gamma", "\u03b3", "g"), ("epsilon", "\u03b5", "e"),
                                 ("kappa", "\u03ba", "k"), ("lambda", "\u03bb", "l")):
        term = term.replace(word, letter).replace(symbol, letter)
    return re.sub(r"[^a-z0-9]", "", term)


def numeric_match(portal, local, effect=False):
    value = float(portal)
    if effect:
        return math.isclose(value, float(local), rel_tol=1e-8, abs_tol=1e-8)
    # Older cached tables contain displayed decimals; respect their last
    # printed digit without treating every tiny P value as an absolute zero.
    quantum = 0.0 if value in (0.0, 1.0) else 10.0 ** Decimal(portal).as_tuple().exponent
    return math.isclose(value, float(local), rel_tol=1e-8, abs_tol=0.51 * quantum)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_dir", type=Path)
    parser.add_argument("--local-dir", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--analysis", choices=("cytokine_response", "cell_polarization"))
    parser.add_argument("--method", choices=("score", "hypergeometric", "projection"))
    parser.add_argument("--require-concordance", action="store_true")
    args = parser.parse_args()
    root = args.reference_dir
    local_dir = args.local_dir or root
    rows = []

    def compare(portal_file, local_file, input_name, analysis, method):
        if args.analysis and analysis != args.analysis:
            return
        if args.method and method != args.method:
            return
        local = {}
        for record in read(local_file):
            key = canonical(record["term"])
            if key in local:
                raise ValueError(f"Ambiguous canonical term: {key}")
            local[key] = record
        for original in read(portal_file):
            label = original.get("Cytokine", original.get("Polarization"))
            match = local.get(canonical(label))
            row = dict(input=input_name, analysis=analysis, method=method,
                       portal_term=label, local_term=match["term"] if match else "")
            for short, portal_key, local_key in (("effect", "Enrichment Score", "effect"),
                                                ("p_value", "pval", "p_value"),
                                                ("fdr", "padj", "fdr")):
                row[f"portal_{short}"] = original[portal_key]
                row[f"local_{short}"] = match[local_key] if match else ""
                row[f"{short}_match"] = bool(match) and numeric_match(
                    original[portal_key], match[local_key], effect=short == "effect")
            row["significance_match"] = bool(match) and (
                (float(original["padj"]) < 0.05) == (float(match["fdr"]) < 0.05))
            row["status"] = "term_alias_unresolved" if not match else (
                "concordant" if all(row[f"{key}_match"] for key in
                                    ("effect", "p_value", "fdr")) else "numeric_mismatch")
            rows.append(row)

    for analysis in ("cytokine_response", "cell_polarization"):
        for method in ("score", "hypergeometric"):
            stem = f"NK_mouse_genelist_{analysis}_{method}"
            compare(root / f"portal_{stem}.csv", local_dir / f"scop_{stem}.csv",
                    "gene_list", analysis, method)
        for column, index in (("aPD1_to_Ctrl", 0), ("aCTLA4_to_Ctrl", 1)):
            compare(root / f"portal_NK_mouse_example_matrix_{analysis}_{index}.csv",
                    local_dir / f"scop_NK_mouse_matrix_{column}_{analysis}.csv",
                    column, analysis, "projection")
    out = args.output or local_dir / "numeric_comparison.csv"
    with out.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    matched = [row for row in rows if row["status"] != "term_alias_unresolved"]
    concordant = sum(row["status"] == "concordant" for row in rows)
    print(f"Matched {len(matched)}/{len(rows)} terms; full numeric concordance: {concordant}.")
    print(f"Effect matches: {sum(row['effect_match'] for row in rows)}; "
          f"FDR threshold disagreements: {sum(not row['significance_match'] for row in matched)}.")
    print(f"Report: {out}")
    if args.require_concordance and concordant != len(rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
