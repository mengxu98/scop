"""Compare cached portal result tables with local NK demonstration output.

Run from the SCOP checkout: python inst/irea/irea_compare.py D:/scop-dev/datasets/IREA
The portal tables are obtained using its public example input and download button.
"""
import csv
import re
import sys
from pathlib import Path


root = Path(sys.argv[1])
rows = []


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def local_term(portal_term):
    aliases = {
        "4-1BBL": "41BBL", "ADSF": "Resistin", "AdipoQ": "Adiponectin",
        "CT-1": "Cardiotrophin-1", "FGF-��": "FGF-basic",
        "FLT3L": "Flt3l", "IGF-1": "IGF-I", "IL-1Ra": "IL1ra",
        "IL-36Ra": "IL36RA", "IL-Y": "IL-Y", "NP": "Neuropoietin",
        "PRL": "Prolactin", "PSPN": "Persephin",
        "TGF-��1": "TGF-beta-1", "TNF-��": "TNFa",
    }
    if portal_term in aliases:
        return aliases[portal_term]
    if portal_term.startswith("IL-"):
        return portal_term.replace("-", "")
    return portal_term


for analysis in ("cytokine_response", "cell_polarization"):
    for method in ("score", "hypergeometric"):
        stem = f"NK_mouse_genelist_{analysis}_{method}"
        portal = read(root / f"portal_{stem}.csv")
        local = {r["term"]: r for r in read(root / f"scop_{stem}.csv")}
        for original in portal:
            label = original.get("Cytokine", original.get("Polarization"))
            term = local_term(label)
            match = local.get(term)
            rows.append(dict(input="gene_list", analysis=analysis, method=method,
                portal_term=label, local_term=term if match else "",
                portal_effect=original["Enrichment Score"],
                local_effect=match["effect"] if match else "",
                portal_fdr=original["padj"],
                local_fdr=match["fdr"] if match else "",
                status="compared" if match else "term_alias_unresolved"))

for analysis in ("cytokine_response", "cell_polarization"):
    for column, index in (("aPD1_to_Ctrl", 0), ("aCTLA4_to_Ctrl", 1)):
        portal = read(root / f"portal_NK_mouse_example_matrix_{analysis}_{index}.csv")
        local = {r["term"]: r for r in read(root / f"scop_NK_mouse_matrix_{column}_{analysis}.csv")}
        for original in portal:
            label = original.get("Cytokine", original.get("Polarization"))
            term = local_term(label)
            match = local.get(term)
            rows.append(dict(input=column, analysis=analysis, method="projection",
                portal_term=label, local_term=term if match else "",
                portal_effect=original["Enrichment Score"],
                local_effect=match["effect"] if match else "",
                portal_fdr=original["padj"],
                local_fdr=match["fdr"] if match else "",
                status="compared" if match else "term_alias_unresolved"))

out = root / "numeric_comparison.csv"
with out.open("w", encoding="utf-8", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

comparable = [r for r in rows if r["status"] == "compared"]
close = [r for r in comparable if abs(float(r["portal_effect"])-float(r["local_effect"])) < 1e-8]
print(f"Compared {len(comparable)}/{len(rows)} terms; exact effect matches: {len(close)}.")
print(f"Report: {out}")
