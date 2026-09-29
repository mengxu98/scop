"""Cache public IREA example output for numerical comparison.

Usage: python inst/irea/irea_capture_portal.py D:/scop-dev/datasets/IREA
Requires requests and beautifulsoup4; this script is not needed by SCOP itself.
"""
import csv
import re
import sys
from pathlib import Path

import requests
from bs4 import BeautifulSoup


dest = Path(sys.argv[1])
dest.mkdir(parents=True, exist_ok=True)
url = "https://www.immune-dictionary.org/irea/"
session = requests.Session()
page = session.get(url, timeout=60)
page.raise_for_status()
token = re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', page.text).group(1)
genes = "\n".join(("Ifitm3", "Isg15", "Ifit3", "Bst2", "Slfn5",
                   "Isg20", "Phf11b", "Zbp1", "Rtp4"))


def write_table(table, stem):
    table = BeautifulSoup(str(table), "html.parser")
    rows = [[cell.get("data-sort-value", cell.get_text(" ", strip=True))
             for cell in tr.find_all(("th", "td"))]
            for tr in table.select("tr")]
    with (dest / f"{stem}.csv").open("w", encoding="utf-8", newline="") as stream:
        csv.writer(stream).writerows(rows)


for analysis in ("cytokine_response", "cell_polarization"):
    for method in ("score", "hypergeometric"):
        data = dict(csrfmiddlewaretoken=token, irea_form="gene-list", genes=genes,
                    Celltype="NK_cell", Species="Mouse", gene_list_method=method,
                    analysis=analysis)
        response = session.post(url, data=data, timeout=600)
        response.raise_for_status()
        stem = f"portal_NK_mouse_genelist_{analysis}_{method}"
        (dest / f"{stem}.html").write_bytes(response.content)
        table = BeautifulSoup(response.content, "html.parser").select_one("table")
        if table is None:
            raise RuntimeError(f"No results table in {stem}")
        write_table(table, stem)

    data = dict(csrfmiddlewaretoken=token, irea_form="gene-matrix",
                use_example="true", Celltype="NK_cell", Species="Mouse",
                gene_diff_cutoff="0.25", analysis=analysis)
    response = session.post(url, data=data, timeout=600)
    response.raise_for_status()
    stem = f"portal_NK_mouse_example_matrix_{analysis}"
    (dest / f"{stem}.html").write_bytes(response.content)
    soup = BeautifulSoup(response.content, "html.parser")
    templates = soup.select("template[data-irea-sample-template]")
    if not templates:
        raise RuntimeError(f"No matrix sample results in {stem}")
    for template in templates:
        table = template.select_one("table")
        if table is None:
            raise RuntimeError(f"No table in {stem} sample {template['data-irea-sample-template']}")
        write_table(table, f"{stem}_{template['data-irea-sample-template']}")
