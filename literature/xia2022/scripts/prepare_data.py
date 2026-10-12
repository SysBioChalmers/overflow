"""Write the chemostat series in the overflow pipeline's data formats.

Proteomics in copies/cell are converted per sample to mmol/gDW by scaling the
sample's total protein mass to Ptot:  mmol/gDW_i = 1000 * Ptot * c_i / sum_j(c_j * MW_j).
Protein groups ("A:B") are assigned to the first member with a mass in
data/uniprot.tsv, isoform suffixes stripped. Byproducts not in the flux sheet are NA.

usage: prepare_data.py PROTEOMICS.xlsx FLUXES.xlsx OVERFLOW_REPO OUT_DIR PTOT_010 PTOT_028 PTOT_033
"""
import sys
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

prot_xlsx, flux_xlsx, repo, out = sys.argv[1:5]
ptot = dict(zip(("mu010", "mu028", "mu033"), map(float, sys.argv[5:8])))
out = Path(out)

# condition -> (row of the flux sheet by mu, proteomics sample columns)
CONDITIONS = {
    "mu010": (0.10225, ["S7", "S8", "S9"]),
    "mu028": (0.28425, ["S19", "S20", "S21"]),
    "mu033": (0.33375, ["S22", "S23", "S24"]),
}


def sheet(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    return list(wb.worksheets[0].iter_rows(values_only=True))


# --- fermentation data -------------------------------------------------
rows = sheet(flux_xlsx)
flux = pd.DataFrame(rows[2:], columns=rows[0]).dropna(subset=["mu"])
records = []
for name, (mu, _) in CONDITIONS.items():
    row = flux.loc[np.isclose(flux["mu"].astype(float), mu)]
    assert len(row) == 1, (name, mu)
    row = row.iloc[0]
    etoh = float(row["qEthanol"])
    records.append({
        "Condition": name, "Ptot": ptot[name], "D": float(row["mu"]),
        "Glucose": float(row["qGluc"]), "CO2": float(row["qCO2"]),
        "Oxygen": float(row["qO2"]), "Glycerol": np.nan, "Acetate": np.nan,
        "Ethanol": etoh if etoh > 0 else np.nan, "Formate": np.nan,
    })
pd.DataFrame(records).to_csv(out / "fermentationData.txt", sep="\t", index=False, na_rep="NA")

# --- proteomics --------------------------------------------------------
rows = sheet(prot_xlsx)
prot = pd.DataFrame(rows[2:], columns=rows[0]).dropna(subset=["Majority protein IDs"])
uni = pd.read_csv(Path(repo) / "data" / "uniprot.tsv", sep="\t")
mass = dict(zip(uni["Entry"], uni["Mass"]))
gene = dict(zip(uni["Entry"], uni["Gene Names (ordered locus)"])) if "Gene Names (ordered locus)" in uni else {}


def resolve(group):
    for member in group.split(":"):
        member = member.split("-")[0]
        if member in mass and not pd.isna(mass[member]):
            return member
    return None


prot["id"] = prot["Majority protein IDs"].map(resolve)
unresolved = prot[prot["id"].isna()]
samples = [s for _, cols in CONDITIONS.values() for s in cols]
prot[samples] = prot[samples].apply(pd.to_numeric, errors="coerce")
prot = prot.dropna(subset=["id"])
# isoform groups that resolve onto an already listed protein are summed
n_merged = prot["id"].duplicated().sum()
table = prot.groupby("id", sort=False).agg(
    {"Gene Name": "first", **{s: lambda v: v.sum(min_count=1) for s in samples}}
)
mw = table.index.map(mass).to_numpy(dtype=float)

orf = [gene.get(p) if isinstance(gene.get(p), str) else g for p, g in zip(table.index, table["Gene Name"])]
out_table = pd.DataFrame({"Protein.IDs": table.index, "Gene": [o.split()[0] for o in orf]})
for name, (_, cols) in CONDITIONS.items():
    for k, s in enumerate(cols, start=1):
        copies = table[s].to_numpy(dtype=float)
        total_mg = np.nansum(copies * mw)  # per cell, in copies * mg/mmol
        out_table[f"{name}_{k}_abs"] = 1000.0 * ptot[name] * copies / total_mg
out_table.to_csv(out / "abs_proteomics.txt", sep="\t", index=False, na_rep="NA", float_format="%.6E")

print(f"proteins: {len(rows) - 2} rows, {len(unresolved)} without a mass dropped, "
      f"{n_merged} isoform rows merged, {len(out_table)} written")
if len(unresolved):
    print("dropped:", ", ".join(unresolved["Majority protein IDs"].astype(str)))
