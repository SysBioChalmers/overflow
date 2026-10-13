"""Write the series rescaled to the proteome coverage of this study.

The series attributes all of Ptot to measured proteins; in this study the
measured proteins add up to less. This scales every abundance by this
study's mean measured fraction of Ptot, so the unmeasured share is the same.

usage: scale_coverage.py OUT_DIR
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

uniprot = pd.read_csv(REPO / "data" / "uniprot.tsv", sep="\t")
mass = dict(zip(uniprot["Entry"], uniprot["Mass"]))  # g/mol, so mmol/gDW * mass = mg/gDW


def measured_fraction(proteomics: pd.DataFrame, ptot: dict[str, float]) -> pd.Series:
    """Measured protein mass over Ptot, per condition (mean of its replicates)."""
    table = proteomics.drop_duplicates("Protein.IDs").set_index("Protein.IDs")
    mw = table.index.map(mass).to_numpy(dtype=float)
    fractions = {}
    for column in [c for c in table.columns if c.endswith("_abs")]:
        condition = column.rsplit("_", 2)[0].removeprefix("CN_")
        grams = np.nansum(table[column].to_numpy(dtype=float) * mw) / 1000.0
        fractions.setdefault(condition, []).append(grams / ptot[condition])
    return pd.Series({c: float(np.mean(v)) for c, v in fractions.items()})


study_ferm = pd.read_csv(REPO / "data" / "fermentationData.txt", sep="\t")
study = measured_fraction(pd.read_csv(REPO / "data" / "abs_proteomics.txt", sep="\t"),
                         dict(zip(study_ferm["Condition"], study_ferm["Ptot"])))
coverage = float(study.mean())

xia_ferm = pd.read_csv(HERE / "data" / "fermentationData.txt", sep="\t")
xia = pd.read_csv(HERE / "data" / "abs_proteomics.txt", sep="\t")
before = measured_fraction(xia, dict(zip(xia_ferm["Condition"], xia_ferm["Ptot"])))
columns = [c for c in xia.columns if c.endswith("_abs")]
xia[columns] = xia[columns] * coverage
after = measured_fraction(xia, dict(zip(xia_ferm["Condition"], xia_ferm["Ptot"])))

xia.to_csv(out / "abs_proteomics.txt", sep="\t", index=False, na_rep="NA", float_format="%.6E")
for name in ("fermentationData.txt", "ribosome_core.tsv"):
    shutil.copy(HERE / "data" / name, out / name)

print("this study, measured fraction of Ptot:", study.round(3).to_dict())
print(f"coverage applied: {coverage:.4f}")
print("series before:", before.round(3).to_dict(), "after:", after.round(3).to_dict())
