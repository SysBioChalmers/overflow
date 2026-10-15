"""Figures and tables: this study next to the Xia et al. (2022) chemostat series."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from overflow.plots import FONTS, LEGACY_COLOURS, capacity_usage_figure, _POINTS_PER_SIZE

HERE = Path(__file__).resolve().parent
FIG, TAB, OUT = HERE / "figures", HERE / "tables", HERE / "respiratory_limit"
REPO = HERE.parents[2]
XIA = HERE.parent
SETS = {
    "ThisStudy": dict(root=REPO, conds=["CN4", "CN22", "CN38", "CN75", "hGR"],
                      label="This study", colour=LEGACY_COLOURS[1], marker="o"),
    "Xia": dict(root=XIA, conds=["mu010", "mu028", "mu033"],
                label="Xia et al. 2022", colour=LEGACY_COLOURS[2], marker="^"),
}
DISPLAY = {"mu010": "μ0.10", "mu028": "μ0.28", "mu033": "μ0.33"}
ORDER = [c for s in SETS.values() for c in s["conds"]]
SYSTEMS = ("Glycolysis", "TCA cycle", "ETC", "Ribosome")
SUPPLEMENT = ("PP shunt", "THF cycle", "Nitrogen metabolism", "Amino acid metabolism")
EXCHANGES = {"Glucose": "r_1714", "CO2": "r_1672", "Oxygen": "r_1992", "Ethanol": "r_1761"}
TEXT, SMALL, LINE = 7.0, 5.6, 0.15 * _POINTS_PER_SIZE
STYLE = {"font.family": FONTS, "pdf.fonttype": 42, "font.size": SMALL, "axes.linewidth": LINE,
         "xtick.major.width": LINE, "ytick.major.width": LINE, "xtick.major.size": 2.75,
         "ytick.major.size": 2.75, "axes.spines.top": False, "axes.spines.right": False,
         "legend.frameon": False, "axes.labelsize": TEXT, "axes.titlesize": SMALL}


def name(c):
    return DISPLAY.get(c, c)


def save(fig, stem):
    fig.savefig(FIG / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def read(dataset, *parts):
    return pd.read_csv(SETS[dataset]["root"].joinpath(*parts), sep="\t")


def fermentation():
    rows = []
    for ds, s in SETS.items():
        f = read(ds, "data", "fermentationData.txt")
        f.insert(0, "dataset", ds)
        rows.append(f[f["Condition"].isin(s["conds"])])
    return pd.concat(rows).rename(columns={"Condition": "condition"}).set_index("condition")


FERM = fermentation()


def overflowing(c):
    return pd.notna(FERM.loc[c, "Ethanol"]) and FERM.loc[c, "Ethanol"] > 0


# ---------------------------------------------------------------- tables

def table_conditions():
    rows = []
    for ds, s in SETS.items():
        for c in s["conds"]:
            flux = read(ds, "results", "modelSimulation", f"allFluxes_{c}.tsv").set_index("rxnID")["flux"]
            gen = lambda kind: s["root"] / "results" / "modelGeneration" / f"{kind}_{c}.tsv"
            relaxed = pd.read_csv(gen("relaxedEnzymes"), sep="\t") if gen("relaxedEnzymes").exists() else pd.DataFrame(columns=["added_mg_gDW"])
            row = dict(dataset=SETS[ds]["label"], condition=name(c), D=FERM.loc[c, "D"], Ptot=FERM.loc[c, "Ptot"])
            for ex, rid in EXCHANGES.items():
                measured = FERM.loc[c, ex]
                row[f"{ex} meas"] = 0.0 if pd.isna(measured) else measured
                row[f"{ex} model"] = abs(flux.get(rid, 0.0))
            row["caps raised to minimum"] = len(pd.read_csv(gen("minimumUsage"), sep="\t"))
            row["caps relaxed for growth"] = len(pd.read_csv(gen("modifiedEnzymes"), sep="\t"))
            row["caps relaxed for rates"] = len(relaxed)
            row["mg/gDW added for rates"] = relaxed["added_mg_gDW"].sum()
            rows.append(row)
    return pd.DataFrame(rows)


def table_usage(r):
    """Capacity usage per system: enzyme mass used over the measured capacity of its enzymes (%)."""
    r = r.assign(usage=100 * r["used"] / r["cap"].where(r["cap"] > 0))
    wide = r.pivot(index="system", columns="condition", values="usage").reindex(list(SYSTEMS) + list(SUPPLEMENT))[ORDER]
    return wide.rename(columns=DISPLAY).reset_index()


def selected_fluxes():
    k = read("ThisStudy", "results", "randomSampling", "selectedFluxes.tsv").set_index("Row")
    x = read("Xia", "results", "randomSampling", "selectedFluxes.tsv").set_index("Row")
    return k.join(x)


def table_atp():
    f = selected_fluxes()
    rows = {
        "glucose uptake": f.loc["rGlu"],
        "ATP from ETC": f.loc["ETC_rATP"],
        "ATP from glycolysis": f.loc["glycolysis_rATP"],
        "ETC share of ATP supply (%)": 100 * f.loc["ETC_rATP"] / (f.loc["ETC_rATP"] + f.loc["glycolysis_rATP"]),
        "ATP per glucose via ETC": f.loc["ETC_YATP_glu"],
        "ATP for growth (GAEC)": f.loc["GAEC_rATP"],
        "ATP not tied to growth (NGAM)": f.loc["NGAM_rATP"],
        "ATP for other metabolism": f.loc["Metabolism_rATP"],
        "pyruvate dehydrogenase": f.loc["rPDH"],
        "isocitrate dehydrogenase (mito)": f.loc["rIDH"],
        "external NADH dehydrogenase": f.loc["rNDE"],
        "NADPH turnover (total)": f.loc["rNADP[tot]"],
    }
    return pd.DataFrame(rows).T[ORDER].rename(columns=DISPLAY).reset_index().rename(columns={"index": "flux (mmol/gDW/h)"})


def system_mass(usage, conditions):
    """Per system and condition: enzyme mass used and the measured cap of its enzymes, in mg/gDW."""
    annotation = read("ThisStudy", "data", "selectedAnnotation.txt")
    system = dict(zip(annotation.iloc[:, 0], annotation.iloc[:, 3]))
    usage = usage.assign(system=usage["protID"].map(system))
    rows = []
    for c in conditions:
        measured = usage[usage[f"UB_{c}"] < 1000]  # unmeasured enzymes draw on the pool, cap 1000
        for sysname, g in measured.groupby("system"):
            rows.append(dict(condition=c, system=sysname, used=g[f"absUse_{c}"].sum(), cap=g[f"UB_{c}"].sum()))
    return pd.DataFrame(rows)


def reserve():
    return pd.concat(system_mass(read(ds, "results", "enzymeUsage", "enzymeUsages.tsv"), s["conds"]).assign(dataset=ds)
                     for ds, s in SETS.items())


def table_reserve(r):
    """Enzyme mass used and measured capacity per system (mg/gDW), the numbers of figure 3."""
    r = r[r["system"].isin(SYSTEMS)].rename(columns={"used": "used", "cap": "measured capacity"})
    long = r.melt(id_vars=["system", "condition"], value_vars=["used", "measured capacity"], var_name="mg/gDW")
    wide = long.pivot(index=["system", "mg/gDW"], columns="condition", values="value")[ORDER]
    wide = wide.reindex(pd.MultiIndex.from_product([list(SYSTEMS), ["used", "measured capacity"]]))
    return wide.rename(columns=DISPLAY).rename_axis(["system", "mg/gDW"]).reset_index()


def table_kcats():
    """The respiratory-chain kcats each dataset's build uses, against the ecModel's."""
    rows = []
    for ds, s in SETS.items():
        k = read(ds, "results", "modelGeneration", "respiratoryKcats.tsv")
        k.insert(0, "dataset", s["label"])
        k["condition"] = k["condition"].map(name)
        rows.append(k)
    return pd.concat(rows).rename(columns={"condition": "set by", "kcat_per_s": "derived kcat (1/s)",
                                           "ecmodel_kcat_per_s": "ecModel kcat (1/s)"})


def limit_table():
    files = sorted(OUT.glob("limit_*.tsv"))
    if not files:
        return None
    t = pd.concat(pd.read_csv(f, sep="\t") for f in files)
    t["condition"] = pd.Categorical(t["condition"], ORDER, ordered=True)
    return t.sort_values(["variant", "condition"])


# ---------------------------------------------------------------- figures

def figure_usage(systems, stem):
    k = read("ThisStudy", "results", "enzymeUsage", "capUsage.tsv")
    x = read("Xia", "results", "enzymeUsage", "capUsage.tsv")
    both = k.drop(columns=["GOterm"]).merge(x[["protID", "mu010", "mu028", "mu033"]], on="protID", how="outer")
    both["GOterm"] = both["protID"].map(dict(zip(pd.concat([k, x])["protID"], pd.concat([k, x])["GOterm"])))
    fig = capacity_usage_figure(both, systems=systems, conditions=ORDER, width=14 / 2.54, height=5 / 2.54)
    for axis in fig.axes:
        axis.axvline(5.5, color="#878787", linewidth=LINE, linestyle=(0, (2, 2)))
        axis.set_xticklabels([name(c) for c in ORDER])
    fig.text(0.985, 0.97, "CN: this study   μ: Xia et al. 2022", ha="right", va="top", fontsize=SMALL, color="#4D4D4D")
    save(fig, stem)


def figure_trend(usage):
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 4, figsize=(16 / 2.54, 8 / 2.54), sharex=True)
        for axis, system in zip(axes.flat, list(SYSTEMS) + list(SUPPLEMENT)):
            row = usage.set_index("system").loc[system]
            for ds, s in SETS.items():
                xs = [FERM.loc[c, "Glucose"] for c in s["conds"]]
                ys = [row[name(c)] for c in s["conds"]]
                order = np.argsort(xs)
                axis.plot(np.array(xs)[order], np.array(ys)[order], color=s["colour"], linewidth=0.6)
                for c, xv, yv in zip(s["conds"], xs, ys):
                    axis.plot(xv, yv, s["marker"], markersize=3.5, color=s["colour"],
                              markerfacecolor=s["colour"] if overflowing(c) else "white", markeredgewidth=0.6)
            axis.set_title(system, color="#1A1A1A")
            axis.set_ylim(-5, 105)
            axis.set_yticks([0, 50, 100])
        for axis in axes[1]:
            axis.set_xlabel("Glucose uptake (mmol/gDW/h)", fontsize=SMALL)
        for axis in axes[:, 0]:
            axis.set_ylabel("Capacity usage (%)", fontsize=SMALL)
        handles = [plt.Line2D([], [], color=s["colour"], marker=s["marker"], markersize=3.5, linewidth=0.6, label=s["label"]) for s in SETS.values()]
        handles.append(plt.Line2D([], [], color="#4D4D4D", marker="o", markersize=3.5, linestyle="none", markerfacecolor="#4D4D4D", label="ethanol secreted"))
        handles.append(plt.Line2D([], [], color="#4D4D4D", marker="o", markersize=3.5, linestyle="none", markerfacecolor="white", label="no ethanol"))
        fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=SMALL, bbox_to_anchor=(0.5, 0.045))
        fig.text(0.5, 0.005, "capacity usage: enzyme mass used over the measured capacity of the system's enzymes",
                 ha="center", va="bottom", fontsize=SMALL - 0.8, color="#4D4D4D")
        fig.subplots_adjust(left=0.08, right=0.99, top=0.93, bottom=0.26, hspace=0.45, wspace=0.3)
        save(fig, "fig2_usage_vs_uptake")


def figure_reserve(r):
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 4, figsize=(16 / 2.54, 5.5 / 2.54))
        for index, (axis, system) in enumerate(zip(axes, SYSTEMS)):
            sub = r[r["system"] == system].set_index("condition").loc[ORDER]
            xs = np.arange(len(ORDER)) + np.array([0] * 5 + [0.6] * 3)
            colour = LEGACY_COLOURS[index]
            axis.bar(xs, sub["cap"], width=0.75, color="white", edgecolor=colour, linewidth=0.35 * _POINTS_PER_SIZE)
            axis.bar(xs, sub["used"], width=0.75, color=colour, edgecolor=colour, linewidth=0.35 * _POINTS_PER_SIZE)
            axis.set_xticks(xs)
            axis.set_xticklabels([name(c) for c in ORDER], rotation=90)
            axis.set_title(system, color="#1A1A1A")
        axes[0].set_ylabel("Enzyme (mg/gDW)", fontsize=TEXT)
        fig.text(0.99, 0.01, "filled: used in the simulation   outline: measured capacity", ha="right", fontsize=SMALL, color="#4D4D4D")
        fig.subplots_adjust(left=0.07, right=0.99, top=0.9, bottom=0.3, wspace=0.35)
        save(fig, "fig3_proteome_reserve")


def figure_atp():
    f = selected_fluxes()[ORDER]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(16 / 2.54, 5.5 / 2.54), gridspec_kw=dict(width_ratios=[1, 1, 0.9]))
        xs = np.arange(len(ORDER)) + np.array([0] * 5 + [0.6] * 3)
        labels = [name(c) for c in ORDER]
        supply = [("ETC", f.loc["ETC_rATP"], LEGACY_COLOURS[2]), ("glycolysis", f.loc["glycolysis_rATP"], LEGACY_COLOURS[0])]
        demand = [("growth (GAEC)", f.loc["GAEC_rATP"], LEGACY_COLOURS[1]), ("other metabolism", f.loc["Metabolism_rATP"], LEGACY_COLOURS[3]),
                  ("not tied to growth", f.loc["NGAM_rATP"], "#E6E0D4")]
        for axis, parts, title in ((axes[0], supply, "ATP supply"), (axes[1], demand, "ATP demand")):
            bottom = np.zeros(len(ORDER))
            for label, values, colour in parts:
                axis.bar(xs, values.to_numpy(), bottom=bottom, width=0.75, color=colour, label=label, edgecolor="white", linewidth=0.2)
                bottom += values.to_numpy()
            axis.set_xticks(xs)
            axis.set_xticklabels(labels, rotation=90)
            axis.set_title(title, color="#1A1A1A")
            axis.legend(fontsize=SMALL - 0.6, loc="upper left", handlelength=1)
        axes[0].set_ylabel("ATP (mmol/gDW/h)", fontsize=TEXT)
        axis = axes[2]
        share = 100 * f.loc["ETC_rATP"] / (f.loc["ETC_rATP"] + f.loc["glycolysis_rATP"])
        for ds, s in SETS.items():
            xv = [FERM.loc[c, "Glucose"] for c in s["conds"]]
            axis.plot(xv, share[s["conds"]], color=s["colour"], linewidth=0.6)
            for c, x in zip(s["conds"], xv):
                axis.plot(x, share[c], s["marker"], markersize=3.5, color=s["colour"],
                          markerfacecolor=s["colour"] if overflowing(c) else "white", markeredgewidth=0.6)
        axis.set_xlabel("Glucose uptake (mmol/gDW/h)", fontsize=SMALL)
        axis.set_ylabel("ATP from ETC (%)", fontsize=SMALL)
        axis.set_title("Respiratory share", color="#1A1A1A")
        axis.set_ylim(60, 100)
        handles = [plt.Line2D([], [], color=s["colour"], marker=s["marker"], markersize=3.5, linewidth=0.6, label=s["label"]) for s in SETS.values()]
        axis.legend(handles=handles, fontsize=SMALL - 0.6, loc="lower left")
        fig.subplots_adjust(left=0.07, right=0.99, top=0.9, bottom=0.27, wspace=0.4)
        save(fig, "fig4_atp_budget")


def figure_limit(limits, variant):
    curves = pd.concat(pd.read_csv(f, sep="\t") for f in OUT.glob(f"curve_*_{variant}.tsv") if f.stat().st_size > 1)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(14 / 2.54, 5.5 / 2.54), sharey=True)
        palette = dict(zip(ORDER, ["#CBBBA0", "#A08F6E", "#6E5F45", "#1D1D1B", "#878787", "#9CC3E6", "#1D71B8", "#0B3B66"]))
        limit = limits[limits["variant"] == variant].set_index("condition")["respiratory_limit"]
        axes[0].axvspan(1.7, 2.3, color="#E6E0D4", linewidth=0, zorder=0)
        axes[0].text(2.35, 13.5, "reported onset of overflow", ha="left", va="center", fontsize=SMALL - 0.6, color="#6E5F45")
        for axis, (ds, s) in zip(axes, SETS.items()):
            for c in s["conds"]:
                cur = curves[(curves["dataset"] == ds) & (curves["condition"] == c)].sort_values("glucose")
                axis.plot(cur["glucose"], cur["forced_overflow_C"], color=palette[c], linewidth=0.8, label=name(c))
                if np.isfinite(limit.get(c, np.nan)):
                    axis.plot(limit[c], -0.6, marker="^", markersize=3, color=palette[c], clip_on=False)
                measured_c = sum(w * (0 if pd.isna(FERM.loc[c, ex]) else FERM.loc[c, ex])
                                 for ex, w in (("Ethanol", 2), ("Acetate", 2), ("Glycerol", 3), ("Formate", 1)))
                axis.plot(FERM.loc[c, "Glucose"], measured_c, "o", markersize=3.5, color=palette[c],
                          markerfacecolor=palette[c] if overflowing(c) else "white", markeredgewidth=0.7)
            axis.set_title(s["label"], color="#1A1A1A")
            axis.set_xlabel("Glucose uptake (mmol/gDW/h)", fontsize=SMALL)
            axis.set_xlim(0.5, 7)
            axis.set_ylim(-1, 14)
            axis.legend(fontsize=SMALL - 0.6, loc="upper left", handlelength=1.2)
        axes[0].set_ylabel("Forced overflow (mmol C/gDW/h)", fontsize=SMALL)
        fig.text(0.99, 0.01, "lines: least overflow the proteome allows at the condition's growth rate   ▲: respiratory limit   points: measured byproduct C",
                 ha="right", fontsize=SMALL - 0.6, color="#4D4D4D")
        fig.subplots_adjust(left=0.09, right=0.99, top=0.9, bottom=0.25, wspace=0.12)
        save(fig, "fig5_respiratory_limit")


def figure_turnover():
    """Apparent turnover of each complex relative to the highest in its dataset; needs no kcat."""
    complexes = ("complex III", "complex IV", "ATP synthase")
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(16 / 2.54, 5.5 / 2.54), sharey=True)
        for axis, cx in zip(axes, complexes):
            for ds, s in SETS.items():
                t = read(ds, "results", "modelGeneration", "respiratoryTurnover.tsv")
                d = t[t["complex"] == cx].set_index("condition").loc[s["conds"]]
                xs = [FERM.loc[c, "Glucose"] for c in s["conds"]]
                axis.plot(xs, d["relative_to_highest"], color=s["colour"], linewidth=0.6)
                for c, x in zip(s["conds"], xs):
                    y = d.loc[c, "relative_to_highest"]
                    axis.plot(x, y, s["marker"], markersize=3.5, color=s["colour"],
                              markerfacecolor=s["colour"] if overflowing(c) else "white", markeredgewidth=0.6)
                    xia = ds == "Xia"
                    axis.annotate(name(c), (x, y), xytext=(-3, 4) if xia else (3, -8), textcoords="offset points",
                                  ha="right" if xia else "left", fontsize=SMALL - 1.4, color=s["colour"])
            axis.set_title(cx[0].upper() + cx[1:], color="#1A1A1A")
            axis.set_ylim(0, 1.15)
            axis.set_xlim(0.4, 6.1)
            axis.set_xlabel("Glucose uptake (mmol/gDW/h)", fontsize=SMALL)
        axes[0].set_ylabel("Apparent turnover / highest in dataset", fontsize=SMALL)
        handles = [plt.Line2D([], [], color=s["colour"], marker=s["marker"], markersize=3.5, linewidth=0.6, label=s["label"])
                   for s in SETS.values()]
        handles += [plt.Line2D([], [], color="#4D4D4D", marker="o", markersize=3.5, linestyle="none", label="ethanol secreted"),
                    plt.Line2D([], [], color="#4D4D4D", marker="o", markersize=3.5, linestyle="none", markerfacecolor="white",
                               label="no ethanol")]
        fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=SMALL, bbox_to_anchor=(0.5, 0.0))
        fig.subplots_adjust(left=0.08, right=0.99, top=0.9, bottom=0.32, wspace=0.12)
        save(fig, "fig6_apparent_turnover")


def markdown(tables):
    lines = []
    for title, t in tables:
        lines += [f"## {title}", "", t.to_markdown(index=False, floatfmt=".3g"), ""]
    (TAB / "tables.md").write_text("\n".join(lines))


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    TAB.mkdir(exist_ok=True)
    r = reserve()
    usage = table_usage(r)
    tables = [
        ("Table 1. Conditions, measured and modelled exchange rates (mmol/gDW/h), and the caps the pipeline changed", table_conditions()),
        ("Table 2. Capacity usage per system: enzyme mass used over the measured capacity of its enzymes (%)", usage),
        ("Table 3. Enzyme mass used and measured capacity per system (mg/gDW)", table_reserve(r)),
        ("Table 4. Sampled ATP and redox budget (mean over samples, mmol/gDW/h)", table_atp()),
    ]
    limits = limit_table()
    if limits is not None:
        tables.append(("Table 5. Respiratory limit: the highest glucose uptake the measured proteome can respire completely at the condition's growth rate (mmol/gDW/h)",
                       limits.replace({"dataset": {k: v["label"] for k, v in SETS.items()}})))
    tables.append(("Table 6. Respiratory-chain kcats derived from each dataset: the highest apparent turnover across its conditions", table_kcats()))
    for i, (title, t) in enumerate(tables, 1):
        t.to_csv(TAB / f"table{i}.tsv", sep="\t", index=False)
    markdown(tables)
    figure_usage(SYSTEMS, "fig1_capacity_usage")
    figure_usage(SUPPLEMENT, "figS1_capacity_usage_supplement")
    figure_trend(usage)
    figure_reserve(r)
    figure_atp()
    if limits is not None:
        figure_limit(limits, "built")
    figure_turnover()
    print(open(TAB / "tables.md").read())
