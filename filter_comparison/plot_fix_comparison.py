"""plot_fix_comparison.py -- paired before/after figure for the synthetic experiments

Reads results_fixed/result-state*/<p>-results.csv (fixed DDMM and all baselines)
and <p>-results_original_core.csv (original DDMM on the same simulated data),
written by compare_new_extend_fixed.R and compare_family_extend_fixed.R with
RUN_ORIGINAL=1. Writes figure_fixed/synthetic_fix_comparison.png and a summary
table figure_fixed/synthetic_fix_summary.csv.

    python3 plot_fix_comparison.py
"""
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figure_fixed")
os.makedirs(OUT, exist_ok=True)


def load(exp):
    rows = []
    for f in glob.glob(os.path.join(HERE, exp, "results_fixed", "result-state*", "*-results*.csv")):
        d = pd.read_csv(f)
        d["state"] = int(re.search(r"state(\d+)", f).group(1))
        d["p"] = int(re.search(r"/(\d+)-results", f).group(1))
        d["arm"] = "original" if "original_core" in f else "fixed"
        rows.append(d)
    d = pd.concat(rows, ignore_index=True)
    d["label"] = d["Method"].str.replace(r"DDMM \((DDMM-[A-Za-z]+), Post-Est\)", r"\1", regex=True)
    return d


NAMES = {"DDMM-E": "DDMM-PCA (E)", "DDMM-O": "DDMM-PCA (O)", "DDMM-Trim": "DDMM-Trim",
         "DDMM-LMM": "DDMM-LMM", "DDMM-PCA": "DDMM-PCA",
         "PCA-E": "PCA-Lasso (E)", "PCA-O": "PCA-Lasso (O)", "LMM": "LMM-Lasso",
         "BH": "LMM-BH", "ggmix": "ggmix", "plmmr": "plmmr"}
DDMM_COL = {"DDMM-PCA (E)": "#882255", "DDMM-PCA (O)": "#009E73", "DDMM-Trim": "#0072B2",
            "DDMM-LMM": "#D55E00", "DDMM-PCA": "#009E73"}
BASE_STYLE = {"PCA-Lasso (E)": ("#E69F00", "^"), "PCA-Lasso (O)": ("#E69F00", "o"),
              "LMM-Lasso": ("#56B4E9", "^"), "LMM-BH": ("#777777", "x"),
              "ggmix": ("#CC79A7", "D"), "plmmr": ("#7A3E9D", "*")}
EXPS = [("admixture", "Low-rank admixture", ["DDMM-PCA (O)", "DDMM-PCA (E)", "DDMM-Trim"],
         ["PCA-Lasso (O)", "PCA-Lasso (E)", "ggmix", "plmmr"]),
        ("family", "Family structure", ["DDMM-Trim", "DDMM-PCA", "DDMM-LMM"],
         ["LMM-Lasso", "LMM-BH", "ggmix", "plmmr"])]

fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))
summary = []
for r, (exp, title, ddmm, base) in enumerate(EXPS):
    d = load(exp)
    d["name"] = d["label"].map(NAMES)
    ns = d["state"].nunique()
    g = (d.groupby(["name", "arm", "p"])[["FDR", "Power"]]
           .agg(["mean", "sem"]).reset_index())
    g.columns = ["name", "arm", "p", "FDR", "FDR_se", "Power", "Power_se"]
    g.insert(0, "experiment", exp)
    g["states"] = ns
    summary.append(g)
    for c, metric in enumerate(["FDR", "Power"]):
        ax = axes[r, c]
        for nm in base:
            s = g[(g.name == nm) & (g.arm == "fixed")].sort_values("p")
            if s.empty:
                continue
            col, mk = BASE_STYLE[nm]
            ax.errorbar(s.p, s[metric], yerr=s[metric + "_se"], color=col, marker=mk, ms=5,
                        lw=1, alpha=0.55, capsize=2, label=nm)
        for nm in ddmm:
            col = DDMM_COL[nm]
            mk = "s" if "Trim" in nm else ("^" if "(E)" in nm else "o")
            for arm, ls, fill, lab in [("original", "--", "none", " original"),
                                       ("fixed", "-", col, " fixed")]:
                s = g[(g.name == nm) & (g.arm == arm)].sort_values("p")
                ax.errorbar(s.p, s[metric], yerr=s[metric + "_se"], color=col, ls=ls, marker=mk,
                            ms=7, mfc=fill, lw=2 if arm == "fixed" else 1.4, capsize=3,
                            label=nm + lab)
        if metric == "FDR":
            ax.axhline(0.1, color="#333333", ls=":", lw=1)
        ax.set_ylim(0, 1)
        ps = sorted(g.p.unique())
        ax.set_xticks(ps)
        ax.set_xticklabels([f"{p:,}" for p in ps])
        ax.set_xlabel("Number of tested variants, p")
        ax.set_ylabel("Empirical FDR" if metric == "FDR" else "Power")
        ax.set_title(f"{title}: {'FDR' if metric == 'FDR' else 'power'}  ({ns} states)",
                     loc="left", fontweight="bold")
        ax.grid(axis="y", color="#e6e6e6")
        ax.spines[["top", "right"]].set_visible(False)
    h, l = axes[r, 0].get_legend_handles_labels()
    axes[r, 1].legend(h, l, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=9)

fig.suptitle("Synthetic experiments: original vs fixed DDMM on the same simulated data\n"
             "dashed/open = original nodewise (one joint lambda), solid/filled = per-variant BIC "
             "nodewise (FIX 1); baselines from the same run; error bars = 1 SE",
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, "synthetic_fix_comparison.png"), dpi=200, bbox_inches="tight")
pd.concat(summary).round(4).to_csv(os.path.join(OUT, "synthetic_fix_summary.csv"), index=False)
print("wrote", os.path.join(OUT, "synthetic_fix_comparison.png"))
