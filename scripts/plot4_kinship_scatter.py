#!/usr/bin/env python3
"""
Plot 4: Pairwise Kinship Scatter — Method vs Method

Three subpanels comparing kinship estimates between each pair of methods.
Only pairs detected by BOTH methods in a given panel are plotted.

Sources:
  PLINK:       KINSHIP from king_3chr_1000.kin0
  GERMLINE2:   kinship_phi_from_P from pairwise_ibd012.csv
  Refined IBD: kinship_phi_from_P from per-chr refinedibd_pairwise_ibd012.csv

Usage:
  python3 scripts/plot4_kinship_scatter.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

PROJDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLOTDIR = os.path.join(PROJDIR, "plots")
os.makedirs(PLOTDIR, exist_ok=True)

# ── Helper: normalize pair to sorted tuple ────────────────────────────────
def norm_pair(a, b):
    return tuple(sorted([str(a), str(b)]))

# ── Load PLINK KING kinship ──────────────────────────────────────────────
king = pd.read_csv(
    os.path.join(PROJDIR, "results", "king_3chr_1000.kin0"),
    sep="\t",
)
plink_kin = {}
for _, row in king.iterrows():
    pair = norm_pair(row["IID1"], row["IID2"])
    plink_kin[pair] = float(row["KINSHIP"])
print(f"PLINK:       {len(plink_kin)} total pairs loaded")

# ── Load GERMLINE2 kinship ───────────────────────────────────────────────
g2 = pd.read_csv(os.path.join(PROJDIR, "results", "GERMLINE2", "pairwise_ibd012.csv"))
g2_kin = {}
for _, row in g2.iterrows():
    pair = norm_pair(row["id1"], row["id2"])
    g2_kin[pair] = float(row["kinship_phi_from_P"])
print(f"GERMLINE2:   {len(g2_kin)} pairs with kinship estimates")

# ── Load Refined IBD kinship (aggregate across chromosomes) ──────────────
ribd_frames = []
for chrom in [20, 21, 22]:
    path = os.path.join(
        PROJDIR, "results", f"chr{chrom}_refined_IBD", "refinedibd_pairwise_ibd012.csv"
    )
    if os.path.exists(path):
        df = pd.read_csv(path)
        ribd_frames.append(df[["id1", "id2", "kinship_phi_from_P"]])

ribd = pd.concat(ribd_frames, ignore_index=True)
ribd["pair"] = ribd.apply(lambda r: norm_pair(r["id1"], r["id2"]), axis=1)
ribd_agg = ribd.groupby("pair", as_index=False)["kinship_phi_from_P"].mean()
ribd_kin = {}
for _, row in ribd_agg.iterrows():
    ribd_kin[row["pair"]] = float(row["kinship_phi_from_P"])
print(f"Refined IBD: {len(ribd_kin)} pairs with kinship estimates")

# ── Find overlapping pairs for each comparison ───────────────────────────
def get_overlap(dict_a, dict_b):
    common = set(dict_a.keys()) & set(dict_b.keys())
    vals_a = [dict_a[p] for p in common]
    vals_b = [dict_b[p] for p in common]
    return np.array(vals_a), np.array(vals_b), len(common)

plink_g2_a, plink_g2_b, n_pg = get_overlap(plink_kin, g2_kin)
plink_ribd_a, plink_ribd_b, n_pr = get_overlap(plink_kin, ribd_kin)
g2_ribd_a, g2_ribd_b, n_gr = get_overlap(g2_kin, ribd_kin)

print(f"\nOverlap: PLINK ∩ GERMLINE2 = {n_pg},  PLINK ∩ Refined IBD = {n_pr},  GERMLINE2 ∩ Refined IBD = {n_gr}")

# ── Correlations ─────────────────────────────────────────────────────────
def safe_corr(a, b):
    if len(a) < 2:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])

r_pg = safe_corr(plink_g2_a, plink_g2_b)
r_pr = safe_corr(plink_ribd_a, plink_ribd_b)
r_gr = safe_corr(g2_ribd_a, g2_ribd_b)

# ── Three-panel scatter plot ─────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

panels = [
    (axes[0], plink_g2_a, plink_g2_b, "PLINK (KING)", "GERMLINE2", n_pg, r_pg, "#2ca02c"),
    (axes[1], plink_ribd_a, plink_ribd_b, "PLINK (KING)", "Refined IBD", n_pr, r_pr, "#ff7f0e"),
    (axes[2], g2_ribd_a, g2_ribd_b, "GERMLINE2", "Refined IBD", n_gr, r_gr, "#9467bd"),
]

for ax, x, y, xlabel, ylabel, n, r, color in panels:
    if n > 0:
        ax.scatter(x, y, alpha=0.4, s=20, c=color, edgecolors="none")
        # Add diagonal reference line
        lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
                max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, "--", color="gray", alpha=0.5, linewidth=1, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
    
    ax.set_xlabel(f"{xlabel} kinship", fontsize=11)
    ax.set_ylabel(f"{ylabel} kinship", fontsize=11)
    ax.set_title(f"{xlabel} vs {ylabel}\n(n={n}, r={r:.3f})", fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.3)

fig.suptitle("Pairwise Kinship Comparison Across Methods", fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()

outpath = os.path.join(PLOTDIR, "plot4_kinship_scatter.png")
fig.savefig(outpath, dpi=200, bbox_inches="tight")
print(f"\nSaved: {outpath}")
plt.close()
