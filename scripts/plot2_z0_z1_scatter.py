#!/usr/bin/env python3
"""
Plot 2: Two-panel IBD comparison

Left panel:  PLINK Z0 vs Z1 scatter — relationship classification from
             whole-genome method-of-moments IBD estimates.
Right panel: GERMLINE2 vs Refined IBD — total shared cM per detected pair,
             showing how much IBD each segment-based method finds.

Sources:
  PLINK:       KING kinship (king_3chr_1000.kin0) + Z0/Z1 (genome_3chr_1000.genome)
  GERMLINE2:   shared_cm_total from pairwise_ibd.csv
  Refined IBD: total_shared_cm from refinedibd_3chr_pairwise_summary.tsv

Usage:
  python3 scripts/plot2_z0_z1_scatter.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

PROJDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLOTDIR = os.path.join(PROJDIR, "plots")
os.makedirs(PLOTDIR, exist_ok=True)

# ── Load PLINK: KING pairs (kinship > 0.0442), then look up Z0/Z1 ─────────
king = pd.read_csv(
    os.path.join(PROJDIR, "results", "king_3chr_1000.kin0"),
    sep="\t",
)
king_related = king[king["KINSHIP"] > 0.0442].copy()
plink_related_pairs = set()
for _, row in king_related.iterrows():
    a, b = sorted([str(row["IID1"]), str(row["IID2"])])
    plink_related_pairs.add((a, b))

genome = pd.read_csv(
    os.path.join(PROJDIR, "results", "genome_3chr_1000.genome"),
    sep=r"\s+",
)
genome["pair"] = genome.apply(
    lambda r: tuple(sorted([str(r["IID1"]), str(r["IID2"])])), axis=1
)
plink_plot = genome[genome["pair"].isin(plink_related_pairs)][["Z0", "Z1"]].copy()
plink_plot.columns = ["z0", "z1"]
print(f"PLINK:       {len(plink_plot)} pairs (KING kinship > 0.0442)")

# ── Load GERMLINE2: total shared cM per pair ──────────────────────────────
g2_pairs_df = pd.read_csv(os.path.join(PROJDIR, "results", "GERMLINE2", "pairwise_ibd.csv"))
print(f"GERMLINE2:   {len(g2_pairs_df)} pairs with shared segments")

# ── Load Refined IBD: total shared cM per pair ───────────────────────────
ribd_summary = pd.read_csv(
    os.path.join(PROJDIR, "results", "refinedibd_3chr_pairwise_summary.tsv"),
    sep="\t",
)
print(f"Refined IBD: {len(ribd_summary)} pairs with shared segments")

# ══════════════════════════════════════════════════════════════════════════
# Two-panel figure
# ══════════════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# ── Left panel: PLINK Z0 vs Z1 ───────────────────────────────────────────
ax1.scatter(
    plink_plot["z0"], plink_plot["z1"],
    alpha=0.5, s=30, c="#1f77b4", edgecolors="none",
)

# Mark expected positions for reference relationships
ref_points = {
    "Parent-child":  (0.0, 1.0),
    "Full siblings":  (0.25, 0.5),
    "Half-siblings":  (0.5, 0.5),
    "First cousins":  (0.75, 0.25),
    "Unrelated":      (1.0, 0.0),
}
for label, (x, y) in ref_points.items():
    ax1.scatter(x, y, marker="*", s=200, c="red", edgecolors="black",
                linewidths=0.8, zorder=5)
    ax1.annotate(
        label, (x, y),
        textcoords="offset points",
        xytext=(8, 8),
        fontsize=8,
        fontstyle="italic",
        color="darkred",
    )

ax1.set_xlabel("Z0  —  P(IBD = 0)", fontsize=12)
ax1.set_ylabel("Z1  —  P(IBD = 1)", fontsize=12)
ax1.set_title(f"PLINK IBD States ({len(plink_plot)} related pairs)", fontsize=13, fontweight="bold")
ax1.set_xlim(-0.05, 1.05)
ax1.set_ylim(-0.05, 1.05)
ax1.grid(True, alpha=0.3)

# ── Right panel: Total shared cM comparison ───────────────────────────────
bins = np.linspace(0, max(g2_pairs_df["shared_cm_total"].max(),
                          ribd_summary["total_shared_cm"].max()) + 0.5, 40)

ax2.hist(
    g2_pairs_df["shared_cm_total"], bins=bins,
    alpha=0.6, color="#2ca02c", label=f"GERMLINE2 ({len(g2_pairs_df)} pairs)",
    edgecolor="white", linewidth=0.5,
)
ax2.hist(
    ribd_summary["total_shared_cm"], bins=bins,
    alpha=0.6, color="#ff7f0e", label=f"Refined IBD ({len(ribd_summary)} pairs)",
    edgecolor="white", linewidth=0.5,
)

ax2.set_xlabel("Total shared IBD per pair (cM)", fontsize=12)
ax2.set_ylabel("Number of pairs", fontsize=12)
ax2.set_title("Segment-Based Methods: Shared IBD per Pair", fontsize=13, fontweight="bold")
ax2.legend(fontsize=10, framealpha=0.9)
ax2.grid(True, alpha=0.3, axis="y")

# Add summary stats as text
g2_median = g2_pairs_df["shared_cm_total"].median()
ribd_median = ribd_summary["total_shared_cm"].median()
stats_text = (
    f"GERMLINE2:   median = {g2_median:.2f} cM\n"
    f"Refined IBD: median = {ribd_median:.2f} cM"
)
ax2.text(0.97, 0.95, stats_text, transform=ax2.transAxes,
         fontsize=9, va="top", ha="right",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.8))

plt.tight_layout()
outpath = os.path.join(PLOTDIR, "plot2_z0_z1_scatter.png")
fig.savefig(outpath, dpi=200, bbox_inches="tight")
print(f"\nSaved: {outpath}")
plt.close()
