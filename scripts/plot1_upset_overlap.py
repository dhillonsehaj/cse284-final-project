#!/usr/bin/env python3
"""
Plot 1: IBD Pair Detection Overlap — UpSet Plot

Compares which sample pairs are detected as related by each of the three
IBD methods (PLINK, GERMLINE2, Refined IBD). A pair is considered "detected"
if it appears in the method's output:

  - PLINK KING:  kinship > 0.0442 (standard 3rd-degree threshold)
  - GERMLINE2:   any pair in pairwise_ibd.csv (has ≥1 shared segment)
  - Refined IBD: any pair in refinedibd_3chr_pairwise_summary.tsv

Usage:
  python3 scripts/plot1_upset_overlap.py
"""

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships

PROJDIR = __file__.rsplit("/scripts/", 1)[0]

# ── Load PLINK KING pairs (kinship > 0.0442) ──────────────────────────────
king = pd.read_csv(
    f"{PROJDIR}/results/king_3chr_1000.kin0",
    sep="\t",
)
king_related = king[king["KINSHIP"] > 0.0442].copy()
plink_pairs = set()
for _, row in king_related.iterrows():
    a, b = sorted([str(row["IID1"]), str(row["IID2"])])
    plink_pairs.add((a, b))

print(f"PLINK:       {len(plink_pairs)} pairs (KING kinship > 0.0442)")

# ── Load GERMLINE2 pairs ──────────────────────────────────────────────────
g2 = pd.read_csv(f"{PROJDIR}/results/GERMLINE2/pairwise_ibd.csv")
g2_pairs = set()
for _, row in g2.iterrows():
    a, b = sorted([str(row["id1"]), str(row["id2"])])
    g2_pairs.add((a, b))

print(f"GERMLINE2:   {len(g2_pairs)} pairs (any shared segment)")

# ── Load Refined IBD pairs ────────────────────────────────────────────────
ribd = pd.read_csv(
    f"{PROJDIR}/results/refinedibd_3chr_pairwise_summary.tsv",
    sep="\t",
)
ribd_pairs = set()
for _, row in ribd.iterrows():
    a, b = sorted([str(row["ID1"]), str(row["ID2"])])
    ribd_pairs.add((a, b))

print(f"Refined IBD: {len(ribd_pairs)} pairs (any shared segment)")

# ── Build membership list for UpSet ───────────────────────────────────────
all_pairs = plink_pairs | g2_pairs | ribd_pairs
print(f"Union:       {len(all_pairs)} unique pairs total")

memberships = []
for pair in all_pairs:
    cats = []
    if pair in plink_pairs:
        cats.append("PLINK")
    if pair in g2_pairs:
        cats.append("GERMLINE2")
    if pair in ribd_pairs:
        cats.append("Refined IBD")
    memberships.append(cats)

data = from_memberships(memberships)

# ── Print overlap statistics ──────────────────────────────────────────────
only_plink = plink_pairs - g2_pairs - ribd_pairs
only_g2 = g2_pairs - plink_pairs - ribd_pairs
only_ribd = ribd_pairs - plink_pairs - g2_pairs
all_three = plink_pairs & g2_pairs & ribd_pairs
plink_g2 = (plink_pairs & g2_pairs) - ribd_pairs
plink_ribd = (plink_pairs & ribd_pairs) - g2_pairs
g2_ribd = (g2_pairs & ribd_pairs) - plink_pairs

print(f"\nOverlap breakdown:")
print(f"  PLINK only:              {len(only_plink)}")
print(f"  GERMLINE2 only:          {len(only_g2)}")
print(f"  Refined IBD only:        {len(only_ribd)}")
print(f"  PLINK ∩ GERMLINE2 only:  {len(plink_g2)}")
print(f"  PLINK ∩ Refined IBD only:{len(plink_ribd)}")
print(f"  GERMLINE2 ∩ Refined IBD: {len(g2_ribd)}")
print(f"  All three:               {len(all_three)}")

# ── Generate UpSet plot ───────────────────────────────────────────────────
fig = plt.figure(figsize=(10, 6))
upset = UpSet(
    data,
    subset_size="count",
    show_counts=True,
    sort_by="cardinality",
    sort_categories_by="cardinality",
    facecolor="steelblue",
)
upset.plot(fig=fig)
fig.suptitle(
    "IBD Pair Detection Overlap Across Methods",
    fontsize=14,
    fontweight="bold",
    y=1.02,
)

outpath = f"{PROJDIR}/plots/plot1_upset_overlap.png"
fig.savefig(outpath, dpi=200, bbox_inches="tight")
print(f"\nSaved: {outpath}")
plt.close()
