#!/usr/bin/env python3
"""
Plot 3 – IBD Segment Length Distribution
Compares GERMLINE2 and Refined IBD segment lengths (in cM).
PLINK does not produce individual segment data so it is excluded.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip, glob, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ── Load GERMLINE2 segments ─────────────────────────────────────────────────
g2_files = sorted(glob.glob(os.path.join(ROOT, "results/GERMLINE2/chr*_1000_g2_out")))
g2_rows = []
for f in g2_files:
    chr_label = os.path.basename(f).split("_")[0]  # e.g. "chr20"
    with open(f) as fh:
        for line in fh:
            cols = line.strip().split()
            # cols: id1, id2, start_bp, end_bp, length_cM, n_markers, ...
            g2_rows.append({
                "chr": chr_label,
                "id1": cols[0],
                "id2": cols[1],
                "start_bp": int(cols[2]),
                "end_bp": int(cols[3]),
                "length_cM": float(cols[4]),
            })
g2 = pd.DataFrame(g2_rows)
print(f"GERMLINE2 segments loaded: {len(g2)}")

# ── Load Refined IBD segments ───────────────────────────────────────────────
ribd_path = os.path.join(ROOT, "results/refinedibd_3chr_all.ibd.gz")
ribd_rows = []
with gzip.open(ribd_path, "rt") as fh:
    for line in fh:
        cols = line.strip().split()
        # cols: id1, hap1, id2, hap2, chr, start_bp, end_bp, LOD, length_cM
        ribd_rows.append({
            "chr": f"chr{cols[4]}",
            "id1": cols[0],
            "id2": cols[2],
            "start_bp": int(cols[5]),
            "end_bp": int(cols[6]),
            "length_cM": float(cols[8]),
        })
ribd = pd.DataFrame(ribd_rows)
print(f"Refined IBD segments loaded: {len(ribd)}")

# ── Also compute segment physical length (Mb) ──────────────────────────────
g2["length_Mb"] = (g2["end_bp"] - g2["start_bp"]) / 1e6
ribd["length_Mb"] = (ribd["end_bp"] - ribd["start_bp"]) / 1e6

# ── Summary stats ───────────────────────────────────────────────────────────
for label, df in [("GERMLINE2", g2), ("Refined IBD", ribd)]:
    print(f"\n{label}:")
    print(f"  Segments: {len(df)}")
    print(f"  cM  – mean={df['length_cM'].mean():.2f}, median={df['length_cM'].median():.2f}, "
          f"min={df['length_cM'].min():.2f}, max={df['length_cM'].max():.2f}")
    print(f"  Mb  – mean={df['length_Mb'].mean():.2f}, median={df['length_Mb'].median():.2f}, "
          f"min={df['length_Mb'].min():.2f}, max={df['length_Mb'].max():.2f}")

# ── Figure: 2-panel (cM distribution + Mb distribution) ────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# --- Panel A: cM ---
bins_cm = np.linspace(0, max(g2["length_cM"].max(), ribd["length_cM"].max()) * 1.05, 50)
axes[0].hist(g2["length_cM"], bins=bins_cm, alpha=0.55, label=f"GERMLINE2 (n={len(g2)})",
             color="#1f77b4", edgecolor="white", linewidth=0.4)
axes[0].hist(ribd["length_cM"], bins=bins_cm, alpha=0.55, label=f"Refined IBD (n={len(ribd)})",
             color="#ff7f0e", edgecolor="white", linewidth=0.4)
axes[0].axvline(g2["length_cM"].median(), color="#1f77b4", ls="--", lw=1.2,
                label=f"G2 median={g2['length_cM'].median():.2f} cM")
axes[0].axvline(ribd["length_cM"].median(), color="#ff7f0e", ls="--", lw=1.2,
                label=f"RIBD median={ribd['length_cM'].median():.2f} cM")
axes[0].set_xlabel("Segment Length (cM)")
axes[0].set_ylabel("Number of Segments")
axes[0].set_title("A) Segment Length Distribution (cM)")
axes[0].legend(fontsize=8)

# --- Panel B: Mb ---
bins_mb = np.linspace(0, max(g2["length_Mb"].max(), ribd["length_Mb"].max()) * 1.05, 50)
axes[1].hist(g2["length_Mb"], bins=bins_mb, alpha=0.55, label=f"GERMLINE2 (n={len(g2)})",
             color="#1f77b4", edgecolor="white", linewidth=0.4)
axes[1].hist(ribd["length_Mb"], bins=bins_mb, alpha=0.55, label=f"Refined IBD (n={len(ribd)})",
             color="#ff7f0e", edgecolor="white", linewidth=0.4)
axes[1].axvline(g2["length_Mb"].median(), color="#1f77b4", ls="--", lw=1.2,
                label=f"G2 median={g2['length_Mb'].median():.2f} Mb")
axes[1].axvline(ribd["length_Mb"].median(), color="#ff7f0e", ls="--", lw=1.2,
                label=f"RIBD median={ribd['length_Mb'].median():.2f} Mb")
axes[1].set_xlabel("Segment Length (Mb)")
axes[1].set_ylabel("Number of Segments")
axes[1].set_title("B) Segment Length Distribution (Mb)")
axes[1].legend(fontsize=8)

fig.suptitle("IBD Segment Length Distribution – GERMLINE2 vs Refined IBD",
             fontsize=13, fontweight="bold", y=1.02)
plt.tight_layout()

out = os.path.join(ROOT, "plots", "plot3_segment_length.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=200, bbox_inches="tight")
print(f"\nSaved: {out}")
