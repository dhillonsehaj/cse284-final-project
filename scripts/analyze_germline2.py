#!/usr/bin/env python3
"""
Analyze GERMLINE2 (g2) output for chr22 project.

Expected repo layout (based on your screenshot):
  data_raw/
    chr22_g2_out          (GERMLINE2 output, tab-delimited segments)
    chr22_g2_out_hap      (optional: output from `./g2 -h ...`)
    chr22_g2.map          (pos  cm/Mb  cM)
    chr22_g2.sample       (IMPUTE/SHAPEIT sample file; 2 header lines)
  results/                (will be created if missing)

Usage:
  python3 analyze_germline2.py \
    --prefix data_raw/chr22_g2_out \
    --map data_raw/chr22_g2.map \
    --sample data_raw/chr22_g2.sample \
    --outdir results

If you have haplotype mode output, also pass:
  --hap_prefix data_raw/chr22_g2_out_hap
"""

from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------
# Utilities
# ----------------------------

def ensure_outdir(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

def read_g2_segments(path: Path) -> pd.DataFrame:
    """
    Read GERMLINE2 output file (e.g., data_raw/chr22_g2_out) which is a TSV with:
      ID1 ID2 P0 P1 cM #words #gaps
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["id1", "id2", "p0", "p1", "cm", "words", "gaps"],
        dtype={"id1": str, "id2": str, "p0": np.int64, "p1": np.int64, "cm": float, "words": np.int64, "gaps": np.int64},
    )
    # Some versions may output spaces; handle if someone opened/saved weirdly
    if df.shape[1] != 7:
        df = pd.read_csv(path, sep=r"\s+", header=None,
                         names=["id1","id2","p0","p1","cm","words","gaps"])
    return df

def read_map(map_path: Path) -> pd.DataFrame:
    """
    chr22_g2.map has 3 columns: pos  cm/Mb  cM
    We'll keep pos and cM.
    """
    m = pd.read_csv(map_path, sep=r"\s+|\t", engine="python", header=None, names=["pos", "cM_per_Mb", "cM"])
    m = m.sort_values("pos").reset_index(drop=True)
    if (m["pos"].diff().fillna(1) < 0).any():
        m = m.sort_values("pos").reset_index(drop=True)
    return m[["pos", "cM"]]

def read_samples(sample_path: Path) -> List[str]:
    """
    IMPUTE2/SHAPEIT sample file:
      line1: ID_1 ID_2 missing
      line2: 0 0 0
      line3+: sample rows, use column 2 as sample ID
    """
    ids = []
    with open(sample_path, "r", encoding="utf-8", errors="ignore") as f:
        for i, line in enumerate(f, start=1):
            if i <= 2:
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                ids.append(parts[1])
    return ids

def strip_hap_suffix(sample_id: str) -> Tuple[str, Optional[int]]:
    """
    If haploid mode is on, g2 appends '.0' or '.1' to IDs.
    Return (base_id, hap_index or None).
    """
    m = re.match(r"^(.*)\.(0|1)$", sample_id)
    if not m:
        return sample_id, None
    return m.group(1), int(m.group(2))

def bp_to_cM(pos_bp: int, map_df: pd.DataFrame) -> float:
    """
    Convert bp position to cM by linear interpolation on the provided map.
    map_df: columns ['pos','cM'] sorted by pos.
    """
    pos_arr = map_df["pos"].to_numpy()
    cm_arr = map_df["cM"].to_numpy()

    if pos_bp <= pos_arr[0]:
        return float(cm_arr[0])
    if pos_bp >= pos_arr[-1]:
        return float(cm_arr[-1])

    idx = np.searchsorted(pos_arr, pos_bp, side="right")
    left = idx - 1
    right = idx

    x0, x1 = pos_arr[left], pos_arr[right]
    y0, y1 = cm_arr[left], cm_arr[right]
    if x1 == x0:
        return float(y0)
    t = (pos_bp - x0) / (x1 - x0)
    return float(y0 + t * (y1 - y0))

def intervals_union(intervals: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """
    Union of half-open intervals [start, end] in cM space.
    """
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ms, me = merged[-1]
        if s <= me:  # overlap/touch
            merged[-1] = (ms, max(me, e))
        else:
            merged.append((s, e))
    return merged

def intervals_length(intervals: List[Tuple[float, float]]) -> float:
    return float(sum(max(0.0, e - s) for s, e in intervals))

def intervals_intersection(a: List[Tuple[float, float]], b: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """
    Intersection of two union-merged interval lists.
    """
    out: List[Tuple[float, float]] = []
    i = j = 0
    while i < len(a) and j < len(b):
        s1, e1 = a[i]
        s2, e2 = b[j]
        s = max(s1, s2)
        e = min(e1, e2)
        if s < e:
            out.append((s, e))
        if e1 < e2:
            i += 1
        else:
            j += 1
    return out


# ----------------------------
# Main analysis
# ----------------------------

def segment_level_outputs(df: pd.DataFrame, outdir: Path) -> None:
    # Basic segment stats
    summary = {
        "n_segments": int(len(df)),
        "mean_cm": float(df["cm"].mean()) if len(df) else 0.0,
        "median_cm": float(df["cm"].median()) if len(df) else 0.0,
        "min_cm": float(df["cm"].min()) if len(df) else 0.0,
        "max_cm": float(df["cm"].max()) if len(df) else 0.0,
    }
    pd.DataFrame([summary]).to_csv(outdir / "segments_summary.csv", index=False)

    # Plot segment length distribution
    if len(df):
        plt.figure()
        plt.hist(df["cm"].values, bins=50)
        plt.xlabel("Segment length (cM)")
        plt.ylabel("Count")
        plt.title("GERMLINE2 segment length distribution")
        plt.tight_layout()
        plt.savefig(outdir / "segment_lengths.png", dpi=200)
        plt.close()

def pairwise_ibd_ge1(df: pd.DataFrame, chr_len_cm: float, outdir: Path) -> pd.DataFrame:
    """
    Compute per-pair:
      shared_cm_total = sum of segment cM
      P_IBD_ge1 ≈ shared_cm_total / chr_len_cm  (clipped to [0,1])
    NOTE: This can exceed chr_len due to overlapping segments being double-counted.
          For a more precise fraction, use haploid mode computation below which unions
          intervals in cM space.
    """
    pair = df.groupby(["id1","id2"], as_index=False)["cm"].sum().rename(columns={"cm":"shared_cm_total"})
    pair["P_IBD_ge1_approx"] = (pair["shared_cm_total"] / chr_len_cm).clip(0.0, 1.0)

    pair.to_csv(outdir / "pairwise_ibd.csv", index=False)

    # plot shared cm distribution
    if len(pair):
        plt.figure()
        plt.hist(pair["shared_cm_total"].values, bins=50)
        plt.xlabel("Total shared length per pair (sum of cM)")
        plt.ylabel("Count")
        plt.title("Per-pair total shared cM (GERMLINE2)")
        plt.tight_layout()
        plt.savefig(outdir / "shared_cm_hist.png", dpi=200)
        plt.close()

    return pair

def pairwise_ibd012_from_hap(
    hap_df: pd.DataFrame,
    map_df: pd.DataFrame,
    chr_len_cm: float,
    outdir: Path,
) -> pd.DataFrame:
    """
    Compute P(IBD=0,1,2) per pair using haploid-mode segments.

    Approach:
      - hap_df contains IDs like NA12878.0 NA12891.1 ...
      - For each base pair (A,B):
          collect intervals for hap0 side and hap1 side separately (regardless of which hap matched which)
          union them on chr22 in cM space
          P2 = length(intersection(union(hap0_intervals), union(hap1_intervals))) / chr_len_cm
          Pge1 = length(union(hap0_intervals ∪ hap1_intervals)) / chr_len_cm
          P1 = Pge1 - P2
          P0 = 1 - Pge1
    Notes:
      - This is a practical approximation that treats overlapping coverage on both haplotypes as IBD2.
      - Works well for project-level reporting.
    """
    # Parse hap suffix
    base1, hap1 = zip(*[strip_hap_suffix(x) for x in hap_df["id1"].tolist()])
    base2, hap2 = zip(*[strip_hap_suffix(x) for x in hap_df["id2"].tolist()])
    hap_df = hap_df.copy()
    hap_df["base1"] = list(base1)
    hap_df["base2"] = list(base2)
    hap_df["hap1"] = list(hap1)
    hap_df["hap2"] = list(hap2)

    # Convert bp endpoints to cM endpoints for interval union/intersection
    # Do it vectorized-ish with apply (fine for project scale)
    hap_df["cM0"] = hap_df["p0"].apply(lambda x: bp_to_cM(int(x), map_df))
    hap_df["cM1"] = hap_df["p1"].apply(lambda x: bp_to_cM(int(x), map_df))
    # Ensure start<=end
    hap_df["cM_start"] = hap_df[["cM0","cM1"]].min(axis=1)
    hap_df["cM_end"] = hap_df[["cM0","cM1"]].max(axis=1)

    # Group by unordered pair (A,B)
    def normalize_pair(a: str, b: str) -> Tuple[str, str]:
        return (a, b) if a <= b else (b, a)

    hap_df["A"], hap_df["B"] = zip(*[normalize_pair(a,b) for a,b in zip(hap_df["base1"], hap_df["base2"])])

    # For each record, decide whether it contributes to "hap0 coverage" or "hap1 coverage"
    # We don't care which individual's hap it is; we care whether there are 0/1/2 haplotypes worth of coverage.
    # A segment is between two haplotypes; classify by whether it involves hap0 or hap1 on each side.
    # We'll add its interval to bucket 0 if either endpoint hap==0; similarly for bucket 1.
    # This heuristic is common for project approximations.
    def buckets(row) -> List[int]:
        bs = []
        if row["hap1"] == 0 or row["hap2"] == 0:
            bs.append(0)
        if row["hap1"] == 1 or row["hap2"] == 1:
            bs.append(1)
        return bs or [0]  # fallback

    # Build interval lists
    pair_to_int0: Dict[Tuple[str,str], List[Tuple[float,float]]] = {}
    pair_to_int1: Dict[Tuple[str,str], List[Tuple[float,float]]] = {}

    for _, r in hap_df.iterrows():
        key = (r["A"], r["B"])
        interval = (float(r["cM_start"]), float(r["cM_end"]))
        for b in buckets(r):
            if b == 0:
                pair_to_int0.setdefault(key, []).append(interval)
            else:
                pair_to_int1.setdefault(key, []).append(interval)

    rows = []
    for key in sorted(set(pair_to_int0.keys()) | set(pair_to_int1.keys())):
        a, b = key
        u0 = intervals_union(pair_to_int0.get(key, []))
        u1 = intervals_union(pair_to_int1.get(key, []))
        u_all = intervals_union(u0 + u1)
        inter = intervals_intersection(u0, u1)

        L0 = intervals_length(u0)
        L1 = intervals_length(u1)
        Lall = intervals_length(u_all)
        L2 = intervals_length(inter)

        P_ge1 = min(1.0, max(0.0, Lall / chr_len_cm))
        P2 = min(P_ge1, max(0.0, L2 / chr_len_cm))
        P1 = max(0.0, P_ge1 - P2)
        P0 = max(0.0, 1.0 - P_ge1)

        kinship_phi = 0.5 * P2 + 0.25 * P1

        rows.append({
            "id1": a,
            "id2": b,
            "covered_cm_union": Lall,
            "covered_cm_ibd2_overlap": L2,
            "P_IBD0": P0,
            "P_IBD1": P1,
            "P_IBD2": P2,
            "kinship_phi_from_P": kinship_phi,
            "covered_cm_bucket0": L0,
            "covered_cm_bucket1": L1,
        })

    out = pd.DataFrame(rows)
    out.to_csv(outdir / "pairwise_ibd012.csv", index=False)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, help="Path to GERMLINE2 output prefix file (e.g., data_raw/chr22_g2_out).")
    ap.add_argument("--hap_prefix", default=None, help="Optional haploid-mode output prefix (e.g., data_raw/chr22_g2_out_hap).")
    ap.add_argument("--map", required=True, help="Genetic map file (pos cm/Mb cM), e.g., data_raw/chr22_g2.map")
    ap.add_argument("--sample", required=True, help="Sample file, e.g., data_raw/chr22_g2.sample")
    ap.add_argument("--outdir", default="results", help="Output directory for CSVs/plots")
    ap.add_argument("--chr_len_cm", default=None, type=float,
                    help="Override chromosome length in cM. If omitted, computed from map as max(cM)-min(cM).")
    args = ap.parse_args()

    prefix = Path(args.prefix)
    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    map_df = read_map(Path(args.map))
    samples = read_samples(Path(args.sample))
    print(f"[info] Loaded {len(samples)} samples from {args.sample}")

    chr_len_cm = args.chr_len_cm
    if chr_len_cm is None:
        chr_len_cm = float(map_df["cM"].iloc[-1] - map_df["cM"].iloc[0])
    print(f"[info] chr length (cM) used: {chr_len_cm:.4f}")

    # Read main segments
    if not prefix.exists():
        raise FileNotFoundError(f"Could not find g2 output file at: {prefix}")

    seg = read_g2_segments(prefix)
    print(f"[info] Loaded {len(seg)} segments from {prefix}")

    # Segment-level summary + plots
    segment_level_outputs(seg, outdir)

    # Pairwise approx P(IBD>=1) by summing cM (fast, may double-count overlaps)
    pair = pairwise_ibd_ge1(seg, chr_len_cm, outdir)
    print(f"[info] Wrote pairwise_ibd.csv with {len(pair)} pairs")

    # If hap file present, compute P(IBD0/1/2) by union/intersection on cM intervals
    if args.hap_prefix:
        hap_path = Path(args.hap_prefix)
        if hap_path.exists():
            hap_seg = read_g2_segments(hap_path)
            print(f"[info] Loaded {len(hap_seg)} hap segments from {hap_path}")
            out012 = pairwise_ibd012_from_hap(hap_seg, map_df, chr_len_cm, outdir)
            print(f"[info] Wrote pairwise_ibd012.csv with {len(out012)} pairs")
        else:
            print(f"[warn] hap_prefix provided but file not found: {hap_path}")

    print(f"[done] Outputs written to: {outdir.resolve()}")


if __name__ == "__main__":
    main()