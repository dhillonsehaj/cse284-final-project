#!/usr/bin/env python3
"""
Analyze Refined IBD output (.ibd.gz) for chr22 project.

Input (Refined IBD .ibd format, whitespace-delimited):
  sample1 hap1 sample2 hap2 chr start_bp end_bp LOD length_cM

Example row:
  NA18923 1 NA18507 2 22 28470074 29055755 13.28 0.586

Outputs (written to --outdir):
  - refinedibd_segments_summary.csv
  - refinedibd_pairwise_ibd.csv
  - refinedibd_pairwise_ibd012.csv
  - refinedibd_segment_lengths.png
  - refinedibd_shared_cm_hist.png

What it computes:
  1) Segment-level summary stats (n, min/median/mean/max cM; LOD stats)
  2) Pair-level aggregation (ignoring hap columns):
       n_segments, total_shared_cm, max/mean segment length, mean/max LOD
  3) Approximate IBD-state probabilities per pair using hap columns:
       P(IBD>=1), P(IBD2), P(IBD1), P(IBD0), and kinship_phi = 0.5*P2 + 0.25*P1

IBD-state approximation details:
  - We convert bp to cM with the same assumption used when no genetic map is provided:
        cM = bp / 1,000,000  (i.e., 1 cM ≈ 1 Mb)
  - For each pair, we collect coverage intervals in cM space into two buckets:
        bucket0: segments where hap1==0 OR hap2==0
        bucket1: segments where hap1==1 OR hap2==1
    Then:
        P_ge1 = length(union(bucket0 ∪ bucket1)) / chr_len_cm
        P2    = length(intersection(union(bucket0), union(bucket1))) / chr_len_cm
        P1    = P_ge1 - P2
        P0    = 1 - P_ge1
  - This is a practical, project-level approximation (similar in spirit to your GERMLINE hap-mode approach).

Usage:
  python3 scripts/analyze_refinedibd.py \
    --ibd results/chr20_refinedibd.ibd.gz \
    --outdir results/chr20 \
    --min_lod 1.0 \
    --min_cm 0.5

  python3 scripts/analyze_refinedibd.py \
    --ibd results/chr21_refinedibd.ibd.gz \
    --outdir results/chr21 \
    --min_lod 1.0 \
    --min_cm 0.5

  python3 scripts/analyze_refinedibd.py \
    --ibd results/chr22_refinedibd.ibd.gz \
    --outdir results/chr22 \
    --min_lod 1.0 \
    --min_cm 0.5


Optional:
  --chr_len_cm <float>   Override chromosome length in cM used for probabilities.
                         If omitted, inferred from min/max bp in the IBD file (using bp/1e6).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ----------------------------
# Helpers
# ----------------------------

def ensure_outdir(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

def read_refinedibd_ibd(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["id1", "hap1", "id2", "hap2", "chr", "start_bp", "end_bp", "lod", "cm"],
        dtype={
            "id1": str,
            "hap1": np.int64,
            "id2": str,
            "hap2": np.int64,
            "chr": str,
            "start_bp": np.int64,
            "end_bp": np.int64,
            "lod": float,
            "cm": float,
        },
        engine="python",
    )
    return df

def normalize_pair(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a <= b else (b, a)

def bp_to_cM_approx(bp: int) -> float:
    # consistent with "No genetic map specified: using 1 cM = 1 Mb"
    return float(bp) / 1_000_000.0

def intervals_union(intervals: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
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
# Outputs
# ----------------------------

def segment_level_outputs(df: pd.DataFrame, outdir: Path) -> None:
    summary = {
        "n_segments": int(len(df)),
        "mean_cm": float(df["cm"].mean()) if len(df) else 0.0,
        "median_cm": float(df["cm"].median()) if len(df) else 0.0,
        "min_cm": float(df["cm"].min()) if len(df) else 0.0,
        "max_cm": float(df["cm"].max()) if len(df) else 0.0,
        "mean_lod": float(df["lod"].mean()) if len(df) else 0.0,
        "median_lod": float(df["lod"].median()) if len(df) else 0.0,
        "min_lod": float(df["lod"].min()) if len(df) else 0.0,
        "max_lod": float(df["lod"].max()) if len(df) else 0.0,
    }
    pd.DataFrame([summary]).to_csv(outdir / "refinedibd_segments_summary.csv", index=False)

    if len(df):
        plt.figure()
        plt.hist(df["cm"].values, bins=30)
        plt.xlabel("Segment length (cM)")
        plt.ylabel("Count")
        plt.title("Refined IBD segment length distribution")
        plt.tight_layout()
        plt.savefig(outdir / "refinedibd_segment_lengths.png", dpi=200)
        plt.close()

def pairwise_aggregation(df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    if len(df) == 0:
        empty = pd.DataFrame(columns=[
            "id1", "id2", "n_segments", "total_shared_cm", "max_segment_cm",
            "mean_segment_cm", "mean_lod", "max_lod"
        ])
        empty.to_csv(outdir / "refinedibd_pairwise_ibd.csv", index=False)
        return empty

    pairs = [normalize_pair(a, b) for a, b in zip(df["id1"].tolist(), df["id2"].tolist())]
    df2 = df.copy()
    df2["A"] = [p[0] for p in pairs]
    df2["B"] = [p[1] for p in pairs]

    agg = df2.groupby(["A", "B"]).agg(
        n_segments=("cm", "size"),
        total_shared_cm=("cm", "sum"),
        max_segment_cm=("cm", "max"),
        mean_segment_cm=("cm", "mean"),
        mean_lod=("lod", "mean"),
        max_lod=("lod", "max"),
    ).reset_index().rename(columns={"A": "id1", "B": "id2"})

    agg.to_csv(outdir / "refinedibd_pairwise_ibd.csv", index=False)

    plt.figure()
    plt.hist(agg["total_shared_cm"].values, bins=30)
    plt.xlabel("Total shared length per pair (sum of cM)")
    plt.ylabel("Count")
    plt.title("Per-pair total shared cM (Refined IBD)")
    plt.tight_layout()
    plt.savefig(outdir / "refinedibd_shared_cm_hist.png", dpi=200)
    plt.close()

    return agg

def pairwise_ibd012(df: pd.DataFrame, chr_len_cm: float, outdir: Path) -> pd.DataFrame:
    """
    Approximate P(IBD0/1/2) per pair using hap columns + interval coverage in cM space.
    """
    if len(df) == 0:
        empty = pd.DataFrame(columns=[
            "id1", "id2",
            "covered_cm_union", "covered_cm_ibd2_overlap",
            "P_IBD0", "P_IBD1", "P_IBD2", "P_IBD_ge1",
            "kinship_phi_from_P",
            "covered_cm_bucket0", "covered_cm_bucket1",
        ])
        empty.to_csv(outdir / "refinedibd_pairwise_ibd012.csv", index=False)
        return empty

    # normalize pairs
    pairs = [normalize_pair(a, b) for a, b in zip(df["id1"].tolist(), df["id2"].tolist())]
    df2 = df.copy()
    df2["A"] = [p[0] for p in pairs]
    df2["B"] = [p[1] for p in pairs]

    # Convert bp endpoints -> cM endpoints (approx)
    cM0 = df2["start_bp"].apply(lambda x: bp_to_cM_approx(int(x)))
    cM1 = df2["end_bp"].apply(lambda x: bp_to_cM_approx(int(x)))
    df2["cM_start"] = np.minimum(cM0, cM1)
    df2["cM_end"] = np.maximum(cM0, cM1)

    # Build bucket intervals per pair
    pair_to_int0: Dict[Tuple[str, str], List[Tuple[float, float]]] = {}
    pair_to_int1: Dict[Tuple[str, str], List[Tuple[float, float]]] = {}

    for _, r in df2.iterrows():
        key = (r["A"], r["B"])
        interval = (float(r["cM_start"]), float(r["cM_end"]))

        # Put interval in bucket0 if either hap index == 0
        if int(r["hap1"]) == 0 or int(r["hap2"]) == 0:
            pair_to_int0.setdefault(key, []).append(interval)

        # Put interval in bucket1 if either hap index == 1
        if int(r["hap1"]) == 1 or int(r["hap2"]) == 1:
            pair_to_int1.setdefault(key, []).append(interval)

    rows = []
    keys = sorted(set(pair_to_int0.keys()) | set(pair_to_int1.keys()))
    for (a, b) in keys:
        u0 = intervals_union(pair_to_int0.get((a, b), []))
        u1 = intervals_union(pair_to_int1.get((a, b), []))
        u_all = intervals_union(u0 + u1)
        inter = intervals_intersection(u0, u1)

        L0 = intervals_length(u0)
        L1 = intervals_length(u1)
        Lall = intervals_length(u_all)
        L2 = intervals_length(inter)

        # Probabilities (clipped)
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
            "P_IBD_ge1": P_ge1,
            "kinship_phi_from_P": kinship_phi,
            "covered_cm_bucket0": L0,
            "covered_cm_bucket1": L1,
        })

    out = pd.DataFrame(rows)
    out.to_csv(outdir / "refinedibd_pairwise_ibd012.csv", index=False)
    return out


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ibd", required=True, help="Path to Refined IBD .ibd.gz output")
    ap.add_argument("--outdir", default="results", help="Output directory for CSVs/plots")
    ap.add_argument("--min_lod", type=float, default=None, help="Optional filter: keep segments with LOD >= this")
    ap.add_argument("--min_cm", type=float, default=None, help="Optional filter: keep segments with cM >= this")
    ap.add_argument("--chr_len_cm", type=float, default=None,
                    help="Override chromosome length (cM) used for IBD probabilities. "
                         "If omitted, inferred from min/max bp in the IBD file using bp/1e6.")
    args = ap.parse_args()

    ibd_path = Path(args.ibd)
    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    if not ibd_path.exists():
        raise FileNotFoundError(f"Could not find Refined IBD output at: {ibd_path}")

    df = read_refinedibd_ibd(ibd_path)
    print(f"[info] Loaded {len(df)} segments from {ibd_path}")

    # Optional filters
    if args.min_lod is not None:
        before = len(df)
        df = df[df["lod"] >= args.min_lod].copy()
        print(f"[info] Filtered by min_lod={args.min_lod}: {before} -> {len(df)}")

    if args.min_cm is not None:
        before = len(df)
        df = df[df["cm"] >= args.min_cm].copy()
        print(f"[info] Filtered by min_cm={args.min_cm}: {before} -> {len(df)}")

    # Determine chr length in cM for probability normalization
    chr_len_cm = args.chr_len_cm
    if chr_len_cm is None:
        if len(df) == 0:
            chr_len_cm = 0.0
        else:
            min_bp = int(min(df["start_bp"].min(), df["end_bp"].min()))
            max_bp = int(max(df["start_bp"].max(), df["end_bp"].max()))
            chr_len_cm = bp_to_cM_approx(max_bp) - bp_to_cM_approx(min_bp)

    print(f"[info] chr length (cM) used for probabilities: {chr_len_cm:.4f}")

    segment_level_outputs(df, outdir)
    pair = pairwise_aggregation(df, outdir)

    out012 = pairwise_ibd012(df, chr_len_cm, outdir) if chr_len_cm > 0 else pd.DataFrame()
    print(f"[info] Wrote refinedibd_pairwise_ibd.csv with {len(pair)} pairs")
    if chr_len_cm > 0:
        print(f"[info] Wrote refinedibd_pairwise_ibd012.csv with {len(out012)} pairs")

    # Small console summary
    if len(df):
        print("[summary] segments:", len(df))
        print("[summary] cM min/median/max:",
              float(df["cm"].min()), float(df["cm"].median()), float(df["cm"].max()))
        print("[summary] unique pairs (any segment):", len(pair))
        if chr_len_cm > 0 and len(out012):
            top = out012.sort_values("P_IBD_ge1", ascending=False).head(5)
            print("[summary] top pairs by P(IBD>=1):")
            print(top[["id1", "id2", "P_IBD_ge1", "P_IBD2", "P_IBD1", "kinship_phi_from_P"]].to_string(index=False))

    print(f"[done] Outputs written to: {outdir.resolve()}")


if __name__ == "__main__":
    main()