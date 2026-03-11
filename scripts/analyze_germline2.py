#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ----------------------------
# Utilities
# ----------------------------

def ensure_outdir(outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)


def infer_chr_label(path: Path) -> str:
    m = re.search(r"chr(\d+)", path.name)
    if m:
        return f"chr{m.group(1)}"
    return path.stem


def read_g2_segments(path: Path, chr_label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Could not find g2 output file: {path}")

    if path.stat().st_size == 0:
        return pd.DataFrame(columns=["id1", "id2", "p0", "p1", "cm", "words", "gaps", "chr"])

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["id1", "id2", "p0", "p1", "cm", "words", "gaps"],
            dtype={
                "id1": str,
                "id2": str,
                "p0": np.int64,
                "p1": np.int64,
                "cm": float,
                "words": np.int64,
                "gaps": np.int64,
            },
        )
        if df.shape[1] != 7:
            raise ValueError("Unexpected number of columns with tab separator")
    except Exception:
        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["id1", "id2", "p0", "p1", "cm", "words", "gaps"],
        )

    df["chr"] = chr_label
    return df


def read_map(map_path: Path) -> pd.DataFrame:
    m = pd.read_csv(
        map_path,
        sep=r"\s+|\t",
        engine="python",
        header=None,
        names=["pos", "cM_per_Mb", "cM"],
    )
    m = m.sort_values("pos").reset_index(drop=True)
    return m[["pos", "cM"]]


def read_samples(sample_path: Path) -> List[str]:
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
    m = re.match(r"^(.*)\.(0|1)$", sample_id)
    if not m:
        return sample_id, None
    return m.group(1), int(m.group(2))


def bp_to_cM(pos_bp: int, map_df: pd.DataFrame) -> float:
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
    if not intervals:
        return []

    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]

    for s, e in intervals[1:]:
        ms, me = merged[-1]
        if s <= me:
            merged[-1] = (ms, max(me, e))
        else:
            merged.append((s, e))
    return merged


def intervals_length(intervals: List[Tuple[float, float]]) -> float:
    return float(sum(max(0.0, e - s) for s, e in intervals))


def intervals_intersection(a: List[Tuple[float, float]], b: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    out = []
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


def parse_elapsed_to_seconds(s: str) -> float:
    """
    Convert /usr/bin/time -v elapsed time string to seconds.
    Examples:
      0:12.41 -> 12.41
      1:02.03 -> 62.03
      1:02:03 -> 3723
    """
    s = s.strip()
    parts = s.split(":")
    if len(parts) == 2:
        minutes = int(parts[0])
        seconds = float(parts[1])
        return minutes * 60 + seconds
    elif len(parts) == 3:
        hours = int(parts[0])
        minutes = int(parts[1])
        seconds = float(parts[2])
        return hours * 3600 + minutes * 60 + seconds
    else:
        return float("nan")


def parse_time_log(log_path: Path, mode: str) -> Dict[str, object]:
    """
    Parse a /usr/bin/time -v log file.
    """
    if not log_path.exists():
        raise FileNotFoundError(f"Runtime log not found: {log_path}")

    chr_label = infer_chr_label(log_path)

    record = {
        "chromosome": chr_label,
        "mode": mode,
        "elapsed_time_raw": None,
        "elapsed_seconds": np.nan,
        "user_time_s": np.nan,
        "system_time_s": np.nan,
        "max_memory_kb": np.nan,
        "max_memory_mb": np.nan,
        "exit_status": np.nan,
    }

    with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if line.startswith("Elapsed (wall clock) time"):
                value = line.split(":", 1)[1].strip()
                record["elapsed_time_raw"] = value
                record["elapsed_seconds"] = parse_elapsed_to_seconds(value)

            elif line.startswith("User time (seconds)"):
                record["user_time_s"] = float(line.split(":", 1)[1].strip())

            elif line.startswith("System time (seconds)"):
                record["system_time_s"] = float(line.split(":", 1)[1].strip())

            elif line.startswith("Maximum resident set size (kbytes)"):
                kb = float(line.split(":", 1)[1].strip())
                record["max_memory_kb"] = kb
                record["max_memory_mb"] = kb / 1024.0

            elif line.startswith("Exit status"):
                record["exit_status"] = int(line.split(":", 1)[1].strip())

    return record


def runtime_outputs(
    diploid_logs: List[Path],
    haploid_logs: List[Path],
    outdir: Path,
) -> None:
    """
    Parse runtime/memory logs and write summary tables.
    """
    rows = []

    for path in diploid_logs:
        rows.append(parse_time_log(path, mode="diploid"))

    for path in haploid_logs:
        rows.append(parse_time_log(path, mode="haploid"))

    if not rows:
        return

    df = pd.DataFrame(rows)
    df = df.sort_values(["mode", "chromosome"]).reset_index(drop=True)
    df.to_csv(outdir / "runtime_summary.csv", index=False)

    by_mode = (
        df.groupby("mode", as_index=False)
        .agg(
            n_runs=("mode", "count"),
            total_elapsed_seconds=("elapsed_seconds", "sum"),
            mean_elapsed_seconds=("elapsed_seconds", "mean"),
            max_elapsed_seconds=("elapsed_seconds", "max"),
            mean_user_time_s=("user_time_s", "mean"),
            mean_system_time_s=("system_time_s", "mean"),
            max_memory_mb=("max_memory_mb", "max"),
            mean_memory_mb=("max_memory_mb", "mean"),
        )
    )
    by_mode.to_csv(outdir / "runtime_summary_by_mode.csv", index=False)


def filter_self_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove self-pairs, including haploid IDs like HG00096.0 vs HG00096.1.
    """
    if len(df) == 0:
        return df

    def base_id(x: str) -> str:
        return re.sub(r"\.(0|1)$", "", str(x))

    keep = df.apply(lambda row: base_id(row["id1"]) != base_id(row["id2"]), axis=1)
    return df.loc[keep].reset_index(drop=True)


# ----------------------------
# Plot helpers
# ----------------------------

def save_hist(data, xlabel: str, ylabel: str, title: str, outfile: Path, bins: int = 50) -> None:
    if len(data) == 0:
        return
    plt.figure()
    plt.hist(data, bins=bins)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()


# ----------------------------
# Segment summaries
# ----------------------------

def segment_level_outputs(df: pd.DataFrame, outdir: Path) -> None:
    summary = {
        "n_segments": int(len(df)),
        "mean_cm": float(df["cm"].mean()) if len(df) else 0.0,
        "median_cm": float(df["cm"].median()) if len(df) else 0.0,
        "min_cm": float(df["cm"].min()) if len(df) else 0.0,
        "max_cm": float(df["cm"].max()) if len(df) else 0.0,
    }
    pd.DataFrame([summary]).to_csv(outdir / "segments_summary.csv", index=False)

    if len(df):
        chr_summary = (
            df.groupby("chr")["cm"]
            .agg(
                n_segments="count",
                mean_cm="mean",
                median_cm="median",
                min_cm="min",
                max_cm="max",
            )
            .reset_index()
        )
        chr_summary.to_csv(outdir / "segments_summary_by_chr.csv", index=False)

        save_hist(
            df["cm"].values,
            "Segment length (cM)",
            "Count",
            "GERMLINE2 segment length distribution",
            outdir / "segment_lengths.png",
        )

        plt.figure()
        for chr_label, sub in df.groupby("chr"):
            plt.hist(sub["cm"].values, bins=50, alpha=0.5, label=chr_label)
        plt.xlabel("Segment length (cM)")
        plt.ylabel("Count")
        plt.title("GERMLINE2 segment length distribution by chromosome")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "segment_lengths_by_chr.png", dpi=200)
        plt.close()


# ----------------------------
# Pairwise summaries
# ----------------------------

def pairwise_ibd_ge1(df: pd.DataFrame, total_len_cm: float, outdir: Path) -> pd.DataFrame:
    if len(df) == 0:
        empty = pd.DataFrame(columns=[
            "id1", "id2", "shared_cm_total", "n_segments",
            "max_segment_cm", "mean_segment_cm", "P_IBD_ge1_approx"
        ])
        empty.to_csv(outdir / "pairwise_ibd.csv", index=False)
        return empty

    pair = (
        df.groupby(["id1", "id2"], as_index=False)
        .agg(
            shared_cm_total=("cm", "sum"),
            n_segments=("cm", "count"),
            max_segment_cm=("cm", "max"),
            mean_segment_cm=("cm", "mean"),
        )
    )
    pair["P_IBD_ge1_approx"] = (pair["shared_cm_total"] / total_len_cm).clip(0.0, 1.0)
    pair.to_csv(outdir / "pairwise_ibd.csv", index=False)

    save_hist(
        pair["shared_cm_total"].values,
        "Total shared length per pair (cM)",
        "Count",
        "Per-pair total shared cM (GERMLINE2)",
        outdir / "shared_cm_hist.png",
    )

    save_hist(
        pair["n_segments"].values,
        "Number of segments per pair",
        "Count",
        "Per-pair segment counts (GERMLINE2)",
        outdir / "segment_count_hist.png",
    )

    save_hist(
        pair["max_segment_cm"].values,
        "Maximum segment length per pair (cM)",
        "Count",
        "Per-pair maximum segment length (GERMLINE2)",
        outdir / "max_segment_hist.png",
    )

    return pair


def pairwise_ibd_by_chr(df: pd.DataFrame, chr_len_by_chr: Dict[str, float], outdir: Path) -> pd.DataFrame:
    if len(df) == 0:
        empty = pd.DataFrame(columns=[
            "chr", "id1", "id2", "shared_cm_total", "n_segments",
            "max_segment_cm", "mean_segment_cm", "P_IBD_ge1_approx_chr"
        ])
        empty.to_csv(outdir / "pairwise_ibd_by_chr.csv", index=False)
        return empty

    rows = []
    for chr_label, sub in df.groupby("chr"):
        chr_len = chr_len_by_chr[chr_label]
        pair_chr = (
            sub.groupby(["id1", "id2"], as_index=False)
            .agg(
                shared_cm_total=("cm", "sum"),
                n_segments=("cm", "count"),
                max_segment_cm=("cm", "max"),
                mean_segment_cm=("cm", "mean"),
            )
        )
        pair_chr["P_IBD_ge1_approx_chr"] = (pair_chr["shared_cm_total"] / chr_len).clip(0.0, 1.0)
        pair_chr["chr"] = chr_label
        rows.append(pair_chr)

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(outdir / "pairwise_ibd_by_chr.csv", index=False)
    return out


# ----------------------------
# Hap-based P(IBD=0,1,2)
# ----------------------------

def pairwise_ibd012_from_hap(
    hap_df: pd.DataFrame,
    map_by_chr: Dict[str, pd.DataFrame],
    total_len_cm: float,
    outdir: Path,
) -> pd.DataFrame:
    if len(hap_df) == 0:
        out = pd.DataFrame(columns=[
            "id1", "id2", "covered_cm_union", "covered_cm_ibd2_overlap",
            "P_IBD0", "P_IBD1", "P_IBD2", "kinship_phi_from_P",
            "covered_cm_bucket0", "covered_cm_bucket1"
        ])
        out.to_csv(outdir / "pairwise_ibd012.csv", index=False)
        return out

    base1, hap1 = zip(*[strip_hap_suffix(x) for x in hap_df["id1"].tolist()])
    base2, hap2 = zip(*[strip_hap_suffix(x) for x in hap_df["id2"].tolist()])

    hap_df = hap_df.copy()
    hap_df["base1"] = list(base1)
    hap_df["base2"] = list(base2)
    hap_df["hap1"] = list(hap1)
    hap_df["hap2"] = list(hap2)

    cM_start_list = []
    cM_end_list = []

    for _, row in hap_df.iterrows():
        map_df = map_by_chr[row["chr"]]
        cm0 = bp_to_cM(int(row["p0"]), map_df)
        cm1 = bp_to_cM(int(row["p1"]), map_df)
        cM_start_list.append(min(cm0, cm1))
        cM_end_list.append(max(cm0, cm1))

    hap_df["cM_start"] = cM_start_list
    hap_df["cM_end"] = cM_end_list

    def normalize_pair(a: str, b: str) -> Tuple[str, str]:
        return (a, b) if a <= b else (b, a)

    hap_df["A"], hap_df["B"] = zip(*[normalize_pair(a, b) for a, b in zip(hap_df["base1"], hap_df["base2"])])

    def buckets(row) -> List[int]:
        bs = []
        if row["hap1"] == 0 or row["hap2"] == 0:
            bs.append(0)
        if row["hap1"] == 1 or row["hap2"] == 1:
            bs.append(1)
        return bs or [0]

    pair_to_int0: Dict[Tuple[str, str], List[Tuple[float, float]]] = {}
    pair_to_int1: Dict[Tuple[str, str], List[Tuple[float, float]]] = {}

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

        P_ge1 = min(1.0, max(0.0, Lall / total_len_cm))
        P2 = min(P_ge1, max(0.0, L2 / total_len_cm))
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


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--prefix",
        action="append",
        required=True,
        help="GERMLINE2 diploid output file. Use once per chromosome.",
    )
    ap.add_argument(
        "--hap_prefix",
        action="append",
        default=[],
        help="GERMLINE2 haploid output file. Use once per chromosome.",
    )
    ap.add_argument(
        "--map",
        action="append",
        required=True,
        help="Chromosome-specific map file. Use once per chromosome.",
    )
    ap.add_argument("--sample", required=True, help="Sample file")
    ap.add_argument("--outdir", default="results", help="Output directory")
    ap.add_argument(
        "--runtime_log",
        action="append",
        default=[],
        help="Diploid /usr/bin/time -v log file. Use once per chromosome.",
    )
    ap.add_argument(
        "--hap_runtime_log",
        action="append",
        default=[],
        help="Haploid /usr/bin/time -v log file. Use once per chromosome.",
    )

    args = ap.parse_args()

    outdir = Path(args.outdir)
    ensure_outdir(outdir)

    sample_path = Path(args.sample)
    samples = read_samples(sample_path)
    print(f"[info] Loaded {len(samples)} samples from {sample_path}")

    map_by_chr: Dict[str, pd.DataFrame] = {}
    chr_len_by_chr: Dict[str, float] = {}

    for map_str in args.map:
        map_path = Path(map_str)
        chr_label = infer_chr_label(map_path)
        map_df = read_map(map_path)
        map_by_chr[chr_label] = map_df
        chr_len_by_chr[chr_label] = float(map_df["cM"].iloc[-1] - map_df["cM"].iloc[0])

    total_len_cm = sum(chr_len_by_chr.values())
    print(f"[info] Chromosomes loaded: {sorted(map_by_chr.keys())}")
    print(f"[info] Total genetic length used: {total_len_cm:.4f} cM")

    seg_frames = []
    for prefix_str in args.prefix:
        prefix_path = Path(prefix_str)
        chr_label = infer_chr_label(prefix_path)
        seg_df = read_g2_segments(prefix_path, chr_label)
        seg_frames.append(seg_df)
        print(f"[info] Loaded {len(seg_df)} diploid segments from {prefix_path}")

    seg = pd.concat(seg_frames, ignore_index=True) if seg_frames else pd.DataFrame(
        columns=["id1", "id2", "p0", "p1", "cm", "words", "gaps", "chr"]
    )

    n_before = len(seg)
    seg = filter_self_pairs(seg)
    print(f"[info] Filtered self-pairs from diploid segments: {n_before - len(seg)} removed")

    segment_level_outputs(seg, outdir)

    pair = pairwise_ibd_ge1(seg, total_len_cm, outdir)
    print(f"[info] Wrote pairwise_ibd.csv with {len(pair)} pairs")

    pair_chr = pairwise_ibd_by_chr(seg, chr_len_by_chr, outdir)
    print(f"[info] Wrote pairwise_ibd_by_chr.csv with {len(pair_chr)} rows")

    if args.hap_prefix:
        hap_frames = []
        for hap_str in args.hap_prefix:
            hap_path = Path(hap_str)
            chr_label = infer_chr_label(hap_path)
            hap_df = read_g2_segments(hap_path, chr_label)
            hap_frames.append(hap_df)
            print(f"[info] Loaded {len(hap_df)} haploid segments from {hap_path}")

        hap_seg = pd.concat(hap_frames, ignore_index=True) if hap_frames else pd.DataFrame(
            columns=["id1", "id2", "p0", "p1", "cm", "words", "gaps", "chr"]
        )

        n_before_hap = len(hap_seg)
        hap_seg = filter_self_pairs(hap_seg)
        print(f"[info] Filtered self-pairs from haploid segments: {n_before_hap - len(hap_seg)} removed")

        out012 = pairwise_ibd012_from_hap(hap_seg, map_by_chr, total_len_cm, outdir)
        print(f"[info] Wrote pairwise_ibd012.csv with {len(out012)} pairs")

    diploid_logs = [Path(x) for x in args.runtime_log]
    haploid_logs = [Path(x) for x in args.hap_runtime_log]

    if diploid_logs or haploid_logs:
        runtime_outputs(diploid_logs, haploid_logs, outdir)
        print("[info] Wrote runtime_summary.csv and runtime_summary_by_mode.csv")

    print(f"[done] Outputs written to {outdir.resolve()}")


if __name__ == "__main__":
    main()