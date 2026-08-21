#!/usr/bin/env python3
"""
VCF (PostprocessGermlineCNVCalls出力) をもとにサンプルごとにCNV数をバッチ毎にボックスプロットする

前提：
- ALTが<DEL> or <DUP>
- 10列目がGT:CN:etc.
"""
import argparse
import gzip
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

METRICS = [
    "DEL_or_DUP", "DEL", "DUP",
    "DEL_CN0", "DEL_CN1",
    "DUP_CN3", "DUP_CN4", "DUP_CN5",
]

MEDIAN_METRICS = [m + "_ChrMedian" for m in METRICS]
ALL_METRICS = METRICS + MEDIAN_METRICS

# 常染色体 (chr1~chr22 / 1~22。"chr"有無どちらでも受け付ける)
AUTOSOMES = [str(i) for i in range(1, 23)]


def normalize_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def count_sample(vcf_path: str) -> dict:
    # 染色体ごとの内訳をとっておき、常染色体だけで 合計 と 染色体内中央値 を出す
    per_chr = {c: {m: 0 for m in METRICS} for c in AUTOSOMES}
    #counts = {m: 0 for m in METRICS}

    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom = normalize_chrom(fields[0])
            if chrom not in per_chr:
                continue  # 性染色体・その他コンティグは除外

            alt = fields[4]
            if alt not in ("<DEL>", "<DUP>"):
                continue

            # FORMAT = GT:CN:NP:etc. 前提なので、サンプル列の2番目(index=1)がCN
            sample_vals = fields[9].split(":")
            try:
                cn = int(sample_vals[1])
            except (ValueError, IndexError):
                print(f"[WARN] {vcf_path}: CN取得失敗の行をスキップ -> {line.strip()[:80]}", file=sys.stderr)
                continue

            c = per_chr[chrom]
            c["DEL_or_DUP"] += 1
            if alt == "<DEL>":
                c["DEL"] += 1
                if cn == 0:
                    c["DEL_CN0"] += 1
                elif cn == 1:
                    c["DEL_CN1"] += 1
            else:  # <DUP>
                c["DUP"] += 1
                if cn == 3:
                    c["DUP_CN3"] += 1
                elif cn == 4:
                    c["DUP_CN4"] += 1
                elif cn == 5:
                    c["DUP_CN5"] += 1
 
    counts = {}
    for m in METRICS:
        per_chr_values = [per_chr[c][m] for c in AUTOSOMES]
        counts[m] = sum(per_chr_values)
        counts[m + "_ChrMedian"] = float(np.median(per_chr_values))
    return counts


def natural_sort_key(s: str):
    """'110_111' や '150-2' のようなバッチ名を、数字部分は数値として、
    それ以外は文字列として比較する自然順ソート用のキーを返す。"""
    return [int(tok) if tok.isdigit() else tok.lower() for tok in re.split(r"(\d+)", str(s))]


def _process_one(args):
    sample_id, vcf_path = args
    counts = count_sample(vcf_path)
    return sample_id, counts


def main():
    ap = argparse.ArgumentParser(description="")
    ap.add_argument("--sample-sheet", required=True,
                     help="Sample_ID<TAB>Batch<TAB>VCF_path の3列、へっだーなし")
    ap.add_argument("--outdir", default="results/02/")
    ap.add_argument("--workers", type=int, default=8, help="VCF解析の並列数")
    args = ap.parse_args()


    # サンプルシート読む
    sheet_df = pd.read_csv(args.sample_sheet, sep="\t", header=None, names=["Sample_ID", "Batch", "VCF_path"], dtype=str)
    sheet_df["Batch"] = sheet_df["Batch"].astype(str).str.strip()


    # 並列に各サンプルを処理する
    results = {}
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {
            ex.submit(_process_one, (row.Sample_ID, row.VCF_path)): row.Sample_ID
            for row in sheet_df.itertuples()
        }
        for i, fut in enumerate(as_completed(futures), 1):
            sid, counts = fut.result()
            results[sid] = counts
            if i % 100 == 0 or i == len(sheet_df):
                print(f"[INFO] 処理済み {i}/{len(sheet_df)}", file=sys.stderr)

    rows = []
    missing = []
    for sid in sheet_df["Sample_ID"]:
        counts = results.get(sid)
        if counts is None:
            missing.append(sid)
            continue
        row = {"Sample_ID": sid}
        row.update(counts)
        rows.append(row)

    if missing:
        print(f"[WARN] VCF未取得/解析失敗サンプル: {len(missing)}件 -> {missing[:20]}{' ...' if len(missing) > 20 else ''}",
              file=sys.stderr)


    # dfに統合する？
    count_df = pd.DataFrame(rows, columns=["Sample_ID"] + ALL_METRICS)
    count_df = count_df.merge(sheet_df[["Sample_ID", "Batch"]], on="Sample_ID", how="left")

    # サンプル別集計
    per_sample_tsv = os.path.join(args.outdir, "per_sample_counts.tsv")
    count_df.to_csv(per_sample_tsv, sep="\t", index=False)

    # バッチの並び順は自然順(数値部分は数値として、それ以外は文字列として比較)
    all_batches = sorted(count_df["Batch"].unique(), key=natural_sort_key)

    # バッチ別要約統計
    batch_summary = count_df.groupby("Batch")[ALL_METRICS].describe().reindex(all_batches)
    batch_summary_tsv = os.path.join(args.outdir, "batch_summary_stats.tsv")
    batch_summary.to_csv(batch_summary_tsv, sep="\t")

    # 指標ごとにボックスプロット
    tick_values = [0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
    for metric in ALL_METRICS:
        data = [
            np.log10(count_df.loc[count_df["Batch"] == b, metric].dropna().values + 1)
            for b in all_batches
        ]
        fig_w = max(12, len(all_batches) * 0.15)
        fig, ax = plt.subplots(figsize=(fig_w, 6))
        ax.boxplot(data, positions=range(1, len(all_batches) + 1), showfliers=True,
            patch_artist=True,
            boxprops=dict(facecolor="blue", edgecolor="black"),
            flierprops=dict(markerfacecolor="white", markeredgecolor="black")
        )
        ax.set_xticks(range(1, len(all_batches) + 1))
        ax.set_xticklabels(all_batches, rotation=90, fontsize=6)
        ax.set_xlabel("Sequence batch")
        ax.set_yticks(np.log10(np.array(tick_values) + 1))
        ax.set_yticklabels([str(v) for v in tick_values])
        ax.set_ylabel("Number of segments")
        ax.set_title(metric)
        fig.tight_layout()
        out_png = os.path.join(args.outdir, f"boxplot_{metric}.png")
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        print(f"[INFO] done.")


if __name__ == "__main__":
    main()
