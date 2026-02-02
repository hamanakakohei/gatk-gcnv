#!/usr/bin/env python3
import argparse
import h5py
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.decomposition import PCA
import umap
import matplotlib.pyplot as plt


def load_counts_hdf5(path):
    with h5py.File(path, "r") as f:
        counts = f["/counts/values"][0]  # (N,)
        intervals = f["/intervals/transposed_index_start_end"][:]  # (3, N)
        contig_names = f["/intervals/indexed_contig_names"][:].astype(str)

    contig_idx = intervals[0].astype(int)
    start = intervals[1].astype(int)
    end = intervals[2].astype(int)

    autosome_mask = contig_idx < 22  # chr1–22

    df = pd.DataFrame({
        "chr": [contig_names[i] for i in contig_idx[autosome_mask]],
        "start": start[autosome_mask],
        "end": end[autosome_mask],
        "read_count": counts[autosome_mask]
    })

    return df


def load_multiple_samples(paths):
    matrices = []
    sample_names = []

    for p in paths:
        df = load_counts_hdf5(p)
        matrices.append(df["read_count"].values)
        sample_names.append(Path(p).stem)

    X = np.vstack(matrices)
    return X, sample_names


def bin_intervals(X, bin_size):
    if bin_size <= 1:
        return X
    n_bins = X.shape[1] // bin_size
    X = X[:, :n_bins * bin_size]
    return X.reshape(X.shape[0], n_bins, bin_size).sum(axis=2)


def transform(X, method):
    if method == "log":
        return np.log1p(X)
    if method == "sqrt":
        return np.sqrt(X)
    return X


def variance_filter(X, low, high):
    var = X.var(axis=0)
    lo = np.quantile(var, low)
    hi = np.quantile(var, high)
    mask = (var > lo) & (var < hi)
    return X[:, mask]


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", required=True, type=str, help="hdf5ファイルのリストが書かれたファイル")
    parser.add_argument("--bin-size", type=int, default=1, help="Merge N consecutive intervals")
    parser.add_argument("--transform", choices=["log", "sqrt"], default="log")
    parser.add_argument("--var-low", type=float, default=0.01, help="Lower quantile for variance filter")
    parser.add_argument("--var-high", type=float, default=0.99, help="Upper quantile for variance filter")
    parser.add_argument("--pca-dim", type=int, default=20, help="Number of PC for UMAP")
    parser.add_argument("--umap-neighbors", type=int, default=15)
    parser.add_argument("--umap-min-dist", type=float, default=0.3)
    parser.add_argument("--sample-class", help="サンプルメタデータファイル: sample<TAB>class")
    parser.add_argument("--out-prefix", default="out")
    args = parser.parse_args()

    # 準備（bin集約 --> トータルで割る --> log/sqrt変換 --> 分散でフィルター）
    with open(args.input) as f:
        hdf5_paths = [l.strip() for l in f]

    X, sample_names = load_multiple_samples(hdf5_paths)
    X = bin_intervals(X, args.bin_size)
    X = X / X.sum(axis=1, keepdims=True)
    X = transform(X, args.transform)
    X = variance_filter(X, args.var_low, args.var_high)

    if args.sample_class:
        batch_df = pd.read_csv(
            args.sample_class, sep="\t", header=None,
            names=["sample", "batch"]
        )
        batch_map = dict(zip(batch_df.sample, batch_df.batch))
        batch_labels = [batch_map.get(s, "NA") for s in sample_names]
    else:
        batch_labels = ["NA"] * len(sample_names)

    # PCA --> UMAP
    pca = PCA(n_components=args.pca_dim)
    X_pca = pca.fit_transform(X)

    reducer = umap.UMAP(
        n_neighbors=args.umap_neighbors,
        min_dist=args.umap_min_dist,
        random_state=0
    )
    embedding = reducer.fit_transform(X_pca)

    # Plot
    plt.figure(figsize=(6, 6))
    scatter = plt.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=pd.factorize(batch_labels)[0],
        cmap="tab10",
        s=30
    )
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("Coverage UMAP")
    plt.tight_layout()
    plt.savefig(f"{args.out_prefix}.umap.png", dpi=300)


if __name__ == "__main__":
    main()
