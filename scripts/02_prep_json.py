#!/usr/bin/env python
#
# 注：なんかたまにサンプル数が元のファイルからcase,cohortファイルになった時に減ってしまう事がある気がする、、、バグ？
import argparse
import json
import pandas as pd
from pathlib import Path
import glob
import random


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample_clu_map", required=True, default="results/01_umap/out.3434.bin1.highThr1.00.lowThr0.9.pcaDim20.leiden0.25.umap.tsv")
    p.add_argument("--hdf_list", required=True, default="inputs/01/h5.3434.list")
    p.add_argument("--cluster", type=str, required=True)
    
    p.add_argument("--n_sample", type=int, default=200)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--training_samples", type=None, help="デフォルトではランダムに200サンプル選ばれるが、ファイルでサンプルを指定できる、サンプル名は.cram.countsつけなくて良い")

    p.add_argument("--template_cohort", required=True, help="gatk-gcnv cohort modeのinput jsonのtemplate")
    p.add_argument("--template_case", required=True, help="gatk-gcnv case modeのinput jsonのtemplate")

    p.add_argument("--out_prefix", required=True)

    return p.parse_args()


def load_data(sample_clu_map, hdf_list):
    df = pd.read_csv(sample_clu_map, sep="\t", skiprows=1) # header=None)
    df.columns = ["col1", "col2", "sample", "cluster"]
    df["sample"] = df["sample"].apply(lambda x: x.split(".")[0])

    hdf = pd.read_csv(hdf_list, sep="\t", header=None, usecols=[0])
    hdf.columns = ["path"]
    hdf["sample"] = hdf["path"].apply(lambda x: Path(x).name.split(".")[0])

    return df, hdf


def pick_samples(df, clu, n, seed, training_samples=None):
    samples = df[df["cluster"] == clu]["sample"]

    if training_samples is not None:
        samples_200 = pd.read_csv(training_samples, header=None)[0]
        samples_200 = samples_200[samples_200.isin(samples)]
    else:
        samples_200 = samples.sample(n=min(n, len(samples)), random_state=seed)

    samples_oth = samples[~samples.isin(samples_200)]

    return samples, samples_200, samples_oth


def map_hdfs_ordered(hdf_df, samples):
    return (
        pd.DataFrame({"sample": samples})
        .merge(hdf_df, on="sample", how="left")
        ["path"]
        .tolist()
    )


def load_template(path):
    with open(path) as f:
        return json.load(f)


def update_json(template, hdfs, samples, wf_type="CNVGermlineCohortWorkflow_02"):
    j = template.copy()
    j[f"{wf_type}.CollectCounts"]["counts"] = hdfs
    j[f"{wf_type}.CollectCounts"]["entity_id"] = list(samples)

    return j


def main():
    args = parse_args()

    df, hdf = load_data(args.sample_clu_map, args.hdf_list)

    # --- sample & hdf ---
    samples, samples_200, samples_oth = pick_samples(
        df, args.cluster, args.n_sample, args.seed, args.training_samples
    )

    hdfs_all = map_hdfs_ordered(hdf, samples)
    hdfs_200 = map_hdfs_ordered(hdf, samples_200)
    hdfs_oth = map_hdfs_ordered(hdf, samples_oth)


    # --- template ---
    tmpl_cohort = load_template(args.template_cohort)
    tmpl_case = load_template(args.template_case)


    # --- cohort ---
    json_200 = update_json(tmpl_cohort, hdfs_200, samples_200, wf_type="CNVGermlineCohortWorkflow_02")

    with open(f"{args.out_prefix}.cohort.json", "w") as f:
        json.dump(json_200, f, indent=2)


    # --- case ---
    json_oth = update_json(tmpl_case, hdfs_oth, samples_oth, wf_type="CNVGermlineCaseWorkflow")

    with open(f"{args.out_prefix}.case.json", "w") as f:
        json.dump(json_oth, f, indent=2)


if __name__ == "__main__":
    main()
