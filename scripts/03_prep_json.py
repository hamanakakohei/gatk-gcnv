#!/usr/bin/env python
import argparse
import json
import pandas as pd
from pathlib import Path
import glob
import random


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--template_case", required=True, help="gatk-gcnv case modeのinput jsonのtemplate")
    p.add_argument("--model_tars", required=True)
    p.add_argument("--ploidy_model_tar", required=True)
    p.add_argument("--intervals", required=True)
    p.add_argument("--out_prefix", required=True)
    return p.parse_args()


def load_data(model_tars):
    tar = pd.read_csv(model_tars, sep="\t", header=None)
    tar.columns = ["path"]
    return list(tar["path"])


def load_template(path):
    with open(path) as f:
        return json.load(f)


def update_case_json(template, ploidy_model_tar, intervals, tar_paths):
    j = template.copy()
    j["CNVGermlineCaseWorkflow.contig_ploidy_model_tar"] = ploidy_model_tar
    j["CNVGermlineCaseWorkflow.filtered_intervals"] = intervals
    j["CNVGermlineCaseWorkflow.gcnv_model_tars"] = tar_paths
    return j


def main():
    args = parse_args()

   
    # --- model tar ---
    tar_paths = load_data(args.model_tars)


    # --- template ---
    tmpl_case = load_template(args.template_case)


    # --- case ---
    json_oth = update_case_json(
        tmpl_case,
        ploidy_model_tar = args.ploidy_model_tar,
        intervals = args.intervals,
        tar_paths = tar_paths
    )

    #
    with open(f"{args.out_prefix}.json", "w") as f:
        json.dump(json_oth, f, indent=2)


if __name__ == "__main__":
    main()
