#!/usr/bin/env python3
import json
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter samples in CNVGermlineCaseWorkflow JSON"
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-s", "--samples", required=True,
        help="Sample list file (one sample ID per line)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # サンプルリストを読み込む（1行1サンプル、空行・コメント行は無視）
    with open(args.samples) as f:
        target_samples = set(
            line.strip()
            for line in f
            if line.strip() and not line.startswith("#")
        )

    # 入力JSONを読み込む
    with open(args.input) as f:
        data = json.load(f)

    # entity_idとcountsを取得
    collect_counts = data["CNVGermlineCaseWorkflow.CollectCounts"]
    entity_ids = collect_counts["entity_id"]
    counts     = collect_counts["counts"]

    # entity_idとcountsをペアにして、target_samplesに含まれるもののみ抽出
    filtered = [
        (eid, cnt)
        for eid, cnt in zip(entity_ids, counts)
        if eid in target_samples
    ]

    # 結果を分解
    filtered_entity_ids = [x[0] for x in filtered]
    filtered_counts     = [x[1] for x in filtered]

    # 元のJSONを上書き
    data["CNVGermlineCaseWorkflow.CollectCounts"]["entity_id"] = filtered_entity_ids
    data["CNVGermlineCaseWorkflow.CollectCounts"]["counts"]    = filtered_counts

    # 出力JSONに書き込む
    with open(args.output, "w") as f:
        json.dump(data, f, indent=2)

    print("Done.")


if __name__ == "__main__":
    main()
