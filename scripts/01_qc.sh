#!/bin/bash
set -euo pipefail

RUN_ID=$1
PREFIX="../../cromwell-executions/CNVGermlineCohortWorkflow_01/${RUN_ID}/"

outdir="results/01_qc/"
mkdir -p "$outdir"
cd "$outdir"


# 2. hdf5のをリスト作成
ls ${PREFIX}call-CollectCounts/shard-*/execution/*.hdf5 > hdf5.${RUN_ID}.list


# 3. 全サンプル x 全タスクの成功したかの一覧を得る
collect_rc() {
    local exec_glob="$1"
    local input_glob="${exec_glob/execution\/rc/inputs}"

    shopt -s nullglob

    # rc ファイルを配列で取得
    local rc_files=( $exec_glob )
    local input_dirs=( $input_glob )

    shopt -u nullglob

    local n_expected=${#input_dirs[@]}
    local n_found=${#rc_files[@]}

    if (( n_found > 0 )); then
        local vals=()
        for rc in "${rc_files[@]}"; do
            vals+=( "$(tr -d '\n' < "$rc")" )
        done
        printf "%s/%d" "$(IFS=';'; echo "${vals[*]}")" "$n_expected"
    else
        printf "NA/%d" "$n_expected"
    fi
}

echo -e "PreprocessIntervals\t\
AnnotateIntervals\t\
CollectCounts" > qc_simple.${RUN_ID}.out

collect_rc "${PREFIX}call-PreprocessIntervals/execution/rc"   >> qc_simple.${RUN_ID}.out; printf "\t" >> qc_simple.${RUN_ID}.out
collect_rc "${PREFIX}call-AnnotateIntervals/execution/rc"     >> qc_simple.${RUN_ID}.out; printf "\t" >> qc_simple.${RUN_ID}.out
collect_rc "${PREFIX}call-CollectCounts/shard-*/execution/rc" >> qc_simple.${RUN_ID}.out; printf "\n" >> qc_simple.${RUN_ID}.out
