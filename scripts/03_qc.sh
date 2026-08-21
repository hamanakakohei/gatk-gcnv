#!/bin/bash
set -euo pipefail

RUN_ID=$1
PREFIX="../../cromwell-executions/CNVGermlineCaseWorkflow/${RUN_ID}/"

outdir="results/03_qc/"
mkdir -p "$outdir"
cd "$outdir"


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

echo -e "call-DetermineGermlineContigPloidyCaseMode\t\
call-GermlineCNVCallerCaseMode\t\
call-PostprocessGermlineCNVCalls\t\
call-ScatterIntervals\t\
call-ScatterPloidyCallsBySample" > qc_simple.${RUN_ID}.out

collect_rc "${PREFIX}call-DetermineGermlineContigPloidyCaseMode/execution/rc"   >> qc_simple.${RUN_ID}.out; printf "\t" >> qc_simple.${RUN_ID}.out
collect_rc "${PREFIX}call-GermlineCNVCallerCaseMode/shard-*/execution/rc"       >> qc_simple.${RUN_ID}.out; printf "\t" >> qc_simple.${RUN_ID}.out
collect_rc "${PREFIX}call-PostprocessGermlineCNVCalls/shard-*/execution/rc"     >> qc_simple.${RUN_ID}.out; printf "\t" >> qc_simple.${RUN_ID}.out
collect_rc "${PREFIX}call-ScatterIntervals/execution/rc"                        >> qc_simple.${RUN_ID}.out; printf "\t" >> qc_simple.${RUN_ID}.out
collect_rc "${PREFIX}call-ScatterPloidyCallsBySample/execution/rc"              >> qc_simple.${RUN_ID}.out; printf "\t" >> qc_simple.${RUN_ID}.out
           
