#!/usr/bin/env bash
set -euo pipefail


# 1
#KIT=twist_v2
KIT=twist_v1
INPUT_JSON=inputs/01/cnv_germline_cohort_workflow.wdl.minimum.wes.json.${KIT}.01
WDL=workflows/cnv_germline_cohort_workflow.edit.01.wdl

qsub \
    -N gcnv_01 \
    -l s_vmem=312G \
    -o logs/01/${KIT}.o.log \
    -e logs/01/${KIT}.e.log \
    qsubs/cromwell.qsub \
    $WDL \
    $INPUT_JSON


# 1-qc
uuids=(
xxx
xxx
)

for uuid in ${uuids[@]}; do
    bash scripts/01_qc.sh $uuid
    bash scripts/01_cleanup_after_wdl.sh $uuid 
done


# 02-1
SAMPLE_CLU_MAP=inputs/02/sample_batch.list
HDFS=inputs/01/h5.list
jSON_TEMP_CASE=inputs/03/cnv_germline_case_workflow.edit.wdl.json.template 

# twist v1
jSON_TEMP_COHORT=inputs/02/cnv_germline_cohort_workflow.wdl.minimum.json.02.wes.template_twist_v1
while read CLU; do
    scripts/02_prep_json.py \
        --sample_clu_map $SAMPLE_CLU_MAP \
        --hdf_list $HDFS \
        --cluster $CLU \
        --template_cohort $jSON_TEMP_COHORT \
        --template_case $jSON_TEMP_CASE \
        --out_prefix results/02/cluster${CLU}
done < <(tail -n+160 inputs/02/batch.list | head -n25)

# twist v2
jSON_TEMP_COHORT=inputs/02/cnv_germline_cohort_workflow.wdl.minimum.json.02.wes.template_twist_v2
while read CLU; do
    scripts/02_prep_json.py \
        --sample_clu_map $SAMPLE_CLU_MAP \
        --hdf_list $HDFS \
        --cluster $CLU \
        --template_cohort $jSON_TEMP_COHORT \
        --template_case $jSON_TEMP_CASE \
        --out_prefix results/02/cluster${CLU}
done < <(tail -n+185 inputs/02/batch.list)


# 02-2
INPUT_JSON_LIST=results/02/input_jsons.list
STA=1
END=196

qsub -t ${STA}-${END}:1 -tc 1 \
    -N gcnv_02 \
    -l s_vmem=64G \
    qsubs/cromwell.array.qsub \
    $WDL \
    $INPUT_JSON_LIST


# 02-2-qc
UUIDS=(
xxx
xxx
)

for UUID in ${UUIDS[@]}; do
    bash scripts/02_qc.sh $UUID
    bash scripts/02_cleanup_after_wdl.sh $UUID
done


# 02-2-qc-2
scripts/02_qc_count_cnv.py \
    --sample-sheet inputs/02_qc/sample_batch_vcfPath.list \
    --outdir results/02_qc/ \
    --workers 8
