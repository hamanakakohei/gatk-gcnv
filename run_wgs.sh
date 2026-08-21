#!/usr/bin/env bash
set -euo pipefail


# 1
INPUT_JSON=inputs/01/cnv_germline_cohort_workflow.wdl.minimum.json
#WDL=workflows/cnv_germline_cohort_workflow.edit.wdl
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
UUIDS=(
xxx
)

for UUID in ${UUIDS[@]}; do
    bash scripts/01_qc.sh $UUID
    bash scripts/01_cleanup_after_wdl.sh $UUID
done


# 1-FilterIntervals
INPUT_JSON=inputs/01patch/cnv_germline_cohort_workflow.wdl.minimum.json.01patch.random5000
WDL=workflows/cnv_germline_cohort_workflow.edit.01patch.wdl

qsub \
    -N gcnv_01patch \
    -l s_vmem=64G \
    -o logs/01patch/random5000.o.log \
    -e logs/01patch/random5000.e.log \
    qsubs/cromwell.qsub \
    $WDL \
    $INPUT_JSON


# 1-UMAP
source /usr/local/genome/python3-venv/env.sh
INTERVALS=cromwell-executions/CNVGermlineCohortWorkflow_01patch/xxx/call-FilterIntervals/execution/gcnv-chrAll-contig.list.preprocessed.filtered.interval_list

for BIN in 1; do
for HIGH_THR in 1.00; do
for LOW_THR in 0.0; do
for DIM in 25; do
#for LEIDEN_RES in 0.25 0.50 0.75 1.00 1.25 1.50 1.75 2.00; do
python3 scripts/hdf5_to_umap.py \
    --input inputs/01/h5.list \
    --subset-intervals $INTERVALS \
    --load-pickle 19317.pkl \
    --bin-size $BIN \
    --transform log \
    --var-low $LOW_THR \
    --var-high $HIGH_THR \
    --pca-dim $DIM \
    --out-prefix results/01_umap/out.19317.bin${BIN}.highThr${HIGH_THR}.lowThr${LOW_THR}.pcaDim${DIM}
    #--out-prefix results/01_umap/out.19317.bin${BIN}.highThr${HIGH_THR}.lowThr${LOW_THR}.pcaDim${DIM}.leiden${LEIDEN_RES}
    #--n-neighbors 15 \
    #--umap-min-dist 0.3 \
    #--leiden-resolution $LEIDEN_RES \

    #--dump-pickle 19317.pkl \
    #--new-label inputs/01/new_labels.txt \
    #--new-label inputs/01/new_labels2.txt \
    #--plot-clustering \
    #--dump-pickle 3434.pkl \
    #--plot-interval-variance \
    #--plot-pca-variance \
done
done
done
done
#done


# 02-1
SAMPLE_CLU_MAP=inputs/02/sample_batch.list
HDFS=inputs/01/h5.list
jSON_TEMP_COHORT=inputs/02/cnv_germline_cohort_workflow.wdl.minimum.json.02.template
jSON_TEMP_CASE=inputs/03/cnv_germline_case_workflow.edit.wdl.json.template 

while read CLU; do
    scripts/02_prep_json.py \
        --sample_clu_map $SAMPLE_CLU_MAP \
        --hdf_list $HDFS \
        --cluster $CLU \
        --template_cohort $jSON_TEMP_COHORT \
        --template_case $jSON_TEMP_CASE \
        --out_prefix results/02/cluster${CLU}
        #--out_prefix results/02/cnv_germline_cohort_workflow.wdl.minimum.json.02.leiden${CLU}
        #--training_samples tmp.${CLU}.200 \
done < <(tail -n+1 inputs/02/batch.list)


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
)

for UUID in ${UUIDS[@]}; do
    bash scripts/02_qc.sh $UUID
    bash scripts/02_cleanup_after_wdl.sh $UUID
done


# 02-2-qc-2
scripts/02_qc_count_cnv.py \
    --sample-sheet inputs/02_qc/sample_batch_vcfPath.list.17317 \
    --outdir results/02_qc/ \
    --workers 8


# 03-1
# 02-1で作った残りものサンプルのjsonをさらに02-2の情報を加えて編集する
UUIDS=(
xxx
)

for CLU in {0..21}; do
    UUID=${UUIDS[${CLU}]}
    jSON_CASE_02=results/02/cluster${CLU}.case.json
    MODEL_TARS=results/03/cluster${CLU}_model_tars.list
    
    ls cromwell-executions/CNVGermlineCohortWorkflow_02/${UUID}/call-GermlineCNVCallerCohortMode/shard-*/execution/MyPractice-gcnv-model-shard-*.tar.gz > $MODEL_TARS
    
    scripts/03_prep_json.py \
        --template_case $jSON_CASE_02 \
        --model_tars $MODEL_TARS \
        --ploidy_model_tar cromwell-executions/CNVGermlineCohortWorkflow_02/${UUID}/call-DetermineGermlineContigPloidyCohortMode/execution/MyPractice-contig-ploidy-model.tar.gz \
        --intervals cromwell-executions/CNVGermlineCohortWorkflow_02/${UUID}/call-FilterIntervals/execution/gcnv-chrAll-contig.list.preprocessed.filtered.interval_list \
        --out_prefix results/03/cluster$CLU
done


# 03-2
WDL=workflows/cnv_germline_case_workflow.edit.wdl
INPUT_JSON_LIST=results/03/input_jsons.list
STA=8
END=12

qsub -t ${STA}-${END}:1 -tc 1 \
    -N gcnv_03 \
    -l s_vmem=192G \
    qsubs/cromwell.array.qsub \
    $WDL \
    $INPUT_JSON_LIST


# 03-qc
UUIDS=(
xxx
)

for UUID in ${UUIDS[@]}; do
    bash scripts/03_qc.sh $UUID
    bash scripts/03_cleanup_after_wdl.sh $UUID
done 


# 4
ls cromwell-executions/CNVGermlineCaseWorkflow/*/call-PostprocessGermlineCNVCalls/shard-*/execution/genotyped-segments-*.vcf.gz > $VCF_CASE_LIST
ls cromwell-executions/CNVGermlineCohortWorkflow_02/*/call-PostprocessGermlineCNVCalls/shard-*/execution/genotyped-segments-*.vcf.gz > $VCF_COHORT_LIST

SAMPLE_CNV_LIST=inputs/04/cnv_posicon_details_simple.list
VCF_CASE_LIST=results/04/case_vcf.list
VCF_COHORT_LIST=results/04/cohort_vcf.list
VCF_ALL_LIST=results/04/case_and_cohort_vcf.list

while IFS=$'\t:-' read -r SAMPLE PROBAND CHR STA END CNV; do
    VCF=`grep -w $SAMPLE $VCF_ALL_LIST | cut -f2`
    if [ -z $VCF ]; then continue ;fi
    echo $SAMPLE
    echo $CHR
    echo $STA
    echo $END
    echo $CNV
    echo $VCF
    zcat $VCF \
        | bcftools filter \
            --output-type v \
            -i "CHROM=\"$CHR\" && POS < $END && $STA < INFO/END" \
        | awk '$0!~/^#/'
            #-i "CHROM=\"$CHR\" && POS < $END && $STA < INFO/END && GT!=\"0/0\"" \
done < <(tail -n+2 $SAMPLE_CNV_LIST | grep -v ";") > tmp.res
