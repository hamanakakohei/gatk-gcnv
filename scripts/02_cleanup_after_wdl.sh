#!/usr/bin/env bash
set -euo pipefail
uuid=$1


ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/execution/*.preprocessed.annotated.tsv
ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/inputs/*/*.preprocessed.interval_list
ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/inputs/*/human_g1k_v37_fix.fasta
ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/inputs/*/k100.umap.20260509.hg19.bed.gz
#ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-DetermineGermlineContigPloidyCohortMode/execution
ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-DetermineGermlineContigPloidyCohortMode/inputs/*/*.hdf5
ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-DetermineGermlineContigPloidyCohortMode/inputs/*/*.preprocessed.filtered.interval_list
ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/inputs/*/*.hdf5
ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/inputs/*/*.preprocessed.annotated.tsv
#ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/execution/*/*gcnv-calls*.tar.gz
#ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/execution/*/out/*-calls/SAMPLE_*/
ls cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-PostprocessGermlineCNVCalls/shard-*/inputs/*/*

#chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/execution/*.preprocessed.annotated.tsv
#chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/inputs/*/*.preprocessed.interval_list
#chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/inputs/*/human_g1k_v37_fix.fasta
#chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/inputs/*/k100.umap.20260509.hg19.bed.gz
##chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-DetermineGermlineContigPloidyCohortMode/execution
#chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-DetermineGermlineContigPloidyCohortMode/inputs/*/*.hdf5
#chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-DetermineGermlineContigPloidyCohortMode/inputs/*/*.preprocessed.filtered.interval_list
#chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/inputs/*/*.hdf5
#chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/inputs/*/*.preprocessed.annotated.tsv
##chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/execution/*/*gcnv-calls*.tar.gz
##chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/execution/*/out/*-calls/SAMPLE_*/
#chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-PostprocessGermlineCNVCalls/shard-*/inputs/*/*
#
#rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/execution/*.preprocessed.annotated.tsv
#rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/inputs/*/*.preprocessed.interval_list
#rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/inputs/*/human_g1k_v37_fix.fasta
#rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-AnnotateIntervals/inputs/*/k100.umap.20260509.hg19.bed.gz
##rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-DetermineGermlineContigPloidyCohortMode/execution
#rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-DetermineGermlineContigPloidyCohortMode/inputs/*/*.hdf5
#rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-DetermineGermlineContigPloidyCohortMode/inputs/*/*.preprocessed.filtered.interval_list
#rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/inputs/*/*.hdf5
#rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/inputs/*/*.preprocessed.annotated.tsv
##rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/execution/*/*gcnv-calls*.tar.gz
##rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-GermlineCNVCallerCohortMode/shard-*/execution/*/out/*-calls/SAMPLE_*/
#rm cromwell-executions/CNVGermlineCohortWorkflow_02/${uuid}/call-PostprocessGermlineCNVCalls/shard-*/inputs/*/*
