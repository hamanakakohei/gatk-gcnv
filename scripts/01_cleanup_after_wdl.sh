#!/usr/bin/env bash
set -euo pipefail
uuid=$1

ls cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.bam
ls cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.bai
ls cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.cram
ls cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.crai
ls cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.fasta

ls cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.preprocessed.interval_list
chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.bam
chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.bai
chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.cram
chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.crai
chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.fasta
chmod 644 cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.preprocessed.interval_list

rm cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.bam
rm cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.bai
rm cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.cram
rm cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.crai
rm cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.fasta
rm cromwell-executions/CNVGermlineCohortWorkflow_01/${uuid}/call-CollectCounts/shard-*/inputs/*/*.preprocessed.interval_list
