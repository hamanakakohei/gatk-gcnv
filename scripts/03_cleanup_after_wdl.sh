#!/usr/bin/env bash
set -euo pipefail
uuid=$1


#ls cromwell-executions/CNVGermlineCaseWorkflow/${uuid}/call-DetermineGermlineContigPloidyCaseMode/inputs/*/*.hdf5
#ls cromwell-executions/CNVGermlineCaseWorkflow/${uuid}/call-GermlineCNVCallerCaseMode/shard-*/inputs/*/*.hdf5
#ls cromwell-executions/CNVGermlineCaseWorkflow/${uuid}/call-PostprocessGermlineCNVCalls/shard-*/inputs/*/*
#
#chmod 644 cromwell-executions/CNVGermlineCaseWorkflow/${uuid}/call-DetermineGermlineContigPloidyCaseMode/inputs/*/*.hdf5
#chmod 644 cromwell-executions/CNVGermlineCaseWorkflow/${uuid}/call-GermlineCNVCallerCaseMode/shard-*/inputs/*/*.hdf5
#chmod 644 cromwell-executions/CNVGermlineCaseWorkflow/${uuid}/call-PostprocessGermlineCNVCalls/shard-*/inputs/*/*

rm cromwell-executions/CNVGermlineCaseWorkflow/${uuid}/call-DetermineGermlineContigPloidyCaseMode/inputs/*/*.hdf5
rm cromwell-executions/CNVGermlineCaseWorkflow/${uuid}/call-GermlineCNVCallerCaseMode/shard-*/inputs/*/*.hdf5
rm cromwell-executions/CNVGermlineCaseWorkflow/${uuid}/call-PostprocessGermlineCNVCalls/shard-*/inputs/*/*
