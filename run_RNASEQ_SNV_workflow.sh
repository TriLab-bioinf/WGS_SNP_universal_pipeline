#!/bin/bash
set -o errexit

# Set default configuration file to config.txt
# Otherwise you can change the name of the configuration file used to for example "my_custom_configuration_file.txt" by running this scripts like this:
# $ run_wgs_workflow.sh my_custom_configuration_file.txt

CONFIG=$1
if [[ -z $1 ]]; then
    CONFIG=config.txt
fi

# Checking files required for running jobs
if [[ ! -e ${CONFIG} ]]; then echo Configuration file ${CONFIG} not found in current directory!!; echo ""; exit; fi

source ${CONFIG}

if [[ ! -d ${LOG_DIR} ]]; then
    mkdir -m 777 ${LOG_DIR}
fi

# Check input files from config are OK
if [[ ! -e ${INPUTFILE_1} || ! -e ${INPUTFILE_2} ]]; then echo INPUT file with read file paths missing!!; echo ""; exit; fi
if [[ ! -e ${STAR_DB}/SA ]]; then echo STARdb missing!!; echo ""; exit; fi
if [[ ! -e ${GENOME} ]]; then echo Genome file missing!!; echo ""; exit; fi
if [[ ! -e ${KNOWN_SNP_SITES} ]]; then echo Known SNP sites file missing!!; echo ""; exit; fi

JOB_IDS=()

if [[ ${STEP1} == false ]]; then echo RNASEQ_SNV pipeline stops before STEP-1; exit; fi
# STEP-1: FASTQC
echo Submitting RNASEQ_SNV-pipeline-step1-fastqc.swarm
jobid1=$(swarm --logdir ./00-swarm-log --job-name step1-fastqc  --sbatch "--export=CONFIG_FILE=${CONFIG}" -f ./scripts/RNASEQ_SNV-pipeline-step1-fastqc.swarm)
JOB_IDS+=($jobid1)


if [[ ${STEP2} == false ]]; then echo RNASEQ_SNV pipeline stops before STEP-2; exit; fi
# STEP-2: MAPPING
echo Submitting RNASEQ_SNV-pipeline-step2-mapping.swarm
jobid2=$(swarm --logdir ./00-swarm-log --job-name step2-mapping --sbatch "--export=CONFIG_FILE=${CONFIG}" -f ./scripts/RNASEQ_SNV-pipeline-step2-mapping.swarm)
JOB_IDS+=($jobid2)

if [[ ${STEP3} == false ]]; then echo RNASEQ_SNV pipeline stops before STEP-3; exit; fi
# STEP-3: SPLIT CHR
echo Submitting RNASEQ_SNV-pipeline-step3-split.swarm
jobid3=$(swarm --dependency afterok:$jobid2 --logdir ./00-swarm-log --job-name step3-split --sbatch "--export=CONFIG_FILE=${CONFIG}" -f ./scripts/RNASEQ_SNV-pipeline-step3-split.swarm)
JOB_IDS+=($jobid3)

if [[ ${STEP4} == false ]]; then echo RNASEQ_SNV pipeline stops before STEP-4; exit; fi
# STEP-4: DELETE DISCORDANT READS
echo Submitting RNASEQ_SNV-pipeline-step4-cleanBAM.swarm
jobid4=$(swarm --dependency afterok:$jobid3 --logdir ./00-swarm-log --job-name step4-cleanBAM --sbatch "--export=CONFIG_FILE=${CONFIG}" -f ./scripts/RNASEQ_SNV-pipeline-step4-cleanBAM.swarm)
JOB_IDS+=($jobid4)

if [[ ${STEP5} == false ]]; then echo RNASEQ_SNV pipeline stops before STEP-5; exit; fi
# STEP-5: MARK DUPLICATES
echo Submitting RNASEQ_SNV-pipeline-step5-MarkDuplicates.swarm
jobid5=$(swarm --dependency afterok:$jobid4  --logdir ./00-swarm-log --job-name step5-MarkDuplicates --sbatch "--export=CONFIG_FILE=${CONFIG}" -f ./scripts/RNASEQ_SNV-pipeline-step5-MarkDuplicates.swarm)
JOB_IDS+=($jobid5)

if [[ ${STEP6} == false ]]; then echo RNASEQ_SNV pipeline stops before STEP-6; exit; fi
# STEP-6: BASE RECALIBRATION
echo Submitting RNASEQ_SNV-pipeline-step6-recalibration.swarm
jobid6=$(swarm --dependency afterok:$jobid5 --logdir ./00-swarm-log --job-name step6-recalibration --sbatch "--export=CONFIG_FILE=${CONFIG}" -f ./scripts/RNASEQ_SNV-pipeline-step6-recalibration.swarm)
JOB_IDS+=($jobid6)

if [[ ${STEP7} == false ]]; then echo RNASEQ_SNV pipeline stops before STEP-7; exit; fi
# STEP-7: SPLIT RECALIBRATED BAM FILES
echo Submitting RNASEQ_SNV-pipeline-step7-split.swarm
jobid7=$(swarm --dependency afterok:$jobid6 --logdir ./00-swarm-log --job-name step7-split --sbatch "--export=CONFIG_FILE=${CONFIG}" -f ./scripts/RNASEQ_SNV-pipeline-step7-split.swarm)
JOB_IDS+=($jobid7)

if [[ ${STEP8} == false ]]; then echo RNASEQ_SNV pipeline stops before STEP-8; exit; fi
# STEP-8: VARIANT CALLING
echo Submitting RNASEQ_SNV-pipeline-step8-VariantCalling.swarm
jobid8=$(swarm --dependency afterok:$jobid7 --logdir ./00-swarm-log --job-name step8-VariantCalling --sbatch "--export=CONFIG_FILE=${CONFIG}" -f ./scripts/RNASEQ_SNV-pipeline-step8-VariantCalling.swarm)
JOB_IDS+=($jobid8)

if [[ ${STEP9} == false ]]; then echo RNASEQ_SNV pipeline stops before STEP-9; exit; fi
# STEP-9: FILTER VARIANTS
echo Submitting RNASEQ_SNV-pipeline-step9-filterVariants.swarm
jobid9=$(swarm --dependency afterok:$jobid8 --logdir ./00-swarm-log --job-name step9-filterVariants --sbatch "--export=CONFIG_FILE=${CONFIG}" -f ./scripts/RNASEQ_SNV-pipeline-step9-filterVariants.swarm)
JOB_IDS+=($jobid9)

# Print out job ids into workflow_processes.txt
rm -f workflow_processes.txt
for jobid in ${JOB_IDS[@]}; do 
    echo $jobid >> workflow_processes.txt
done    
