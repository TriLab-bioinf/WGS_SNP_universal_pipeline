#!/bin/bash
for dir in 01-fastqc 02-trimming 02-mapping 03-split-bam 04-cleanBAM 05-MarkDuplicates 06-recalibration 07-split-recalibration 08-variant-calling 09-filtered_variants 00-swarm-log 00-workflow-log workflow_processes.txt; do
    echo Deleting ... ${dir}
    rm -Rf ${dir}
done
