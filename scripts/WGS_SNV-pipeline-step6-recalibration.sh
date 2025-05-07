#!/bin/sh
#$ -S /bin/sh
set -o errexit

#Analysis directory
AD=`pwd`

# load modules
module load GATK/4.3.0.0
module load samtools/1.17

# Get parameters
PARFILE=$1
source ./$PARFILE
LOGFILE=${LOG_DIR}/${STEP6_LOG}

# run samples
while IFS= read -r fastq1 && IFS= read -r fastq2 <&3;
do
	# Prepare filenames file 1
	FASTQGZ_1=`basename ${fastq1}`
	FASTQ_1=${FASTQGZ_1%.*}
	SAMPLE_1=${FASTQ_1%.*}

	# Prepare filenames file 2
	FASTQGZ_2=`basename ${fastq2}`
	FASTQ_2=${FASTQGZ_2%.*}
	SAMPLE_2=${FASTQ_2%.*}

	# Get sample ID
	SAMPLE=${FASTQ_1%_R1*}

	# Make output directory
	OUTPUT_DIR=06-recalibration/${SAMPLE}
	if [[ ! -d ${OUTPUT_DIR} ]]; then
		mkdir -p -m 777 ${OUTPUT_DIR}
	fi	

	# Input dir
	INPUT_DIR=05-MarkDuplicates/${SAMPLE}

	# Number of threads
	NT=${SLURM_CPUS_PER_TASK}
	MEM=${SLURM_MEM_PER_NODE}
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`

	# Print header
	echo ""
	echo "#################################"
	echo "#" ${SAMPLE}
	echo "#" MEMORY PER NODE: ${MEM}
	echo "#" CPUS: ${NT}
	echo "# Mark duplicated reads start:"
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""

        ##############################################
        ###    Merge BAM files for recalibration   ###
        ##############################################

		OUTPUT_BAM=${OUTPUT_DIR}/dedup-${SAMPLE}.bam
		OUTPUT_RG_BAM=${OUTPUT_DIR}/dedup-${SAMPLE}-RG_FIX.bam
		SORTED_OUTPUT_RG_BAM=${OUTPUT_DIR}/sorted-dedup-${SAMPLE}.bam

        samtools merge -f -n -r -O BAM -c \
            -@ 8 ${OUTPUT_BAM} \
            ${INPUT_DIR}/dedup-${SAMPLE}-chr*.bam # Regular expression for grabing all input bam files for a given sample

        ##############################################
        ###    Add @RG SM line to output bam file  ###
        ##############################################
		
		echo ""
		echo FIX @RG FIELD
		echo ""

		samtools addreplacerg -r ID:test -r SM:test -r PL:ILLUMINA -o ${OUTPUT_RG_BAM} ${OUTPUT_BAM}
       
        samtools sort -@ 10 -T TMP -o ${SORTED_OUTPUT_RG_BAM} ${OUTPUT_RG_BAM}

        samtools index -@ 10 ${SORTED_OUTPUT_RG_BAM}
    
		rm -f ${OUTPUT_RG_BAM}

        #############################
        ###    Run recalibration  ###
        #############################

		echo ""
		echo "Base recalibration ..."
		echo ""

		echo ">>> BaseRecalibrator ..."
		gatk BaseRecalibrator \
			-R ${GENOME} \
			-I ${SORTED_OUTPUT_RG_BAM} \
			--known-sites ${KNOWN_SNP_SITES} \
			-O ${OUTPUT_DIR}/${SAMPLE}.BaseRecalibrator.table
			#-nct $NT

		echo ""
		echo ">>> RUN ApplyBQSR ..."
		echo ""

		echo ">>> Apply base recalibrator ..."

		RECALIBRATED_BAM=${OUTPUT_DIR}/recalRealignedReads-${SAMPLE}.bam

		gatk ApplyBQSR \
 		     -R ${GENOME} \
 		     -I ${SORTED_OUTPUT_RG_BAM} \
 		  	 --bqsr-recal-file ${OUTPUT_DIR}/${SAMPLE}.BaseRecalibrator.table \
  		     -O ${RECALIBRATED_BAM} \
			 --create-output-bam-index true
	     
		echo ""
		echo ">>> Rerun recalibration on recalibrated bam (After)"
		echo ""

		gatk BaseRecalibrator \
 		    -R ${GENOME} \
 		    -I ${RECALIBRATED_BAM} \
 		  	--known-sites ${KNOWN_SNP_SITES} \
 		  	-O ${OUTPUT_DIR}/${SAMPLE}.BaseRecalibratorAfter.table
  		    #-nct $NT \
			# -BQSR ${SAMPLE}.BaseRecalibrator 
  		    
		echo ""
		echo ">>> GENERATE PLOTS  ..." 
		echo ""

  		gatk AnalyzeCovariates \
   			-before ${OUTPUT_DIR}/${SAMPLE}.BaseRecalibrator.table \
 		    -after  ${OUTPUT_DIR}/${SAMPLE}.BaseRecalibratorAfter.table \
  			-plots  ${OUTPUT_DIR}/${SAMPLE}-realignedRecalPlots.pdf

	# Print footnote
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
	echo ""
	echo "#################################"
	echo "# Finished base recalibration ..."
	echo "#" ${SAMPLE}
	echo "#" $TIMESTAMP
	echo "#################################"
	echo "" 		 
		
done < $INPUTFILE_1 3<$INPUTFILE_2 > $LOGFILE 2>&1
