#!/bin/sh
#$ -S /bin/sh
set -o errexit

#Analysis directory
AD=`pwd`

# load modules 
module load picard/2.9.2
module load samtools/1.17

# Get parameters
PARFILE=$1
source ./$PARFILE
LOGFILE=${LOG_DIR}/${STEP4_LOG}


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
	OUTPUT_DIR=04-cleanBAM/${SAMPLE}
	if [[ ! -d ${OUTPUT_DIR} ]]; then
		mkdir -p -m 777 ${OUTPUT_DIR}
	fi	

	# Input dir
	INPUT_DIR=03-split-bam/${SAMPLE}

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
	echo "# Splitting by chromosome start:"
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""

	# Iterate through chromosomes
	for counter in {10..19} {1..9} M X Y; do
		

		CHR=chr${counter}
		INPUT_BAM=${INPUT_DIR}/sorted-${SAMPLE}-${CHR}.bam
		CLEAN_SAM=${OUTPUT_DIR}/clean-${SAMPLE}-${CHR}.sam
		SORTED_CLEAN_BAM=${OUTPUT_DIR}/sorted-clean-${SAMPLE}-${CHR}.bam

		# -F 256 => removes secondary alignments
		samtools view -F 256 -h ${INPUT_BAM} | awk '($7 == "=" || $1 ~ /^@/) &&  $6 !~ /N/' > ${CLEAN_SAM}

		java -Xmx60g -jar $PICARDJARPATH/picard.jar SortSam INPUT=${CLEAN_SAM} OUTPUT=${SORTED_CLEAN_BAM}  SORT_ORDER=coordinate
				
		samtools index ${SORTED_CLEAN_BAM}
		
		echo "Removing SAM file"
		rm ${CLEAN_SAM}
		
	done
	
	# Print footnote
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
	echo ""
	echo "#################################"
	echo "# Finished cleaning BAM ..."
	echo "#" ${SAMPLE}
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""


 done < $INPUTFILE_1 3<$INPUTFILE_2 > $LOGFILE 2>&1
