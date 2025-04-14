#!/bin/sh
#$ -S /bin/sh
set -o errexit

#Analysis directory
AD=`pwd`

# load modules
module load picard/2.9.2

# Get parameters
PARFILE=$1
CHR=$2
source ./$PARFILE
LOGFILE=${LOG_DIR}/${CHR}_${STEP5_LOG}



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
	OUTPUT_DIR=05-MarkDuplicates/${SAMPLE}
	if [[ ! -d ${OUTPUT_DIR} ]]; then
		mkdir -p -m 777 ${OUTPUT_DIR}
	fi	

	# Input dir
	INPUT_DIR=04-cleanBAM/${SAMPLE}

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

	INPUT_BAM=${INPUT_DIR}/sorted-clean-${SAMPLE}-${CHR}.bam
	OUTPUT_BAM=${OUTPUT_DIR}/dedup-${SAMPLE}-${CHR}.bam

	java -Xmx60g -XX:ParallelGCThreads=20 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
	INPUT=${INPUT_BAM} OUTPUT=${OUTPUT_BAM} \
	METRICS_FILE=metrics-${SAMPLE}.txt

	echo "CREATE BAM INDEX ..."
	java -Xmx60g -jar $PICARDJARPATH/picard.jar BuildBamIndex INPUT=${OUTPUT_BAM}
	
	# Print footnote
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
	echo ""
	echo "#################################"
	echo "# Finished marking duplicated reads ..."
	echo "#" ${SAMPLE}
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""
		
done < $INPUTFILE_1 3<$INPUTFILE_2 > $LOGFILE 2>&1
