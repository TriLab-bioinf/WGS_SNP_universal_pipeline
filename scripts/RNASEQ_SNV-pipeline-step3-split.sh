#!/bin/sh
#$ -S /bin/sh
set -o errexit

#Analysis directory
AD=`pwd`

# load modules 
module load samtools/1.17

# Get parameters
PARFILE=$1
source ./$PARFILE
LOGFILE=${LOG_DIR}/${STEP3_LOG}

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
	OUTPUT_DIR=03-split-bam/${SAMPLE}
	if [[ ! -d ${OUTPUT_DIR} ]]; then
		mkdir -p -m 777 ${OUTPUT_DIR}
	fi	

	# Input dir
	INPUT_DIR=02-mapping/${SAMPLE}

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
	
	INPUT_BAM=${INPUT_DIR}/sorted-${SAMPLE}.bam
	

	for counter in M {1..22} X Y; do

		CHR=chr${counter}
		
		OUTPUT_BAM=${OUTPUT_DIR}/sorted-${SAMPLE}-${CHR}.bam

		echo Processing ${CHR}

		samtools view ${INPUT_BAM} $CHR -hb | samtools sort -@ 15 -T ${OUTPUT_DIR}/TMP_${CHR} -o ${OUTPUT_BAM} -
			
		samtools index ${OUTPUT_BAM}

	done

	# Print footnote
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
	echo ""
	echo "#################################"
	echo "# Finished splitting by chromosome ..."
	echo "#" ${SAMPLE}
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""

 done < $INPUTFILE_1 3<$INPUTFILE_2 > $LOGFILE 2>&1
