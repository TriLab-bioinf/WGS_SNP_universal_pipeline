#!/bin/sh
#$ -S /bin/sh
set -o errexit

# Analysis directory
AD=`pwd`

# Load modules 
module load fastqc/0.11.5

# Get parameters
PARFILE=$1
source ./$PARFILE
LOGFILE=${LOG_DIR}/${STEP1_LOG}

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
	OUTPUT_DIR=01-fastqc/${SAMPLE}
	if [[ ! -d ${OUTPUT_DIR} ]]; then
		mkdir -p -m 777 ${OUTPUT_DIR}
	fi	
	
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
	echo "# Fastqc start:"
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""
	

	fastqc -o ${OUTPUT_DIR} ${fastq1}
	fastqc -o ${OUTPUT_DIR} ${fastq2}
	
	# Print footnote
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
	echo ""
	echo "#################################"
	echo "# Finished fastqc ..."
	echo "#" ${SAMPLE}
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""
	
		
done < $INPUTFILE_1 3<$INPUTFILE_2 > $LOGFILE 2>&1
