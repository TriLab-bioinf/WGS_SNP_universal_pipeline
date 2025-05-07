#!/bin/sh
#$ -S /bin/sh
set -o errexit

#Analysis directory
AD=`pwd`

# load modules 
module load bbtools/39.06 
module load bwa-mem2/2.2.1 
module load GATK/4.3.0.0
module load samtools/1.17 

# Get parameters
PARFILE=$1
source ./$PARFILE
LOGFILE=${LOG_DIR}/${STEP2_LOG}

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
	OUTPUT_DIR=02-trimming/${SAMPLE}
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
	echo "# Splitting by chromosome start:"
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""

	#############################
	###       trimming        ###
	#############################

	echo ""
	echo "Read QC and trimming using BBMAP ..."
	echo ""

	CLEAN_FASTQ1=${OUTPUT_DIR}/clean-${FASTQGZ_1}
	CLEAN_FASTQ2=${OUTPUT_DIR}/clean-${FASTQGZ_2}

	bbtools bbduk in=${fastq1} in2=${fastq2} \
		out1=${CLEAN_FASTQ1} out2=${CLEAN_FASTQ2} \
		trimq=10 maq=10 t=${NT} overwrite=t \
		bhist=${OUTPUT_DIR}/${SAMPLE}.bhist qhist=${OUTPUT_DIR}/${SAMPLE}.qhist \
		lhist=${OUTPUT_DIR}/${SAMPLE}.lhist


	# Print footnote
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
	echo ""
	echo "#################################"
	echo "# Finished trimming ..."
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""
	
	#############################
	###       alignment       ###
	#############################
	
	# Make output directory
	MAP_OUTPUT_DIR=02-mapping/${SAMPLE}
	if [[ ! -d ${MAP_OUTPUT_DIR} ]]; then
		mkdir -p -m 777 ${MAP_OUTPUT_DIR}
	fi	

	INPUT_DIR=02-trimming/${SAMPLE}

	echo ""
	echo " Mapping reads with BWA MEM2 ..."
	echo ""

    RAW_OUTPUT_BAM=${MAP_OUTPUT_DIR}/raw-sorted-${SAMPLE}.out.bam

    bwa-mem2 mem \
        -t 30 \
        ${BWA_DB} ${CLEAN_FASTQ1} ${CLEAN_FASTQ2} | \
        samtools view -hb - | \
        samtools sort -@ 16 \
            -O BAM -T ${SAMPLE}.tmp \
            -o ${RAW_OUTPUT_BAM} - > {log.logfile} 2>&1

        samtools index -@ 8 ${OUTPUT_BAM}   

	echo ""
	echo "Spllitting N cigars ..."
	echo ""
	
	INPUT_BAM=${RAW_OUTPUT_BAM}
	OUTPUT_BAM=${MAP_OUTPUT_DIR}/sorted-${SAMPLE}.bam

	gatk SplitNCigarReads \
		-R ${GENOME} \
		-I ${INPUT_BAM} \
		-O ${OUTPUT_BAM} \
		--process-secondary-alignments true

	rm -f ${RAW_OUTPUT_BAM}

	samtools index ${OUTPUT_BAM}

	# Print footnote
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
	echo ""
	echo "#################################"
	echo "# Finished mapping ..."
	echo "#" ${SAMPLE}
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""
		
done < $INPUTFILE_1 3<$INPUTFILE_2 > $LOGFILE 2>&1
