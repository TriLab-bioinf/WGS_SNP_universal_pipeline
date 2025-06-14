#!/bin/sh
#$ -S /bin/sh
set -o errexit

#Analysis directory
AD=`pwd`

# load modules 
module load GATK/4.3.0.0

# Get parameters
PARFILE=$1
source ./$PARFILE
LOGFILE=${LOG_DIR}/${STEP8_LOG}
CHR=$2

# Save sample names into MY_SAMPLES array 
MY_SAMPLES=()

# Generate gvcf files for each sample
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
	OUTPUT_DIR=08-variant-calling/${SAMPLE}
	if [[ ! -d ${OUTPUT_DIR} ]]; then
		mkdir -p -m 777 ${OUTPUT_DIR}
	fi	

	# Input dir
	INPUT_DIR=07-split-recalibration/${SAMPLE}

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

	echo ""
	echo "VARIANT CALLING  ..."
	echo ""

	INPUT_BAM=${INPUT_DIR}/sorted-recalRealignedReads-${SAMPLE}-${CHR}.bam
	OUTPUT_GVCF=${OUTPUT_DIR}/${SAMPLE}-${CHR}.g.vcf
	MY_SAMPLES+=("-V ${OUTPUT_GVCF}")

	gatk HaplotypeCaller \
		-R ${GENOME} \
		-I ${INPUT_BAM} \
		-O ${OUTPUT_GVCF} \
		-ERC GVCF \
		-L $CHR \
		--standard-min-confidence-threshold-for-calling 30 \
		--native-pair-hmm-threads 4 \
		-dont-use-soft-clipped-bases
		# --genotyping_mode DISCOVERY is default in gatk version > 4.0 

	
	# Print footnote
	TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
	echo ""
	echo "#################################"
	echo "# Finished gvcf ..."
	echo "#" for sample ${SAMPLE} on ${CHR}
	echo "#" $TIMESTAMP
	echo "#################################"
	echo ""	
		
		
done < $INPUTFILE_1 3<$INPUTFILE_2 > $LOGFILE 2>&1

#################################################################
###    Call Variants within chromosomes across all samples    ###
#################################################################

# STEP-1 Combine GVCFs into a single VCF file

OUTPUT_COMBINED_VCF=08-variant-calling/all_samples-${CHR}.g.vcf

gatk CombineGVCFs \
    -R ${GENOME} \
    ${MY_SAMPLES[@]} \
    -O ${OUTPUT_COMBINED_VCF}

echo ""
echo "#################################"
echo "# Finished merging VCF files"
echo "#################################"
echo ""

# STEP-2 Call Variants within chromosomes across all samples


OUTPUT_GenotypeGVCFs=08-variant-calling/all_samples-${CHR}.vcf

gatk --java-options "-Xms2G -Xmx2G -XX:ParallelGCThreads=2" GenotypeGVCFs \
     -R ${GENOME} \
     --variant ${OUTPUT_COMBINED_VCF} \
	 --standard-min-confidence-threshold-for-calling 20 \
     -O ${OUTPUT_GenotypeGVCFs}   

# Print footnote
TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
echo ""
echo "#################################"
echo "# Finished calling variants ..."
echo "#" on ${CHR} across samples
echo "#" $TIMESTAMP
echo "#################################"
echo ""

