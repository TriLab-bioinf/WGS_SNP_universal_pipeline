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
LOGFILE=${LOG_DIR}/${STEP9_LOG}


# Input dir
INPUT_DIR=08-variant-calling

# Output dir
OUTPUT_DIR=09-filtered_variants
if [[ ! -d ${OUTPUT_DIR} ]]; then
	mkdir -m 777 ${OUTPUT_DIR}
fi	

# Merge all per-chr VCF files into a single VCF
PER_CHR_VCF=()
for i in {10..19} {1..9} M X Y; do 
	PER_CHR_VCF+="-I ${INPUT_DIR}/all_samples-chr${i}.vcf "
done	

echo ""
echo Merging vcf files ...
echo ${PER_CHR_VCF} 
echo ""

MERGED_VCF=${INPUT_DIR}/all_samples.vcf

gatk GatherVcfs ${PER_CHR_VCF} -O ${MERGED_VCF} > ${LOGFILE} 2>&1

MERGED_SNP_VCF=${INPUT_DIR}/all_samples-SNP.vcf
MERGED_INDEL_VCF=${INPUT_DIR}/all_samples-INDEL.vcf

# Subset to SNPs-only callset with SelectVariants
echo ""
echo Subseting SNPs from merged vcf file ${MERGED_VCF} ...
echo ""

gatk SelectVariants \
    -V ${MERGED_VCF} \
    -select-type SNP \
    -O ${MERGED_SNP_VCF}  >> ${LOGFILE} 2>&1

# Subset to indels-only callset with SelectVariants
echo ""
echo Subseting INDELs from merged vcf file ...
echo ""

gatk SelectVariants \
    -V ${MERGED_VCF} \
    -select-type INDEL \
    -O ${MERGED_INDEL_VCF} >> ${LOGFILE} 2>&1

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
echo "# Filtering variants start:"
echo "#" $TIMESTAMP
echo "#################################"
echo ""


echo ""
echo "FILTERING SNP VARIANTS ..."
echo ""

FILTERED_SNP_VCF=${OUTPUT_DIR}/all_samples-SNP-filtered.vcf

# Filter based on GATK recomendations:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantFiltration \
    --R ${GENOME} \
    --V ${MERGED_SNP_VCF} \
    --window 35 \
    --cluster 3 \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${FILTERED_SNP_VCF} >> ${LOGFILE} 2>&1

echo ""
echo "FILTERING INDEL VARIANTS ..."
echo ""

FILTERED_INDEL_VCF=${OUTPUT_DIR}/all_samples-INDEL-filtered.vcf

# Filter based on GATK recomendations:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantFiltration \
    --R ${GENOME} \
    --V ${MERGED_INDEL_VCF} \
    --window 35 \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O ${FILTERED_INDEL_VCF} >> ${LOGFILE} 2>&1

e# Print footnote
TIMESTAMP=`date "+%Y-%m-%d %H:%M:%S"`
echo ""
echo "#################################"
echo "# Finished filtering variants ..."
echo "#" across samples
echo "#" $TIMESTAMP
echo "#################################"
echo ""	

