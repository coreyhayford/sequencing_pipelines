#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --output=pipelineWXS_combineFilter_3.out

REFERENCE='/scratch/hayforc/Ref_genome/Homo_sapiens_assembly38.fasta'
MILLS='/scratch/hayforc/known_variants/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
DBSNP='/scratch/hayforc/known_variants/Homo_sapiens_assembly38.dbsnp138.vcf'

module load GATK/3.8-0
# call using java -jar $EBROOTGATK/GenomeAnalysisTK.jar

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $REFERENCE \
-V Sample1_called_vars_named.g.vcf \
-V Sample2_called_vars_named.g.vcf \
-V Sample3_called_vars_named.g.vcf \
-V Sample4_called_vars_named.g.vcf \
-V Sample5_called_vars_named.g.vcf \
-V Sample6_called_vars_named.g.vcf \
-V Sample7_called_vars_named.g.vcf \
-V Sample8_called_vars_named.g.vcf \
-o samples_called_vars_named.vcf 

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REFERENCE \
-V samples_called_vars_named.vcf \
--selectTypeToInclude SNP \
-o samples_raw_snps_named.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REFERENCE \
-V samples_called_vars_named.vcf \
--selectTypeToInclude INDEL \
-o samples_raw_indels_named.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $REFERENCE \
-V samples_raw_snps_named.vcf \
-filter "QD < 2.0" --filterName "QD2" \
-filter "QUAL < 30.0" --filterName "QUAL30" \
-filter "SOR > 3.0" --filterName "SOR3" \
-filter "FS > 60.0" --filterName "FS60" \
-filter "MQ < 40.0" --filterName "MQ40" \
-filter "MQRankSum < -12.5" --filterName "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum-8" \
-o samples_filtered_snps_named.vcf

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $REFERENCE \
-V samples_raw_indels_named.vcf \
-filter "QD < 2.0" --filterName "QD2" \
-filter "QUAL < 30.0" --filterName "QUAL30" \
-filter "FS > 200.0" --filterName "FS200" \
-filter "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum-20" \
-o samples_filtered_indels_named.vcf
