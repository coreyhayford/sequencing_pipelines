#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=32G
#SBATCH --output=pipeline.out

R1=$1
R2=$2
REFERENCE='/scratch/hayforc/Ref_genome/Homo_sapiens_assembly38.fasta'
MILLS='/scratch/hayforc/known_variants/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
DBSNP='/scratch/hayforc/known_variants/Homo_sapiens_assembly38.dbsnp138.vcf'

module load GCC/5.4.0-2.26
module load Intel/2016.3.210
module load STAR/2.5.2b

STAR --genomeDir /scratch/hayforc/STAR_ref_genome --runThreadN 8 \
--readFilesCommand gunzip -c --twopassMode Basic \
--readFilesIn $R1 $R2

module load GCC/6.4.0-2.28
module load Intel/2017.4.196
module load SAMtools/1.6

samtools view -bS Aligned.out.sam | samtools sort -o aligned_sorted.bam

module load picard/2.17.10
# call using java -jar $EBROOTPICARD/picard.jar

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=aligned_sorted.bam \
    O=rg_added_sorted.bam \
    SO=coordinate \
    RGID=id \
    RGLB=library \
    RGPL=illumina \
    RGPU=machine \
    RGSM=sample

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=rg_added_sorted.bam \
    O=deduped.bam  \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M=output.metrics

java -jar $EBROOTPICARD/picard.jar ReorderSam \
    I=deduped.bam \
    O=reordered.bam \
    R=$REFERENCE \
    CREATE_INDEX=TRUE

module load GATK/3.8-0
# call using java -jar $EBROOTGATK/GenomeAnalysisTK.jar

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T SplitNCigarReads \
    -R $REFERENCE \
    -I reordered.bam \
    -o split.bam \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS

samtools view -bS split.bam | samtools sort -o split_sorted.bam
samtools index split_sorted.bam

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R $REFERENCE \
    -I split_sorted.bam \
    -o realigner.intervals

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R $REFERENCE \
    -I split_sorted.bam \
    -known $MILLS \
    -targetIntervals realigner.intervals \
    -o realigned.bam

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R $REFERENCE \
    -I realigned.bam \
    -knownSites $DBSNP \
    -knownSites $MILLS \
    -o recal.table

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T PrintReads  \
    -R $REFERENCE \
    -I realigned.bam \
    -BQSR recal.table \
    -o recal.bam

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R $REFERENCE \
    -I recal.bam \
    -dontUseSoftClippedBases \
    -stand_call_conf 20.0 \
    -o called_vars.vcf \

java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R $REFERENCE \
    -V called_vars.vcf \
    -window 35 \
    -cluster 3 \
    -filterName FS \
    -filter "FS > 30.0" \
    -filterName QD \
    -filter "QD < 2.0" \
    -o filtered_vars.vcf

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