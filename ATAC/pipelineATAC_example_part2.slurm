#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --output=pipelineATAC_sample1_part1_LO.out

cd /data/lola/hayforc/ATACseq/SKMEL5/raw_data

ml GCC/6.4.0-2.28
ml Intel/2017.4.196
ml picard/2.17.10
ml SAMtools/1.6

#java -jar $EBROOTPICARD/picard.jar SortSam I=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned.sam O=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted.sam SO=coordinate VALIDATION_STRINGENCY=LENIENT
#samtools view -S -b trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted.sam > trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted.bam
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted.bam O=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup.bam M=3334-CH-1-TAGGCATG-ACTGCATA_S121_metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics I=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup.bam O=3334-CH-1-TAGGCATG-ACTGCATA_S121_ISM.txt H=3334-CH-1-TAGGCATG-ACTGCATA_S121_metrics.pdf

samtools view -b -q 10 trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup.bam > trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique.bam
#samtools view -S -b trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique.sam > trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique.bam
samtools index trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique.bam


