#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --output=pipelineATAC_sample1_part2.out

cd /data/lola/hayforc/ATACseq/SKMEL5/preprocessed

ml GCC/6.4.0-2.28
ml Intel/2017.4.196
ml SAMtools/1.6

#samtools view trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique.bam | grep -v chrM | samtools sort -o trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique_rmChrM.bam
samtools view -b -f 3 -F 4 -F 8 -F 256 -F 2048 -q 30 trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique.bam > trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique_fullClean.bam

ml GCC/6.4.0-2.28  
ml OpenMPI/2.1.1
ml Intel/2017.4.196
ml IntelMPI/2017.3.196
ml MACS2/2.1.2.1-Python-2.7.14

macs2 callpeak -t trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique_fullClean.bam -f BAMPE -n untreated -g hs -q 0.05 --nomodel --shift -100 --extsize 200 -B --broad --outdir untreated_peaks

