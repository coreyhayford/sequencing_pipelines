#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH	--cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem=64G
#SBATCH --output=ATAC1_LibComp_part2.out

cd /data/lola/hayforc/ATACseq/SKMEL5/raw_data

ml Anaconda2/5.0.1
source activate RNAvelocityEnvReal

preseq lc_extrap -pe -verbose -o trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_lcExtrapResults.txt -bam trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted.bam
preseq c_curve -pe -verbose -o trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_cCurveResults.txt -bam trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted.bam

ml picard/2.17.10

java -jar $EBROOTPICARD/picard.jar EstimateLibraryComplexity I=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted.bam O=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_ELC.txt VALIDATION_STRINGENCY=LENIENT
java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics I=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted.bam H=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_histogram.pdf O=trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_ISM.txt VALIDATION_STRINGENCY=LENIENT
