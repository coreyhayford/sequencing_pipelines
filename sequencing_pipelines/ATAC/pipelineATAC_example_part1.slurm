#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --ntasks=1
#SBATCH --time=36:00:00
#SBATCH --mem=32G
#SBATCH --output=pipelineATAC_sample1.out

cd /data/lola/hayforc/ATACseq/SKMEL5/raw_data

R1=$1
R2=$2
REFERENCE='/data/lola/hayforc/WXS/bwa_ref_genome/Homo_sapiens_assembly38.fasta'

ml GCC/6.4.0-2.28
ml Intel/2017.4.196
ml cutadapt/1.16-Python-3.6.3

cutadapt -a TAGGCATG -A TATGCAGT -o trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_R1_001.fastq.gz -p trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_R2_001.fastq.gz $R1 $R2


ml BWA/0.7.17

bwa mem $REFERENCE trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_R1_001.fastq.gz trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_R2_001.fastq.gz -M -t $SLURM_NTASKS_PER_NODE > trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned.sam
