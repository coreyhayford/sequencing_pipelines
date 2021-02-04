#!/bin/bash

#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --output=RNAseq_test_Sample1_.out

module load Anaconda2/5.0.1
source activate seqEnv

cd /data/lola/hayforc/RNAseq/
samplenum=$5
mkdir ${samplenum}
cd ${samplenum}/

# This script takes a fastq file of RNA-Seq data, runs FastQC and outputs a counts file for it.
# USAGE: sh rnaseq_analysis.sh <name of fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
fq1=$1
fq2=$2

# grab base of filename for naming outputs
base1=`basename $fq1 .subset.fq`
base2='basename $fq2 .subset.fq'
echo "Sample name is $base1"
echo "Sample name is $base2"

# specify the number of cores to use
# cores=6

# directory with genome reference FASTA and index files + name of the gene annotation file
genome=/scratch/hayforc/STAR_ref_genome/
gtf=/data/lola/hayforc/scRNAseq/genes.gtf

# make all of the output directories
# The -p option means mkdir will create the whole path if it
# does not exist and refrain from complaining if it does exist

#mkdir -p ~/data/lola/hayforc/RNAseq/results/fastqc
#mkdir -p ~/data/lola/hayforc/RNAseq/results/STAR
#mkdir -p ~/data/lola/hayforc/RNAseq/results/counts

# set up output filenames and locations
#fastqc_out=~/data/lola/hayforc/RNAseq/results/fastqc
#align_out=~/data/lola/hayforc/RNAseq/results/STAR/${base1}_
#counts_input_bam=~/data/lola/hayforc/RNAseq/results/STAR/${base1}_Aligned.sortedByCoord.out.bam
#counts=~/data/lola/hayforc/RNAseq/results/counts/${base1}_featurecounts.txt
align_out=/data/lola/hayforc/RNAseq/Sample1/${base1}_
counts_input_bam=/data/lola/hayforc/RNAseq/Sample1/${base1}_Aligned.sortedByCoord.out.bam
counts=/data/lola/hayforc/RNAseq/Sample1/${base1}_featurecounts.txt

echo "Processing file $fq1"

# Run FastQC and move output to the appropriate folder
#export PATH=/home/hayforc/.conda/envs/seqEnv/fastqc:$PATH
#fastqc $fq1
#fastqc $fq2
#echo "FASTQC done"

R1=$3
R2=$4

# Run STAR
module load GCC/5.4.0-2.26
module load Intel/2016.3.210
module load STAR/2.5.2b
STAR --genomeDir $genome \
 --readFilesIn $R1 $R2 \
--outFileNamePrefix $align_out \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
 --readFilesCommand gunzip -c \
 --runThreadN 8
echo "Alignment done"

# Create BAM index
module load GCC/6.4.0-2.28
module load Intel/2017.4.196
module load SAMtools/1.6
samtools index $counts_input_bam
echo "Indexing done"

# Count mapped reads
export PATH=/home/hayforc/.conda/envs/seqEnv/bin/featureCounts:$PATH
featureCounts -p -t exon -g gene_id -a $gtf -o $counts -M $counts_input_bam
echo "Gene counting done"



# this `for` loop, will take the fastq files as input and run the script for all of them one after the other. 
# cd /data/lola/hayforc/RNAseq/
# chmod +x rnaSeq_analysis.sh
# sh rnaSeq_analysis.sh 2959-CH-1-AAGACCGT-GTCGATTG_S10_R1_001.fastq.gz 2959-CH-1-AAGACCGT-GTCGATTG_S10_R2_001.fastq.gz

# for fq in ~/unix_lesson/rnaseq/raw_data/*.fq # change to folder with data
# do
#   echo "running analysis on $fq"
#   rnaSeq_analysis.sh $fq
# done
