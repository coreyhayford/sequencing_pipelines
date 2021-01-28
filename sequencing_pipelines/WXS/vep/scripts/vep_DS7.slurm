#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --output=vep6.out
#SBATCH --cpus-per-task=8

module load Anaconda3
source activate RNAvelocityEnvReal

export VEP_PATH=/data/lola/hayforc/WXS/vep
export VEP_DATA=/data/lola/hayforc/WXS/.vep
cd $VEP_PATH
export PERL5LIB=$VEP_PATH:$PERL5LIB

./vep \
--verbose \
--everything \
--species homo_sapiens \
-i /data/lola/hayforc/WXS/PC9/processed_data/all_vars/sample6_vars.vcf.gz \
-o /data/lola/hayforc/WXS/PC9/processed_data/all_vars/vep_sample6.txt \
-sf /data/lola/hayforc/WXS/PC9/processed_data/all_vars/vep_summary_sample6.html \
--warning_file /data/lola/hayforc/WXS/PC9/processed_data/all_vars/vep_warnings_sample6.txt \
--fork $SLURM_CPUS_ON_NODE \
--dir $VEP_DATA \
--assembly GRCh38 \
--offline \
--fasta /data/lola/hayforc/WXS/bwa_ref_genome/Homo_sapiens_assembly38.fasta
