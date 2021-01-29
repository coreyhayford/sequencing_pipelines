#!/bin/bash
#SBATCH --mail-user=corey.hayford@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=velo_untreat.txt

cd /data/lola/hayforc/scRNAseq/1614-1/
module load Anaconda3/5.0.1
module load GCC/6.4.0-2.28
module load Intel/2017.4.196
module load SAMtools/1.6
source activate RNAvelocityEnvReal
velocyto run10x /data/lola/hayforc/scRNAseq/1614-1/ /data/lola/hayforc/scRNAseq/genes.gtf
