#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --output=ATAC_shift_UT.out

cd /data/lola/hayforc/ATACseq/SKMEL5/preprocessed

ml GCC/8.2.0  
ml OpenMPI/3.1.4
ml R-bundle-Bioconductor/3.9

Rscript --no-save ATAC_shift_UT.R
