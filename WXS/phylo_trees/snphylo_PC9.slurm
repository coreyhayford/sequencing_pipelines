#!/bin/bash
#SBATCH --mail-user=<corey.hayford@vanderbilt.edu>
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=snphylo_PC9_chrRenamed.out

ml GCC/6.4.0-2.28 OpenMPI/2.1.1 SNPhylo/20160204-Python-2.7.14
snphylo.sh -v samples_called_vars_named_chrRenamed.vcf -a 3366
