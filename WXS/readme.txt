### How to produce variants from Whole Exome Sequencing data ###
1. Fastq files stored in raw_data on VU1file
2. Initial slurm processing scripts stored in samples (run these on fastqs in raw_data - there is one for each sample)
3. Intermediate files (gVCFs) are stored in each of the SampleX directories (in VU1file)
4. Copy gVCF files and rename for each sample using renameSample.txt (in downstream) command to add sample information 
4.1 Each variant file SampleX_called_vars_... - named SampleX_called_vars_named_... (in VU1file)
5. Use pipelineWXS_combineFilter_3.slurm to make combined VCF file (in VU1file)
6. Use tabix to bgzip (bgzip -c file.vcf > file.vcf.gz) and tabix index (tabix -p vcf file.vcf.gz) the files
7. Use vcftools to separate to samples again for later analyses
	vcftools --gzvcf file.vcf.gz (sample statistics)
	vcf-compare file1.vcf.gz file2.vcf.gz etc.
8. QC results (from MultiQC) is in QC folder
9. Phylogenetic tree prediction is in phylo_trees folder
10. Scripts to annotate variants with Ensembl Variant Effect Predictor (VEP) database are in vep folder 
