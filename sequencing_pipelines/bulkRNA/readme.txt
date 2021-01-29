### How to analyze RNAseq data ###
1. Fastqs in the raw_data folder (VU1file)
2. rnaSeq_analysis_test_S1.slurm is an  example script (in scripts folder)
2.1  rnaSeq_analysis_test_each.slurm were used for other samples independently
3. Data all stored in SampleX folder (VU1file)
4. Final output is ...featureCounts.txt file (in SampleX folder in VU1file)
4.1 Outputs are read into R and cbind together to make count matrix
5. Variant (mutation) calling pipeline (adapted from GATK unpublished RNA-seq variant calling pipeline) is in variantCalling folder
