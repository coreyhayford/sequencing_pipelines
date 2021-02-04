## How to analyze ATAC-seq data ##
1. Run pipelineATAC_example... in order of part1, part2, and then part3
1.1 part 1 = alignment; part2 = trimming; part3 = peak counting
2. To compare library complexities before and after trimming, run ATAC_library complexity.slurm
2.1 Run preseq_curves.R to compare libraries based on actual and projected depth
3. To shift reads for coverage analysis, run ATACseq_shift_example.slurm, which calls ATAC_shift_example.R
4. To profile various sites, run ATACseq_cutSites_examples.R
5. For downstream analyses, run ATAC_analysis.R
5.1 In runConsensusRegions section, run modified soGGi_runConsensusRegions_fixed.R to prevent errors

