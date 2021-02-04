setwd('/data/lola/hayforc/ATACseq/SKMEL5/preprocessed')

.libPaths("~/R/rlib-3.6.0/")
library(S4Vectors)
library(GenomicAlignments)
library(GenomicRanges)
library(Biostrings)
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)

CTCF <- query(MotifDb, c("CTCF"))
ctcfMotif <- CTCF[[1]]
myRes <- matchPWM(ctcfMotif,BSgenome.Hsapiens.UCSC.hg38)
mainChrom <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
               "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", 
               "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")
# toCompare <- GRanges(mainChrom,ranges(myRes))


BAM_UT <- "trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique_fullClean.bam"
atacReads_Open_UT <- readGAlignmentPairs(BAM_UT)
read1_UT <- first(atacReads_Open_UT)
read2_UT <- second(atacReads_Open_UT)
Firsts_UT <- resize(granges(read1_UT),fix="start",1)
First_Pos_toCut_UT <- shift(granges(Firsts_UT[strand(read1_UT) == "+"]),
                            4)
First_Neg_toCut_UT <- shift(granges(Firsts_UT[strand(read1_UT) == "-"]),
                            -5)
Seconds_UT <- resize(granges(read2_UT),fix="start",1)
Second_Pos_toCut_UT <- shift(granges(Seconds_UT[strand(read2_UT) == "+"]),
                             4)
Second_Neg_toCut_UT <- shift(granges(Seconds_UT[strand(read2_UT) == "-"]),
                             -5)
test_toCut_UT <- c(First_Pos_toCut_UT,First_Neg_toCut_UT,
                   Second_Pos_toCut_UT,Second_Neg_toCut_UT)
cutsCoverage_UT <- coverage(test_toCut_UT)
# CTCF_Cuts_open_UT <- regionPlot(cutsCoverage_UT,
#                                 testRanges = myRes,
#                                 style = "point",
#                                 format="rlelist",distanceAround = 500)

save(read1_UT, read2_UT, Firsts_UT, First_Pos_toCut_UT, First_Neg_toCut_UT,
     Seconds_UT, Second_Pos_toCut_UT, Second_Neg_toCut_UT, test_toCut_UT,
     cutsCoverage_UT, myRes,
     file = "ATACseq_shiftinR_UT.RData")
