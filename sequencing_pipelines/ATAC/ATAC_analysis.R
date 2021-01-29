setwd('/Volumes/Transcend/ATACseq/')

library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPseeker)
library(soGGi)
library(MotifDb)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicAlignments)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# UT_bam <- "/Volumes/Transcend/ATACseq/data/trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique_fullClean.bam"
# I_bam <- "/Volumes/Transcend/ATACseq/data/trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_dedup_unique_fullClean.bam"
# 
# UT_mappedReads <- idxstatsBam(UT_bam)
# UT_mappedReads_main <- UT_mappedReads[1:25,]
# I_mappedReads <- idxstatsBam(I_bam)
# I_mappedReads_main <- I_mappedReads[1:25,]
# 
# library(ggplot2)
# ggplot(UT_mappedReads_main,aes(seqnames,mapped,fill=seqnames))+
#   geom_bar(stat="identity")+coord_flip()+
#   ylab("Number of Mapped Reads") + xlab("Chromosome") +
#   ggtitle("Untreated Mapped Read Distribution by Chromosome")+theme_bw()+
#   theme(legend.text = element_text(size = 12), legend.position = "none", 
#         plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
#         legend.title = element_blank(), axis.title=element_text(size=12),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   ggsave("untreated_mappedReadDistribution.pdf")
# 
# ggplot(I_mappedReads_main,aes(seqnames,mapped,fill=seqnames))+
#   geom_bar(stat="identity")+coord_flip()+
#   ylab("Number of Mapped Reads") + xlab("Chromosome") +
#   ggtitle("Idling Mapped Read Distribution by Chromosome")+theme_bw()+
#   theme(legend.text = element_text(size = 12), legend.position = "none", 
#         plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=12),
#         legend.title = element_blank(), axis.title=element_text(size=12),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   ggsave("idling_mappedReadDistribution.pdf")


library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)

blkList <- import.bed("ENCFF356LFX.bed.gz")
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
          "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
          "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

openRegionPeaks_UT <- "data/untreated_peaks_all/untreated_all_peaks.narrowPeak_peaks.narrowPeak"
qcRes_UT <- ChIPQCsample("data/trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_dedup_unique_fullClean.bam",
                      peaks = openRegionPeaks_UT,
                      annotation ="hg38",
                      chromosomes = chrs,
                      blacklist = blkList,
                      verboseT = FALSE)

openRegionPeaks_I <- "data/idling_peaks_all/idling_all_peaks.narrowPeak_peaks.narrowPeak"
qcRes_I <- ChIPQCsample("data/trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_dedup_unique_fullClean.bam",
                         peaks = openRegionPeaks_I,
                         annotation ="hg38",
                         chromosomes = chrs,
                         blacklist = blkList,
                         verboseT = FALSE)

# save(qcRes_UT, qcRes_I, file = "allPeaks_narrow.RData")

load("allPeaks_narrow.RData")

library(reshape2)

qcMets <- data.frame(Untreated = c(qcRes_UT@SSD, qcRes_UT@SSDBL,
                                   qcRes_UT@CountsInPeaks, qcRes_UT@CountsInBlackList,
                                   qcRes_UT@FlagAndTagCounts[2], qcRes_UT@FlagAndTagCounts[6]),
                     Idling = c(qcRes_I@SSD, qcRes_I@SSDBL,
                                   qcRes_I@CountsInPeaks, qcRes_I@CountsInBlackList,
                                   qcRes_I@FlagAndTagCounts[2], qcRes_I@FlagAndTagCounts[6]))


rownames(qcMets) <- c("SSD", "SSDBL", "CountsPeaks", "CountsBL", "NumMapped", "NumDuped")
qcMets_T <- as.data.frame(t(qcMets))

qcMets_T$DupPerc <- (qcMets_T$NumDuped / qcMets_T$NumMapped) * 100 
qcMets_T$nonDupPerc <- 100 - qcMets_T$DupPerc
qcMets_T$RiP <- c(51.55, 34.35)
qcMets_T$NRiP <- 100 - qcMets_T$RiP
qcMets_T$RiBL <- c(0.002695422*100, 0.002269439*100)
qcMets_T$NRiBL <- 100 - qcMets_T$RiBL
qcMets_T$NPeaks <- c(76345, 96827)


qcMets_extras <- melt(qcMets_T[,7:ncol(qcMets_T)])
qcMets_extras$pop <- rep(c("Untreated", "Idling"), rep = 7)

qcMets_dup <- qcMets_extras[1:4,]
ggplot(qcMets_dup, aes(x = factor(pop, levels = c("Untreated", "Idling")),
                       y = value,
                       fill=variable)) +
  geom_bar(stat = "identity", color = "black") + theme_classic() +
  scale_fill_manual(values = c("grey20", "grey80"), 
                    labels = c("Duplicated", "Non-duplicated")) + 
  # scale_color_manual(values = c("red", "blue")) +
  ggtitle("Percentage of Reads Duplicated") +
  labs(x="Population", y="Percentage") + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", 
        plot.title = element_text(size = 14, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_blank(), axis.title=element_text(size=12)) +
  ggsave("PercDup.pdf", width = 5, height = 4)

qcMets_RiP <- qcMets_extras[5:8,]
ggplot(qcMets_RiP, aes(x = factor(pop, levels = c("Untreated", "Idling")),
                       y = value,
                       fill=factor(variable, levels = c("NRiP", "RiP")))) +
  geom_bar(stat = "identity", color = "black") + theme_classic() +
  scale_fill_manual(values = c("grey20", "grey80"),
                    labels = c("Outside", "Inside")) +
  # scale_color_manual(values = c("red", "blue")) +
  ggtitle("Percentage of Reads in Peaks") +
  labs(x="Population", y="Percentage") + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", 
        plot.title = element_text(size = 14, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_blank(), axis.title=element_text(size=12)) +
  ggsave("PercRiP.pdf", width = 5, height = 4)

qcMets_RiBL <- qcMets_extras[9:12,]
ggplot(qcMets_RiBL, aes(x = factor(pop, levels = c("Untreated", "Idling")),
                       y = value,
                       fill=variable)) +
  geom_bar(stat = "identity", color = "black") + theme_classic() +
  scale_fill_manual(values = c("grey20", "grey80"),
                    labels = c("Inside", "Outside")) +
  # scale_color_manual(values = c("red", "blue")) +
  ggtitle("Percentage of Reads in Blacklisted Regions") +
  labs(x="Population", y="Percentage") + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", 
        plot.title = element_text(size = 14, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_blank(), axis.title=element_text(size=12)) +
  ggsave("PercRiBL.pdf", width = 5, height = 4)

qcMets_numPeaks <- qcMets_extras[13:14,]
ggplot(qcMets_numPeaks, aes(x = factor(pop, levels = c("Untreated", "Idling")),
                            y = value,
                            fill=factor(pop, levels = c("Untreated", "Idling")))) + 
  geom_bar(stat = "identity", color = "black") + theme_classic() +
  scale_fill_manual(values = c("red", "blue")) +
  ggtitle("Number of Peaks") +
  labs(x="Population", y="Number of Peaks") + 
  theme(legend.text = element_text(size = 12), legend.position = "none", 
        plot.title = element_text(size = 14, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_blank(), axis.title=element_text(size=12)) +
  ggsave("NumPeaks.pdf", width = 5, height = 4)

# Plot some of these metrics
qcMets_PeakCount <- melt(melt(qcMets, id.vars = colnames(qcMets))["CountsPeaks",])
ggplot(data = qcMets_PeakCount, aes(x=variable, y=value, fill = variable)) +
  geom_bar(stat = "identity", color = "black") + theme_classic() +
  scale_fill_manual(values = c("red", "blue")) + ggtitle("Read Counts in Peaks") +
  labs(x="Population", y="Number of Reads") + theme(legend.position = "none") +
  ggsave("NumReads_Peaks.pdf", width = 5, height = 4)

qcMets_BLCount <- melt(melt(qcMets, id.vars = colnames(qcMets))["CountsBL",])
ggplot(data = qcMets_BLCount, aes(x=variable, y=value, fill = variable)) +
  geom_bar(stat = "identity", color = "black") + theme_classic() +
  scale_fill_manual(values = c("red", "blue")) + ggtitle("Read Counts in Blacklist Regions") +
  labs(x="Population", y="Number of Peaks") + theme(legend.position = "none") +
  ggsave("NumReads_Blacklist.pdf", width = 5, height = 4)

qcMets_MappedCount <- melt(melt(qcMets, id.vars = colnames(qcMets))["NumMapped",])
ggplot(data = qcMets_MappedCount, aes(x=variable, y=value, fill = variable)) +
  geom_bar(stat = "identity", color = "black") + theme_classic() +
  scale_fill_manual(values = c("red", "blue")) + ggtitle("Number of Mapped Reads") +
  labs(x="Population", y="Number of Reads") + theme(legend.position = "none") +
  ggsave("NumReads_Mapped.pdf", width = 5, height = 4)

qcMets_DuplicateCount <- melt(melt(qcMets, id.vars = colnames(qcMets))["NumDuped",])
ggplot(data = qcMets_DuplicateCount, aes(x=variable, y=value, fill = variable)) +
  geom_bar(stat = "identity", color = "black") + theme_classic() +
  scale_fill_manual(values = c("red", "blue")) + ggtitle("Number of Duplicated Reads") +
  labs(x="Population", y="Number of Reads") + theme(legend.position = "none") +
  ggsave("NumReads_Duplicated.pdf", width = 5, height = 4)
  
# Annotation

MacsCalls_UT <- granges(qcRes_UT)
data.frame(Blacklisted=sum(MacsCalls_UT %over% blkList),
           Not_Blacklisted=sum(!MacsCalls_UT %over% blkList))
MacsCalls_UT <- MacsCalls_UT[!MacsCalls_UT %over% blkList]

MacsCalls_I <- granges(qcRes_I)
data.frame(Blacklisted=sum(MacsCalls_I %over% blkList),
           Not_Blacklisted=sum(!MacsCalls_I %over% blkList))
MacsCalls_I <- MacsCalls_I[!MacsCalls_I %over% blkList]

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_UT <- getTagMatrix(MacsCalls_UT, windows=promoter)
tagMatrix_I <- getTagMatrix(MacsCalls_I, windows=promoter)

plotAvgProf(tagMatrix_UT, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix_I, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

MacsCalls_UT_Anno <-  annotatePeak(MacsCalls_UT, 
                                   TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
MacsCalls_I_Anno <-  annotatePeak(MacsCalls_I, 
                                  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)

# MacsCalls_UT_Anno <-  annotatePeak(MacsCalls_UT, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#                                    tssRegion=c(-1000, 1000))
# MacsCalls_I_Anno <-  annotatePeak(MacsCalls_I, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#                                   tssRegion=c(-1000, 1000))

plotAnnoPie(MacsCalls_UT_Anno)
plotAnnoPie(MacsCalls_I_Anno)

vennpie(MacsCalls_UT_Anno)
vennpie(MacsCalls_I_Anno)

upsetplot(MacsCalls_UT_Anno)
upsetplot(MacsCalls_I_Anno)


# Limma venn diagram
peaks <- dir("data/DifferentialAnnotation", pattern = "*.narrowPeak", 
             full.names = TRUE)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)

names(myPeaks) <- c("Idling", "Untreated")
Group <- factor(c("Idling", "Untreated"))

# library(devtools)
# install_github("https://github.com/ColeWunderlich/soGGi.git")

source("soGGi_runConsensusRegions_fixed.R")
consensusToCount <- runConsensusRegions(GRangesList(myPeaks), "none")
consensusToCount

library(limma)

as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(Untreated, Idling) %>% 
  vennDiagram(main = "Overlap for ATAC Open regions")

vd <- as.data.frame(elementMetadata(consensusToCount))[,c("Untreated", "Idling")]
vd_counts <- plyr::count(vd)

library(eulerr)
library(UpSetR)

shared_peaks <- c("Untreated" = vd_counts$freq[2],
                  "Idling" = vd_counts$freq[1],
                  "Untreated&Idling" = vd_counts$freq[3])

venn_sharedPeaks <- euler(shared_peaks)
plot(venn_sharedPeaks, fills = c("red", "blue"),
     shape = "ellipse", quantities = TRUE)

upset(fromExpression(shared_peaks), order.by = "freq",
      sets.bar.color = c("blue", "red"), 
      mb.ratio = c(0.55, 0.45), nintersects = NA,
      show.numbers = "no", point.size = 5, line.size = 0.5, mainbar.y.max = 120000,
      mainbar.y.label = "Peak Intersections", sets.x.label = "Set Size", 
      text.scale = c(2, 1.5, 2, 1.5, 2, 1))

# Annotate unique regions
library(clusterProfiler)
library(ChIPseeker)

CTC_UT <- consensusToCount[(elementMetadata(consensusToCount)[,"Untreated"] == 1) &
                           (elementMetadata(consensusToCount)[,"Idling"] == 0)]
CTC_I <- consensusToCount[(elementMetadata(consensusToCount)[,"Untreated"] == 0) &
                             (elementMetadata(consensusToCount)[,"Idling"] == 1)]
CTC_con <- consensusToCount[(elementMetadata(consensusToCount)[,"Untreated"] == 1) &
                            (elementMetadata(consensusToCount)[,"Idling"] == 1)]


anno_UT <- annotatePeak(CTC_UT, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_I <- annotatePeak(CTC_I, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
anno_con <- annotatePeak(CTC_con, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)

save(anno_UT, anno_I, anno_con, 
     file = "SKMEL5_ATAC_annotatedPeaks_uniqueShared.RData")

plotAnnoBar(anno_UT)
plotAnnoBar(anno_I)
plotAnnoBar(anno_con)

library(ggimage)
upsetplot(anno_UT, vennpie = T) + ggsave("anno_uniqueUT.pdf")
upsetplot(anno_I, vennpie = T) + ggsave("anno_uniqueI.pdf")
upsetplot(anno_con, vennpie = T) + ggsave("anno_uniqueCon.pdf")

# Unique - submit to great, GO, cluster profiler
library(org.Hs.eg.db)
go_UT_BP <- enrichGO(as.data.frame(as.GRanges(anno_UT))$geneId, OrgDb = "org.Hs.eg.db", 
                ont = "BP", maxGSSize = 5000)
go_I_BP <- enrichGO(as.data.frame(as.GRanges(anno_I))$geneId, OrgDb = "org.Hs.eg.db", 
                ont = "BP", maxGSSize = 5000)
go_con_BP <- enrichGO(as.data.frame(as.GRanges(anno_con))$geneId, OrgDb = "org.Hs.eg.db", 
                ont = "BP", maxGSSize = 5000)

go_UT_MF <- enrichGO(as.data.frame(as.GRanges(anno_UT))$geneId, OrgDb = "org.Hs.eg.db", 
                     ont = "MF", maxGSSize = 5000)
go_I_MF <- enrichGO(as.data.frame(as.GRanges(anno_I))$geneId, OrgDb = "org.Hs.eg.db", 
                    ont = "MF", maxGSSize = 5000)
go_con_MF <- enrichGO(as.data.frame(as.GRanges(anno_con))$geneId, OrgDb = "org.Hs.eg.db", 
                      ont = "MF", maxGSSize = 5000)

go_UT_CC <- enrichGO(as.data.frame(as.GRanges(anno_UT))$geneId, OrgDb = "org.Hs.eg.db", 
                     ont = "CC", maxGSSize = 5000)
go_I_CC <- enrichGO(as.data.frame(as.GRanges(anno_I))$geneId, OrgDb = "org.Hs.eg.db", 
                    ont = "CC", maxGSSize = 5000)
go_con_CC <- enrichGO(as.data.frame(as.GRanges(anno_con))$geneId, OrgDb = "org.Hs.eg.db", 
                      ont = "CC", maxGSSize = 5000)

# save(go_UT_BP, go_I_BP, go_con_BP,
#      go_UT_MF, go_I_MF, go_con_MF,
#      go_UT_CC, go_I_CC, go_con_CC,
#      file = "GO_enrichment_ATACunique.RData")

load("GO_enrichment_ATACunique.RData")

dotplot(go_UT_BP) + ggsave("GOenrichment_UT_BP.pdf", width = 8, height = 5)
dotplot(go_I_BP) + ggsave("GOenrichment_I_BP.pdf", width = 8, height = 5)
dotplot(go_con_BP) + ggsave("GOenrichment_con_BP.pdf", width = 8, height = 5)

dotplot(go_UT_MF) + ggsave("GOenrichment_UT_MF.pdf", width = 8, height = 5)
dotplot(go_I_MF) + ggsave("GOenrichment_I_MF.svg", width = 8, height = 5)
dotplot(go_con_MF) + ggsave("GOenrichment_con_MF.pdf", width = 8, height = 5)

dotplot(go_UT_CC) + ggsave("GOenrichment_UT_CC.pdf", width = 8, height = 5)
dotplot(go_I_CC) + ggsave("GOenrichment_I_CC.pdf", width = 8, height = 5)
dotplot(go_con_CC) + ggsave("GOenrichment_con_CC.pdf", width = 8, height = 5)

df <- data.frame(UT_MF = go_UT_MF@result$Description[1:25],
                 I_MF = go_I_MF@result$Description[1:25],
                 con_MF = go_con_MF@result$Description[1:25])

files <- list(CTC_UT, CTC_I, CTC_con)
names(files) <- c("Untreated", "Idling", "Shared")
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList) +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 14, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("ATAC_annotationDistribution_UniqueShared.pdf", width = 8, height = 6)

plotDistToTSS(peakAnnoList) +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 14, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("ATAC_distanceToTSS_UniqueShared.pdf", width = 8, height = 6)

tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000)) +
  theme(legend.text = element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = 14, hjust = 0.5), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggsave("ATAC_averageBindingProfile_UniqueShared.pdf", width = 6, height = 4)

# plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
# tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compEgo <- compareCluster(geneCluster   = genes,
                           fun           = "enrichGO",
                           OrgDb='org.Hs.eg.db',
                          ont = "MF",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compEgo, showCategory = 30, title = "Pathway Enrichment Analysis")

untreatedAnno <- data.frame(peakAnnoList[["Untreated"]]@anno)
idlingAnno <- data.frame(peakAnnoList[["Idling"]]@anno)

# Get the entrez IDs
entrez_UT_MF <- untreatedAnno$geneId
entrez_I_MF <- idlingAnno$geneId

# Return the gene symbol for the set of Entrez IDs
library(AnnotationDbi)
annotations_edb_UT_MF <- AnnotationDbi::select(org.Hs.eg.db,
                                         keys = entrez_UT_MF,
                                         columns = c("GENENAME"),
                                         keytype = "ENTREZID")
annotations_edb_I_MF <- AnnotationDbi::select(org.Hs.eg.db,
                                               keys = entrez_I_MF,
                                               columns = c("GENENAME"),
                                               keytype = "ENTREZID")
# Change IDs to character type to merge
annotations_edb_UT_MF$ENTREZID <- as.character(annotations_edb_UT_MF$ENTREZID)
annotations_edb_I_MF$ENTREZID <- as.character(annotations_edb_I_MF$ENTREZID)

# Write to file
untreatedAnno_gene <- untreatedAnno %>% 
  left_join(annotations_edb_UT_MF, by=c("geneId"="ENTREZID")) 

idlingAnno_gene <- idlingAnno %>% 
  left_join(annotations_edb_I_MF, by=c("geneId"="ENTREZID")) 


testVect <- go_I_MF[1:10,]$geneID
data.frame(as.list(testVect))
as.data.frame(strsplit(as.character(as.list(testVect[[1]])), split = "/"))
Vect <- lapply(testVect, function(x) strsplit(as.character(x), split = "/"))
Vect <- lapply(Vect, function(y) unlist(y))

I_MF_genes <- list()
z <- 0
for (i in seq_along(Vect)) {
  comp <- unlist(Vect[i])
  temp <- AnnotationDbi::select(org.Hs.eg.db,
                        keys = comp,
                        columns = c("SYMBOL"),
                        keytype = "ENTREZID")
  I_MF_genes[[i]] <- temp$SYMBOL
  a = length(comp)
  z = z + a
}


sort(table(unlist(I_MF_genes)))

## Manually curated TF list from genes in 5+ of top 10 GO terms
## from the unique idling molecular function GO type
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(soGGi)
load("ATACseq_shiftinR_UT.RData")
load("ATACseq_shiftinR_I.RData")


## The main TFs include:
# FOXO_
test <- query(MotifDb, c("FOXO", "hsapiens", "Jaspar2018"))[2:4]
# CREB
test <- query(MotifDb, c("CREB", "hsapiens", "Jaspar2018"))
# AP1
test <- query(MotifDb, c("AP1", "hsapiens", "JASPAR_CORE"))
# MZF-1
test <- query(MotifDb, c("MZF1", "hsapiens", "Jaspar2018"))
# NRF1
test <- query(MotifDb, c("NRF", "hsapiens", "Jaspar2018"))
# NFKB
test <- query(MotifDb, c("NFKB", "hsapiens", "Jaspar2018"))
# STAT_
test <- query(MotifDb, c("STAT", "hsapiens", "Jaspar2018"))
# JUN
test <- query(MotifDb, c("JUN", "hsapiens", "Jaspar2018"))[c(5:10, 14:34)]
# ATF2
test <- query(MotifDb, c("ATF2", "hsapiens"))[3:4]
# GATA_
test <- query(MotifDb, c("GATA", "hsapiens", "Jaspar2018"))
# PAX_
test <- query(MotifDb, c("PAX", "hsapiens", "Jaspar2018"))
# PPARgamma
test <- query(MotifDb, c("PPARG", "hsapiens", "Jaspar2018"))
# XBP1
test <- query(MotifDb, c("XBP1", "hsapiens", "Jaspar2018"))
# NFAT
test <- query(MotifDb, c("NFAT", "hsapiens", "Jaspar2018"))[2]
# NFSF
# N/A
# NF1
# N/A

testMotif <- test[[1]]
myRes_test <- matchPWM(testMotif,BSgenome.Hsapiens.UCSC.hg38)

# CREB <- query(MotifDb, c("CREB"))
# crebMotif <- CREB[[1]]
# myRes_CREB <- matchPWM(crebMotif,BSgenome.Hsapiens.UCSC.hg38)
# setwd('/Volumes/Transcend/ATACseq/data/')



# library(ggplot2)
# 
# dir.create("CREB_cutSites")
# setwd("CREB_cutSites/")
# dir.create("Untreated")
# dir.create("Idling")
# setwd('/Volumes/Transcend/ATACseq/data/')

a<-regionPlot(cutsCoverage_UT,
              testRanges = myRes_test,
              style = "point",
              format="rlelist",
              distanceAround = 500,
              samplename = "Untreated")
b<-regionPlot(cutsCoverage_I,
              testRanges = myRes_CREB,
              style = "point",
              format="rlelist",
              distanceAround = 500,
              samplename = "Idling")
# plotRegion(a,outliers = 0.001)+
#   ggtitle("Nucleosome Free Cuts Centred on CREB")+theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5))+
#   ggsave("CREB_I_allChromosomes", width = 6, height = 5)



mainChrom <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
               "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
               "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

for (i in mainChrom){
  a<-regionPlot(cutsCoverage_I,
                testRanges = GRanges(i,ranges(myRes_CREB)),
                style = "point",
                format="rlelist",
                distanceAround = 500,
                samplename = paste("Idling", i))
  plotRegion(a,outliers = 0.001)+
    ggtitle("NucFree Cuts Centred on CREB")+theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    ggsave(paste0("CREB_cutSites/Idling/CREB_I_", i,".pdf"), width = 6, height = 5)
}
  
  
# vennplot(genes)
# ChIPQCreport(qcRes_UT, reportName="ChIP QC report: Untreated All", reportFolder="ChIPQCreport")
# ChIPQCreport(qcRes_I, reportName="ChIP QC report: Idling", reportFolder="ChIPQCreport")


# myMetrics_UT <- QCmetrics(qcRes_UT)
# myMetrics_UT[c("RiBL%","RiP%")]
# 
# myMetrics_UT <- QCmetrics(qcRes_UT)
# myMetrics_UT[c("RiBL%","RiP%")]

# flgCounts_UT <- flagtagcounts(qcRes_UT)
# DupRate_UT <- flgCounts_UT["DuplicateByChIPQC"]/flgCounts_UT["Mapped"]
# DupRate_UT*100

# flgCounts_I <- flagtagcounts(qcRes_I)
# DupRate_I <- flgCounts_I["DuplicateByChIPQC"]/flgCounts_I["Mapped"]
# DupRate_I*100


MacsCalls_UT <- granges(qcRes_UT)
data.frame(Blacklisted=sum(MacsCalls_UT %over% blkList),
           Not_Blacklisted=sum(!MacsCalls_UT %over% blkList))
MacsCalls_UT <- MacsCalls_UT[!MacsCalls_UT %over% blkList]

MacsCalls_I <- granges(qcRes_I)
data.frame(Blacklisted=sum(MacsCalls_I %over% blkList),
           Not_Blacklisted=sum(!MacsCalls_I %over% blkList))
MacsCalls_I <- MacsCalls_I[!MacsCalls_I %over% blkList]

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix_UT <- getTagMatrix(MacsCalls_UT, windows=promoter)
tagMatrix_I <- getTagMatrix(MacsCalls_I, windows=promoter)

plotAvgProf(tagMatrix_UT, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix_I, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

MacsCalls_UT_Anno <-  annotatePeak(MacsCalls_UT, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                   tssRegion=c(-1000, 1000))
MacsCalls_I_Anno <-  annotatePeak(MacsCalls_I, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                                  tssRegion=c(-1000, 1000))

# plotAnnoPie(MacsCalls_UT_Anno)
# plotAnnoPie(MacsCalls_I_Anno)
# 
# vennpie(MacsCalls_UT_Anno)
# vennpie(MacsCalls_I_Anno)

plotAnnoBar(MacsCalls_UT_Anno)
plotAnnoBar(MacsCalls_I_Anno)

library(ggupset)
upsetplot(MacsCalls_UT_Anno)
upsetplot(MacsCalls_I_Anno)

plotDistToTSS(MacsCalls_UT_Anno, title="Distribution of transcription factor-binding loci relative to TSS")
plotDistToTSS(MacsCalls_I_Anno, title="Distribution of transcription factor-binding loci relative to TSS")



library(clusterProfiler)
library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(MacsCalls_UT_Anno)$geneId)
head(pathway1, 5)$Description

pathway2 <- enrichPathway(as.data.frame(MacsCalls_I_Anno)$geneId)
head(pathway2, 5)$Description

dotplot(pathway1)
dotplot(pathway2)

gene_UT <- seq2gene(MacsCalls_UT, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway_TSS_UT <- enrichPathway(gene_UT)
head(pathway_TSS_UT, 5)$Description

gene_I <- seq2gene(MacsCalls_I, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway_TSS_I <- enrichPathway(gene_I)
head(pathway_TSS_I, 5)$Description

dotplot(pathway_TSS_UT)
dotplot(pathway_TSS_I)

MacsGR_Anno_UT <- as.GRanges(MacsCalls_UT_Anno)
MacsGR_TSS_UT <-   MacsGR_Anno_UT[abs(MacsGR_Anno_UT$distanceToTSS) < 500]
MacsGR_TSS_UT[1,]

MacsGR_Anno_I <- as.GRanges(MacsCalls_I_Anno)
MacsGR_TSS_I <-   MacsGR_Anno_I[abs(MacsGR_Anno_I$distanceToTSS) < 500]
MacsGR_TSS_I[1,]

library(rGREAT)
great_Job_UT <- submitGreatJob(MacsCalls_UT, species = "hg38")
great_Job_I <- submitGreatJob(MacsCalls_I, species = "hg38")

great_ResultTable_UT = getEnrichmentTables(great_Job_UT, category = "GO")
# names(great_ResultTable_UT)
great_ResultTable_UT[["GO Biological Process"]][1:20, ]$name

great_ResultTable_I = getEnrichmentTables(great_Job_I, category = "GO")
great_ResultTable_I[["GO Biological Process"]][1:20, ]$name

library(eulerr)
library(UpSetR)

shared_peaks <- c("Untreated" = 70384,
                  "Idling" = 116023,
                  "Untreated&Idling" = 62477)

venn_sharedPeaks <- euler(shared_peaks)
plot(venn_sharedPeaks, fills = c("red", "blue"),
     shape = "ellipse", quantities = TRUE)

upset(fromExpression(shared_peaks), order.by = "freq",
      sets.bar.color = c("blue", "red"), 
      mb.ratio = c(0.55, 0.45), nintersects = NA,
      show.numbers = "no", point.size = 5, line.size = 0.5, mainbar.y.max = 120000,
      mainbar.y.label = "Peak Intersections", sets.x.label = "Set Size", 
      text.scale = c(2, 1.5, 2, 1.5, 2, 1))
