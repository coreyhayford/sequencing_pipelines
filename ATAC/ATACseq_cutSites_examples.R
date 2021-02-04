## Manually curated TF list from genes in 5+ of top 10 GO terms
## from the unique idling molecular function GO type
.libPaths("~/R/rlib-3.6.0/")
library(MotifDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(soGGi)

setwd("/data/lola/hayforc/ATACseq/SKMEL5/preprocessed")
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

testMotif <- test[[1]]
myRes_test <- matchPWM(testMotif,BSgenome.Hsapiens.UCSC.hg38)

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

save(a,b, file = "ATACseq_cutSites_FOXO.RData")