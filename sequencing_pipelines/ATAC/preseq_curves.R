setwd('/Volumes/Transcend/ATACseq/preseq_results/')

library(tidyverse)

UT_lc_res <- read_tsv('trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_lcExtrapResults.txt') %>%
  mutate(Library="Untreated")
I_lc_res <- read_tsv('trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_lcExtrapResults.txt') %>%
  mutate(Library="Idling")

lc_res <- bind_rows(UT_lc_res, I_lc_res)

ggplot(lc_res) + theme_bw() +
  geom_ribbon(data = UT_lc_res, aes(x=TOTAL_READS/10^6, 
                                    ymin=LOWER_0.95CI/10^6, 
                                    ymax=UPPER_0.95CI/10^6), 
              fill="red", alpha = 0.4) +
  geom_ribbon(data = I_lc_res, aes(x=TOTAL_READS/10^6, 
                                   ymin=LOWER_0.95CI/10^6, 
                                   ymax=UPPER_0.95CI/10^6), 
              fill="blue", alpha = 0.4) +
  geom_line(aes(x=TOTAL_READS/10^6, 
                y=EXPECTED_DISTINCT/10^6, 
                color=Library), 
            size=1, linetype="dashed") + 
  geom_vline(xintercept = 100, colour="black", linetype="dotted") +
  coord_cartesian(x=c(0,600), y=c(0,200)) +
  xlab("Sequenced reads (M)") + 
  ylab("Distinct reads (M)") +
  ggtitle("Estimated Complexity of ATAC-seq Libraries") +
  # scale_x_continuous(breaks = seq(0, 750, by = 50)) +
  # scale_y_continuous(breaks = seq(0, 75, by = 5)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("blue", "red")) +
  ggsave("SKMEL5_ATACseq_LCcurve_Untreated-Idling_comparison.pdf",
         width = 6, height = 4)

#####


c_res_UT <- read_tsv('trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_cCurveResults.txt') %>%
  mutate(Library="Untreated")
c_res_I <- read_tsv('trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_cCurveResults.txt') %>%
  mutate(Library="Idling")

c_res <- bind_rows(c_res_UT, c_res_I)

ggplot(c_res) + theme_bw() +
  geom_line(aes(x=total_reads/10^6, 
                y=distinct_reads/10^6, 
                color=Library), size=1, linetype="solid") + 
  # coord_cartesian(x=c(0,75), y=c(0,37.5)) +
  xlab("Sequenced reads (M)") + 
  ylab("Distinct reads (M)") +
  ggtitle("Observed Complexity of ATAC-seq Libraries") +
  # scale_x_continuous(breaks = seq(0, 75, by = 5)) +
  # scale_y_continuous(breaks = seq(0, 37.5, by = 2.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("blue", "red")) +
  ggsave("SKMEL5_ATACseq_Ccurve_Untreated-Idling_comparison.pdf",
         width = 6, height = 4)


####

c_res_UT <- read_tsv('trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_cCurveResults.txt') %>%
  mutate(Type="Observed")
lc_res_UT <- read_tsv('trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_lcExtrapResults.txt') %>%
  select(TOTAL_READS, EXPECTED_DISTINCT) %>%
  mutate(Type="Expected")

colnames(lc_res_UT) <- c("total_reads", "distinct_reads", "Type")
res_UT <- bind_rows(lc_res_UT, c_res_UT)

ggplot(res_UT, aes(x=total_reads/10^6, y=distinct_reads/10^6)) + 
  theme_bw() +
  geom_line(aes(linetype=Type, colour=Type), size=2) + 
  coord_cartesian(x=c(0,200), y=c(0,100)) +
  xlab("Sequenced reads (M)") + 
  ylab("Distinct reads (M)") +
  ggtitle("Observed and Expected Complexity of Untreated") +
  # scale_x_continuous(breaks = seq(0, 75, by = 5)) +
  # scale_y_continuous(breaks = seq(0, 37.5, by = 2.5)) +
  # theme_minimal(base_size = 12, base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("grey", "red")) +
  ggsave("SKMEL5_ATACseq_LC-C-curveComparison_Untreated.pdf",
         width = 6, height = 4)

c_res_I <- read_tsv('trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_cCurveResults.txt') %>%
  mutate(Type="Observed")
lc_res_I <- read_tsv('trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_lcExtrapResults.txt') %>%
  select(TOTAL_READS, EXPECTED_DISTINCT) %>%
  mutate(Type="Expected")

colnames(lc_res_I) <- c("total_reads", "distinct_reads", "Type")
res_I <- bind_rows(lc_res_I, c_res_I)

ggplot(res_I, aes(x=total_reads/10^6, y=distinct_reads/10^6)) + 
  theme_bw() +
  geom_line(aes(linetype=Type, colour=Type), size=2) + 
  coord_cartesian(x=c(0,200), y=c(0,100)) +
  xlab("Sequenced reads (M)") + 
  ylab("Distinct reads (M)") +
  ggtitle("Observed and Expected Complexity of Idling") +
  # scale_x_continuous(breaks = seq(0, 75, by = 5)) +
  # scale_y_continuous(breaks = seq(0, 37.5, by = 2.5)) +
  # theme_minimal(base_size = 12, base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("grey", "blue")) +
  ggsave("SKMEL5_ATACseq_LC-C-curveComparison_Idling.pdf",
         width = 6, height = 4)

####
# ISM

insert_UT <- read_tsv('trimmed_3334-CH-1-TAGGCATG-ACTGCATA_S121_aligned_sorted_ISM.txt', 
                        skip = 10) %>%
  mutate(Library="Untreated")
insert_UT <- mutate(insert_UT, cpm=(All_Reads.fr_count/sum(insert_UT$All_Reads.fr_count))*1e6)

insert_I <- read_tsv('trimmed_3334-CH-2-CGTACTAG-CTCTCTAT_S122_aligned_sorted_ISM.txt', 
                     skip = 10) %>%
  mutate(Library="Idling")
insert_I <- mutate(insert_I, cpm=(All_Reads.fr_count/sum(insert_I$All_Reads.fr_count))*1e6)

insert_all <- bind_rows(insert_UT, insert_I)
insert_all$Library <- factor(insert_all$Library, 
                             levels = c("Untreated", "Idling"))

ggplot(insert_all, aes(x=insert_size, y=cpm, color=Library)) +
  geom_line(size = 0.7) +
  facet_grid(.~Library) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(x=c(0,800)) +
  xlab("Insert Size (bp)") + 
  ylab("Counts per Million") +
  ggtitle("Insert Size Distribution of ATAC-seq Libraries", subtitle = "0-800bp") +
  scale_x_continuous(breaks = seq(0, 900, by = 150)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggsave("SKMEL5_ATACseq_Untreated-Iding-Comparison_notDedup_800bp_ISM.pdf",
         width = 6, height = 4)

ggplot(insert_all, aes(x=insert_size, y=cpm, color = Library)) +
  geom_line(size = 0.7) +
  facet_grid(.~Library) +
  scale_color_manual(values = c("red", "blue")) +
  coord_cartesian(x=c(0,200)) +
  xlab("Insert Size (bp)") + 
  ylab("Counts per Million") +
  ggtitle("Insert Size Distribution of ATAC-seq Libraries", subtitle = "0-200bp") +
  # scale_x_continuous(breaks = seq(0, 900, by = 30)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ggsave("SKMEL5_ATACseq_Untreated-Iding-Comparison_notDedup_200bp_ISM.pdf",
         width = 6, height = 4)
