# memory is expressed in kilobytes

# alignment for most methods
STAR = 34681644
Salmon = 2363164 
STAR_FIL01 = 33568040

# alignment for rats
Salmon_boot = 3052528

# alignment for DEXSeq
DEXSeq_alignment = 770548
DEXSeq_alignment_FILT01 = 371236

# alignment for cjBitSeq
bowtie2 = 842648
bowtie2_FILT01 = 588668

####################################################################################
# DS methods:
####################################################################################
# BANDITS
BANDITS        = 1759124
BANDITS_FILT01 = 943424

# BayesDRIMSeq
BayesDRIMSeq        = 702656
BayesDRIMSeq_FILT01 = 658804

# cjBitSeq
cjBitSeq = max( 385432, 429932, 3424512)
cjBitSeq_FILT01 = max( 195696, 202204, 2175492)

# DEXSeq
DEXSeq        = max(10209716)
DEXSeq_FILT01 = max(5184604)

# DEXSeq ECCs
DEXSeq_ECCs        = max(1359544, 4672512)
DEXSeq_ECCs_FILT01 = max(1359544, 4171804)

# DEXSeq Transcripts:
DEXSeq_TECs        = 1981804
DEXSeq_TECs_FILT01 = 1491988

# DRIMSeq Transcripts:
DRIMSeq        = 749520
DRIMSeq_FILT01 = 672144

# limma 
limma        = 1403832
limma_FILT01 = 1242520

# rats
rats        = 5828580
rats_FILT01 = 3292916

####################################################################################
# Tables S6 and S7:
####################################################################################
individual_costs = data.frame(STAR = c(STAR, STAR_FIL01), 
                              Salmon, 
                              Salmon_boot, 
                              BANDITS = c( BANDITS, BANDITS_FILT01),
                              BayesDRIMSeq = c(BayesDRIMSeq, BayesDRIMSeq_FILT01),
                              DEXSeq_ECCs = c(DEXSeq_ECCs, DEXSeq_ECCs_FILT01),
                              DEXSeq_TECs = c(DEXSeq_TECs, DEXSeq_TECs_FILT01),
                              DRIMSeq = c(DRIMSeq, DRIMSeq_FILT01),
                              rats = c(rats, rats_FILT01),
                              bowtie2 = c(bowtie2, bowtie2_FILT01),
                              cjBitSeq = c(cjBitSeq, cjBitSeq_FILT01),
                              DEXSeq_alignment = c(DEXSeq_alignment, DEXSeq_alignment_FILT01),
                              DEXSeq = c(DEXSeq, DEXSeq_FILT01),
                              limma = c(limma, limma_FILT01) )

individual_costs = individual_costs[ , order(individual_costs[1,], decreasing = TRUE) ]

library(xtable)
xtable( t(individual_costs)/10^6, digits = 1 )
# Table S8

####################################################################################
# plot individual cost UNFILTERED:
####################################################################################
gg_data = data.frame(method = colnames(individual_costs),
                     GB = unlist(individual_costs[1,])/10^6, row.names = NULL)
gg_data$method = factor(gg_data$method, levels = gg_data$method[order(gg_data$GB, decreasing = TRUE)])
gg_data; gg_data$method

library(ggplot2)
ggp_3 = ggplot() +
  geom_bar(data = gg_data, aes_string(x = "method", y = "GB"), stat = "identity", fill = "blue") +
  theme_bw() + 
  xlab("") +
  ylab("GB") + 
  scale_y_sqrt( breaks = c(1, 2.5, 5, 7.5, 10, 15, 20, 30) ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(1.5)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
        aspect.ratio = 1 )
ggp_3

####################################################################################
# plot individual cost UNFILTERED:
####################################################################################
gg_data = data.frame(method = colnames(individual_costs),
                     GB = unlist(individual_costs[2,])/10^6, row.names = NULL)
gg_data$method = factor(gg_data$method, levels = gg_data$method[order(gg_data$GB, decreasing = TRUE)])
gg_data; gg_data$method

library(ggplot2)
ggp_4 = ggplot() +
  geom_bar(data = gg_data, aes_string(x = "method", y = "GB"), stat = "identity", fill = "blue") +
  theme_bw() + 
  xlab("") +
  ylab("GB") + 
  scale_y_sqrt( breaks = c(1, 2.5, 5, 7.5, 10, 15, 20, 30) ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(1.5)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
        aspect.ratio = 1 )
ggp_4

library(ggpubr)
ggarrange(ggp_3, ggp_4, labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom",
          nrow = 2, ncol = 1)
# saved with size: 8 * 12 cm
# name RAM_812, alias Figure S12.
