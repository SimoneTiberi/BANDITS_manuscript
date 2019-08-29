# times are expressed in seconds.

# alignment for most methods
STAR = 16244
Salmon = 2493
STAR_FIL01 = 2376 + 12883

# alignment for rats
Salmon_boot = 8897

# alignment for DEXSeq
DEXSeq_alignment = 112248
DEXSeq_alignment_FILT01 = 104397

# alignment for cjBitSeq
bowtie2 = 1627 + 65055
bowtie2_FILT01 = 652 + 47993

####################################################################################
# DS methods:
####################################################################################
# BANDITS
BANDITS        = 10469.629
BANDITS_FILT01 = 3526.125

1/60^2 * (STAR + Salmon + BANDITS)
1/60^2 * (STAR + Salmon + BANDITS_FILT01)

# BayesDRIMSeq
BayesDRIMSeq        = 1737.231
BayesDRIMSeq_FILT01 = 1081.561

1/60^2 * (STAR + Salmon + BayesDRIMSeq)
1/60^2 * (STAR + Salmon + BayesDRIMSeq_FILT01)

# cjBitSeq
cjBitSeq = 38872 + 860 + 254010
cjBitSeq_FILT01 = 28300 + 147 + 175592

1/60^2 * (bowtie2 + cjBitSeq)
1/60^2 * (STAR + Salmon + bowtie2_FILT01 + cjBitSeq_FILT01)

# DEXSeq
DEXSeq        = 3661.774
DEXSeq_FILT01 = 1823.948

1/60^2 * (STAR + DEXSeq_alignment + DEXSeq)
1/60^2 * (STAR + Salmon + STAR_FIL01 + DEXSeq_alignment_FILT01 + DEXSeq_FILT01)

# DEXSeq ECCs
DEXSeq_ECCs        = 60 + 2509.126
DEXSeq_ECCs_FILT01 = 60 + 2294.913

1/60^2 * (STAR + Salmon + DEXSeq_ECCs)
1/60^2 * (STAR + Salmon + DEXSeq_ECCs_FILT01)

# DEXSeq Transcripts:
DEXSeq_TECs        = 737.685
DEXSeq_TECs_FILT01 = 203.982

1/60^2 * (STAR + Salmon + DEXSeq_TECs)
1/60^2 * (STAR + Salmon + DEXSeq_TECs_FILT01)

# DRIMSeq Transcripts:
DRIMSeq        = 799.250
DRIMSeq_FILT01 = 574.778

1/60^2 * (STAR + Salmon + DRIMSeq)
1/60^2 * (STAR + Salmon + DRIMSeq_FILT01)

# limma 
limma        = 94.234
limma_FILT01 = 57.025

1/60^2 * (STAR + DEXSeq_alignment + limma)
1/60^2 * (STAR + Salmon + STAR_FIL01 + DEXSeq_alignment_FILT01 + limma_FILT01)

# rats
rats        = 1762.992
rats_FILT01 = 1612.762

1/60^2 * (STAR + Salmon_boot + rats)
1/60^2 * (STAR + Salmon_boot + rats_FILT01)

####################################################################################
# Tables S6 and S7:
####################################################################################
individual_costs = data.frame(STAR, 
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

overall_time = data.frame(BANDITS = c( STAR + Salmon + BANDITS, STAR + Salmon + BANDITS_FILT01),
                          BayesDRIMSeq = c(STAR + Salmon + BayesDRIMSeq, STAR + Salmon + BayesDRIMSeq_FILT01),
                          cjBitSeq = c(bowtie2 + cjBitSeq, STAR + Salmon + bowtie2_FILT01 + cjBitSeq_FILT01),
                          DEXSeq = c(STAR + DEXSeq_alignment + DEXSeq, STAR + Salmon + STAR_FIL01 + DEXSeq_alignment_FILT01 + DEXSeq_FILT01),
                          DEXSeq_ECCs = c(STAR + Salmon + DEXSeq_ECCs, STAR + Salmon + DEXSeq_ECCs_FILT01),
                          DEXSeq_TECs = c(STAR + Salmon + DEXSeq_TECs, STAR + Salmon + DEXSeq_TECs_FILT01),
                          DRIMSeq = c(STAR + Salmon + DRIMSeq, STAR + Salmon + DRIMSeq_FILT01),
                          limma = c(STAR + DEXSeq_alignment + limma, STAR + Salmon + STAR_FIL01 + DEXSeq_alignment_FILT01 + limma),
                          rats = c(STAR + Salmon_boot + rats, STAR + Salmon_boot + rats_FILT01))

#individual_costs = individual_costs[ , order(individual_costs[2,], decreasing = TRUE) ]
overall_time = overall_time[ , order(overall_time[2,], decreasing = TRUE) ]

library(xtable)
xtable( t(individual_costs)/60, digits = 0 )
# Table S6

xtable( t(overall_time)/60, digits = 0 )
# Table S7

####################################################################################
# plot overall cost:
####################################################################################

methods = c("BANDITS", "BayesDRIMSeq", "cjBitSeq", "DEXSeq", "DEXSeq_ECCs", "DEXSeq_TECs", "DRIMSeq", "limma", "rats")

overall_time = c(STAR + Salmon + BANDITS,
                 STAR + Salmon + BayesDRIMSeq,
                 bowtie2 + cjBitSeq,
                 STAR + DEXSeq_alignment + DEXSeq,
                 STAR + Salmon + DEXSeq_ECCs,
                 STAR + Salmon + DEXSeq_TECs,
                 STAR + Salmon + DRIMSeq,
                 STAR + DEXSeq_alignment + limma,
                 STAR + Salmon_boot + rats )

alignment_time = c(STAR + Salmon,
                   STAR + Salmon,
                   bowtie2,
                   STAR,
                   STAR + Salmon,
                   STAR + Salmon,
                   STAR + Salmon,
                   STAR,
                   STAR + Salmon_boot )

gg_data = data.frame(method = methods, minutes = overall_time/60, alignment = alignment_time/60)
gg_data$method = factor(gg_data$method, levels = gg_data$method[order(gg_data$minutes, decreasing = TRUE)])
gg_data; gg_data$method

library(ggplot2)
ggp_1 = ggplot() +
  geom_bar(data = gg_data, aes_string(x = "method", y = "minutes"), stat = "identity", fill = "red") +
  geom_bar(data = gg_data, aes_string(x = "method", y = "alignment"), stat = "identity", fill = "blue") +
  theme_bw() + 
  xlab("") +
  ylab("minutes") + 
  scale_y_sqrt( breaks = c(60, 300, 600, 1200, 2400, 3600, 4800) ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(1.5)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
        aspect.ratio = 1 )
ggp_1  

####################################################################################
# plot cost of DS method only (after STAR + Salmon)
####################################################################################
methods_sel = c("BANDITS", "BayesDRIMSeq", "DEXSeq_ECCs", "DEXSeq_TECs", "DRIMSeq")
gg_data_method_only = data.frame(method = methods[methods %in% methods_sel], minutes = (overall_time - alignment_time)[methods %in% methods_sel]/60)
gg_data_method_only$method = factor(gg_data_method_only$method, levels = gg_data_method_only$method[order(gg_data_method_only$minutes, decreasing = TRUE)])
gg_data_method_only; gg_data_method_only$method

library(ggplot2)
ggp_2 = ggplot() +
  geom_bar(data = gg_data_method_only, aes_string(x = "method", y = "minutes"), 
           stat = "identity", fill = "red") +
  theme_bw() + 
  xlab("") +
  ylab("minutes") + 
  scale_y_sqrt( breaks = c(5, 15, 30, 60, 120, 180, 500, 1000, 2000, 3000, 4000) ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(1.5)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
        aspect.ratio = 1 )
ggp_2  



####################################################################################
# plot overall cost FILT01
####################################################################################

methods = c("BANDITS", "BayesDRIMSeq", "cjBitSeq", "DEXSeq", "DEXSeq_ECCs", "DEXSeq_TECs", "DRIMSeq", "limma", "rats")

overall_time = c(STAR + Salmon + BANDITS_FILT01,
                 STAR + Salmon + BayesDRIMSeq_FILT01,
                 STAR + Salmon + bowtie2_FILT01 + cjBitSeq_FILT01,
                 STAR + Salmon + STAR_FIL01 + DEXSeq_alignment_FILT01 + DEXSeq_FILT01,
                 STAR + Salmon + DEXSeq_ECCs_FILT01,
                 STAR + Salmon + DEXSeq_TECs_FILT01,
                 STAR + Salmon + DRIMSeq_FILT01,
                 STAR + Salmon + STAR_FIL01 + DEXSeq_alignment_FILT01 + limma,
                 STAR + Salmon_boot + rats_FILT01 )

alignment_time = c(STAR + Salmon,
                   STAR + Salmon,
                   STAR + Salmon + bowtie2_FILT01,
                   STAR + Salmon + STAR_FIL01,
                   STAR + Salmon,
                   STAR + Salmon,
                   STAR + Salmon,
                   STAR + Salmon + STAR_FIL01,
                   STAR + Salmon_boot )

gg_data = data.frame(method = methods, minutes = overall_time/60, alignment = alignment_time/60)
gg_data$method = factor(gg_data$method, levels = gg_data$method[order(gg_data$minutes, decreasing = TRUE)])
gg_data; gg_data$method

library(ggplot2)
ggp_3 = ggplot() +
  geom_bar(data = gg_data, aes_string(x = "method", y = "minutes"), stat = "identity", fill = "red") +
  geom_bar(data = gg_data, aes_string(x = "method", y = "alignment"), stat = "identity", fill = "blue") +
  theme_bw() + 
  xlab("") +
  ylab("minutes") + 
  scale_y_sqrt( breaks = c(60, 300, 600, 1200, 2400, 3600, 4800) ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(1.5)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
        aspect.ratio = 1 )
ggp_3  
#   scale_y_sqrt()+

####################################################################################
# plot cost of DS method only (after STAR + Salmon)
####################################################################################
methods_sel = c("BANDITS", "BayesDRIMSeq", "DEXSeq_ECCs", "DEXSeq_TECs", "DRIMSeq")
gg_data_method_only = data.frame(method = methods[methods %in% methods_sel], minutes = (overall_time - alignment_time)[methods %in% methods_sel]/60)
gg_data_method_only$method = factor(gg_data_method_only$method, levels = gg_data_method_only$method[order(gg_data_method_only$minutes, decreasing = TRUE)])
gg_data_method_only; gg_data_method_only$method


library(ggplot2)
ggp_4 = ggplot() +
  geom_bar(data = gg_data_method_only, aes_string(x = "method", y = "minutes"), 
           stat = "identity", fill = "red") +
  theme_bw() + 
  xlab("") +
  ylab("minutes") + 
  scale_y_sqrt( breaks = c(5, 15, 30, 60, 120, 180, 500, 1000, 2000, 3000, 4000) ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(1.5)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, color = "grey", linetype = 2),
        aspect.ratio = 1 )
ggp_4  

####################################################################################
# put the 4 plots in 1 image
####################################################################################
library(ggpubr)
ggarrange(ggp_1, ggp_2,
          #ggp_1 + rremove("x.text"), ggp_2 + rremove("x.text"), 
          ggp_3, ggp_4,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
# saved with size: 8 * 8 cm
# name BenchmarkTimesALL, alias Figure 5.
