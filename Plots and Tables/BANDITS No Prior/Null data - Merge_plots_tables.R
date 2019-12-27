####################################################################################################################################
# Load plots and tables
####################################################################################################################################
load("Gene_RES_ggplot_points.RData")
load("Gene_FDR_RES_ggplot_points.RData")
load("Tr_RES_ggplot_points.RData")
load("Tr_FDR_RES_ggplot_points.RData")

####################################################################################################################################
# Merge Tables
####################################################################################################################################
methods = unique(c(as.character(Gene_FDR_RES_ggplot_points$method), as.character(Tr_FDR_RES_ggplot_points$method)))
methods
methods = sort(methods, decreasing = FALSE)
methods

res = matrix(NA, ncol = 4, nrow = length(methods))
rownames(res) = methods

sel = Gene_RES_ggplot_points$threshold == 0.05
res[,1] = Gene_RES_ggplot_points[sel,1][ match(methods, Gene_RES_ggplot_points[sel,3]) ]

sel = Gene_FDR_RES_ggplot_points$threshold == 0.05
res[,2] = Gene_FDR_RES_ggplot_points[sel,1][ match(methods, Gene_FDR_RES_ggplot_points[sel,3]) ]

sel = Tr_RES_ggplot_points$threshold == 0.05
res[,3] = Tr_RES_ggplot_points[sel,1][ match(methods, Tr_RES_ggplot_points[sel,3]) ]

sel = Tr_FDR_RES_ggplot_points$threshold == 0.05
res[,4] = Tr_FDR_RES_ggplot_points[sel,1][ match(methods, Tr_FDR_RES_ggplot_points[sel,3]) ]

library(xtable)
xtable(100 * res, digits = 2)
# Table S12

####################################################################################################################################
# Merge Plots
####################################################################################################################################
# Make plot with GENE (unfilt and filt):
ggp_1  = ggp_1 +
  theme(legend.text=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
        axis.title=element_text(size=rel(1.3)),
        legend.key.width=unit(2, "cm"),
        legend.title=element_text(size=rel(1)) ) +
  guides(colour = guide_legend(order = 1, ncol = 3, byrow = FALSE),
         linetype = guide_legend(order = 1))

ggp_2  = ggp_2 +
  theme(legend.text=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
        axis.title=element_text(size=rel(1.3)),
        legend.key.width=unit(2, "cm"),
        legend.title=element_text(size=rel(1)) ) +
  guides(colour = guide_legend(order = 1, ncol = 3, byrow = FALSE),
         linetype = guide_legend(order = 1))

library(ggpubr)
ggarrange(ggp_1, ggp_2, labels = c("A", "B"),
          common.legend = FALSE, legend = "bottom",
          nrow = 2, ncol = 1)
# saved with size: 8 * 12 cm
# name NoPrior_FPs2in1Gene812, alias Figure S21

####################################################################################################################################
# Make plot with TRANSCRIPT (unfilt and filt):
ggp_3  = ggp_3 +
  theme(legend.text=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
        axis.title=element_text(size=rel(1.3)),
        legend.key.width=unit(2, "cm"),
        legend.title=element_text(size=rel(1)) ) +
  guides(colour = guide_legend(order = 1, ncol = 3, byrow = FALSE),
         linetype = guide_legend(order = 1))

ggp_4  = ggp_4 +
  theme(legend.text=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
        axis.title=element_text(size=rel(1.3)),
        legend.key.width=unit(2, "cm"),
        legend.title=element_text(size=rel(1)) ) +
  guides(colour = guide_legend(order = 1, ncol = 3, byrow = FALSE),
         linetype = guide_legend(order = 1))

ggarrange(ggp_3, ggp_4, labels = c("A", "B"),
          common.legend = FALSE, legend = "bottom",
          nrow = 2, ncol = 1)
# saved with size: 8 * 12 cm
# name NoPrior_FPs2in1Tr812, alias Figure S22
