library(iCOBRA); library(data.table)

####################################################################################################################################
# I add a column of p.vals per method to the truth matrix.
name = "DTU_Sim_3vs3/simulation_details.txt" # for the 3 vs 3 simulation
name = "DTU_Sim_6vs6/simulation_details.txt" # for the 6 vs 6 simulation (filtered and unfiltered)

gene_to_transcript = fread(name, header = TRUE, sep="\t")

RES = unique(data.frame( gene_id = gene_to_transcript$gene_id, truth = gene_to_transcript$gene_ds_status ))

table_gene_id = table(gene_to_transcript$gene_id)
genes_atLeast2Transcripts = names( table_gene_id[table_gene_id > 1] )

# keep genes with at least 2 transcripts:
sel_atLeast2Transcripts = RES$gene_id %in% genes_atLeast2Transcripts
RES = RES[sel_atLeast2Transcripts ,]

# remove genes with < 20 counts:
load("Min_20_gene_counts.RData")
RES = RES[RES$gene_id %in% Min_20_gene_counts, ]

####################################################################################################################################
# BANDITS:
load("BANDITS results.RData")
# load("BANDITS results Filtered.RData") # for 6 vs 6 filtered transcriptome
library(BANDITS)
res = top_genes(results)
fdr = res$adj.p.values
fdr_inv = res$adj.p.values_inverted
gene_id = res$Gene_id

match = match(RES$gene_id, gene_id)

RES$BANDITS = fdr[match]
RES$BANDITS_inv = fdr_inv[match]

rm(fdr); rm(fdr_inv); rm(gene_id); rm(match); 

####################################################################################################################################
# BayesDRIMSeq
load("BayeDRIMSeq_results.RData")
fdr = 1-myRes$FDRraw
fdr_inv = 1-myRes$fdrTrust
gene_id = myRes$geneNames

# FDRRaw is the FDR
# fdrTrust is the corrected FDR for the inversion criterion.

match = match(RES$gene_id, gene_id)

RES$BayesDRIMSeq = fdr[match]
RES$BayesDRIMSeq_inv = fdr_inv[match]

rm(fdr); rm(fdr_inv); rm(gene_id); rm(match); 

####################################################################################################################################
# cjBitSeq
myRes = read.table("withinGeneEstimates.txt", header = T)
fdr = 1-myRes$FDRraw
fdr_inv = 1-myRes$FDR
gene_id = myRes$geneName

# FDRRaw is the FDR
# FDR is the corrected FDR for the inversion criterion.

match = match(RES$gene_id, gene_id)

RES$cjBitSeq = fdr[match]
RES$cjBitSeq_inv = fdr_inv[match]

rm(fdr); rm(fdr_inv); rm(gene_id); rm(match); 

####################################################################################################################################
# DEXSeq
load("DEXSeq results.RData")
fdr = res_gene$qval
gene_id = res_gene$gene

match = match(RES$gene_id, gene_id)

RES$DEXSeq = fdr[match]

rm(fdr); rm(gene_id); rm(match); 

####################################################################################################################################
# DEXSeq ECCs
load("DEXSeq_ECs results.RData")
fdr = res_gene$qval
gene_id = res_gene$gene

match = match(RES$gene_id, gene_id)

RES$DEXSeq_ECCs = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# DEXSeq transcripts
load("DEXSeq_TECs results.RData")
fdr = res_gene$qval
gene_id = res_gene$gene

match = match(RES$gene_id, gene_id)

RES$DEXSeq_TECs = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# DRIMSeq
load("DRIMSeq results.RData")
fdr = results_gene$adj_pvalue
gene_id = results_gene$gene_id

match = match(RES$gene_id, gene_id)

RES$DRIMSeq = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# limma
load("limma results.RData")
fdr = res$FDR
gene_id = res$GeneID

match = match(RES$gene_id, gene_id)

RES$limma = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# rats (with bootstrap)
load("rats results.RData")
fdr = res_boot$pval_corr
gene_id = res_boot$parent_id

match = match(RES$gene_id, gene_id)

RES$rats = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# SUPPA
load("SUPPA2 results.RData")
fdr = res_SUPPA
gene_id = names(res_SUPPA)

match = match(RES$gene_id, gene_id)

RES$SUPPA2 = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# I set to 1 the NA's, Inf and -1
RES[, -c(1,2)][ is.na(RES[,-c(1,2)]) ] = 1
RES[, -c(1,2)][ RES[,-c(1,2)] == Inf ] = 1
RES[, -c(1,2)][ RES[,-c(1,2)] == -1 ] = 1

# sort RES columns by name:
res = RES[, c(1,2, 2+order(colnames(RES)[-c(1,2)]) )]
RES = res

padj_icobra = data.frame(RES[,-c(1,2)])
rownames(padj_icobra) = RES$gene_id

truth_cobra = data.frame(status = RES$truth )
rownames(truth_cobra) = RES$gene_id

# icobra data:
cobra = COBRAData(padj = padj_icobra,
                  truth = truth_cobra,
                  object_to_extend = NULL)

cobraperf <- calculate_performance(cobra, binary_truth = "status")
slotNames(cobraperf)

## Define colour scheme:
library(RColorBrewer); #library(plotly)

colours = c(brewer.pal(3, "Greens")[2:3],
            brewer.pal(6, "BuPu")[3:6],
            brewer.pal(4, "Greys")[2:4],
            brewer.pal(6, "OrRd")[3:6],
            "white")

all_methods = c("BANDITS", "BANDITS_inv", "DRIMSeq", "SUPPA2", 
                "cjBitSeq", "cjBitSeq_inv", "BayesDRIMSeq", "BayesDRIMSeq_inv", 
                "DEXSeq", "DEXSeq_TECs", "DEXSeq_ECCs", "limma", "rats")
all_methods = sort(all_methods)
cols = match( c(colnames(RES)[-c(1,2)], "truth"),  c(all_methods, "truth") )

marker = list(color = c(colours[cols]))

## -------------------------------------------------------------------------------------------------
cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = marker$color, 
                                   facetted = TRUE, incloverall = FALSE,
                                   conditionalfill = FALSE)

# a function to scale the axes with any power transform.
pow_trans <- function(n = 2){
  scales::trans_new("pow",
                    transform = function(x) sign(x) * x^n,
                    inverse = function(x) sign(x) * abs(x)^(1/n))
}

library(ggplot2); library(wesanderson)

# TPR vs FDR plot, points only:
ggp_gene = plot_fdrtprcurve(cobraplot, linewidth = 1, pointsize = 5, plottype = "points") +
  scale_colour_manual(values=colours[cols], aesthetics = "colour", name  = "method") + 
  scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(trans = pow_trans(2)) +
  theme(legend.position="bottom", 
        legend.box = "vertical",
        legend.text=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(2)),
        legend.key.width=unit(2, "cm"),
        axis.text.x = element_text(angle = 0, hjust=0.5), 
        legend.title=element_text(size=rel(1.2)),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.75) +
  guides(linetype = guide_legend(ncol = 3, byrow = FALSE))

ggp_gene_3vs3 = ggp_gene      # for 3 vs 3 simulation;
ggp_gene_6vs6 = ggp_gene      # for 6 vs 6 simulation;
ggp_gene_6vs6_filt = ggp_gene # for 3 vs 3 simulation with transcript pre-filtering;

# TPR vs FDR plot, full curve:
ggp_gene_full = plot_fdrtprcurve(cobraplot, linewidth = 1, pointsize = 5, stripsize = 10) + 
  scale_colour_manual(values=colours[cols], aesthetics = "colour", name  = "method") + 
  scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(trans = pow_trans(2)) +
  theme(legend.position="bottom", 
        legend.box = "vertical",
        legend.text=element_text(size=rel(1.5)),
        axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(2)),
        legend.key.width=unit(3, "cm"),
        axis.text.x = element_text(angle = 0, hjust=0.5), 
        legend.title=element_text(size=rel(1.2)),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.5) +
  guides(linetype = guide_legend(ncol = 3, byrow = FALSE))

ggp_gene_full_3vs3 = ggp_gene_full      # for 3 vs 3 simulation;
ggp_gene_full_6vs6 = ggp_gene_full      # for 6 vs 6 simulation;
ggp_gene_full_6vs6_filt = ggp_gene_full # for 3 vs 3 simulation with transcript pre-filtering;

####################################################################################################################################
# repeat the above code for
# 1) 3 vs 3 simulation (ggp_gene_3vs3)
# 2) 6 vs 6 simulation; 3) 6 vs 6 simulation (filtering transcripts).
# Then merge plots into 1 Figure as follows:
library(ggpubr)

ggarrange(ggp_gene_3vs3, ggp_gene_6vs6, ggp_gene_6vs6_filt,
          labels = c("A", "B", "C"),
          common.legend = TRUE, legend = "bottom",
          nrow = 3, ncol = 1)
# saved with size: 8 * 20 cm
# name Gene3in1820LegendBottom, alias Figure 1.

ggarrange(ggp_gene_full_3vs3, ggp_gene_full_6vs6, ggp_gene_full_6vs6_filt,
          labels = c("A", "B", "C"),
          common.legend = TRUE, legend = "right",
          nrow = 3, ncol = 1)
# saved with size: 12 * 18 cm
# name GeneAllFull1218, alias Figure S1.
