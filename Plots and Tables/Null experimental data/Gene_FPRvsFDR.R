library(iCOBRA); library(data.table)

####################################################################################################################################
# I add a column of p.vals per method to the truth matrix.
load("null_data/gene_transcript_gtf.RData")
gene_to_transcript = data.frame(gene_id = X$gene_id, transcript_id = X$transcript_id, stringsAsFactors = FALSE)

table_gene_id = table(gene_to_transcript$gene_id)
genes_atLeast2Transcripts = names( table_gene_id[table_gene_id > 1] )
length(genes_atLeast2Transcripts)

# create the truth table:
genes = unique(genes_atLeast2Transcripts)
RES = data.frame( gene_id = genes, truth = 0 )
# truth = 0 for all genes/transcripts (this is a null data analysis)

# keep genes with at least 2 transcripts:
sel_atLeast2Transcripts =  RES$gene_id %in% genes_atLeast2Transcripts
RES = RES[sel_atLeast2Transcripts ,]

# remove genes with < 20 counts:
load("Min_20_gene_counts.RData")
RES = RES[RES$gene_id %in% Min_20_gene_counts, ]

####################################################################################################################################
# BANDITS:
load("BANDITS results.RData")
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
# icobra data
# I set to 1 the NA's and Inf
RES[, -c(1,2)][ is.na(RES[,-c(1,2)]) ] = 1
RES[, -c(1,2)][ RES[,-c(1,2)] == Inf ] = 1
RES[, -c(1,2)][ RES[,-c(1,2)] == -1 ] = 1

res <- RES[, c(1,2, 2+order(colnames(RES)[-c(1,2)]) )]
RES = res

# compute the FPs for the grid:
grid = seq(0,1, 0.01)
fps = sapply(3:ncol(RES), function(i){
  sapply(grid, function(x){
    mean(RES[,i] < x )
  })
})
length(fps); length(rep(grid, ncol(RES) - 2)); length(rep(colnames(RES)[-c(1,2)], each = length(grid)) )

RES_ggplot = data.frame( FPs = c(fps), threshold = rep(grid, ncol(RES) - 2),
                         method = rep(colnames(RES)[-c(1,2)], each = length(grid)) )
head(RES_ggplot); tail(RES_ggplot)

# make the data.frame for the 1, 5 and 10 % points.
RES_ggplot_points = RES_ggplot[RES_ggplot$threshold %in% c(0.01, 0.05, 0.1),]

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
cols = match( c(colnames(RES)[-c(1,2)]), c(all_methods) )

# make a ggplot of it:
library(ggplot2); library(wesanderson); library(iCOBRA)

ggp_2 = ggplot(RES_ggplot, aes(threshold, FPs, colour = method)) +
  scale_shape_discrete(name = "FDR") +
  labs(x = "FDR") + 
  labs(y = "FPR") + 
  scale_colour_manual(values=colours[cols], aesthetics = "colour") +
  geom_line(size = 1.5) + # size (line width) must be regulated via geom_line
  scale_x_sqrt( breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1) ) +
  scale_y_sqrt(breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1) ) +
  geom_vline(xintercept=0.01, linetype=2) +
  geom_vline(xintercept=0.05, linetype=2) +
  geom_vline(xintercept=0.1, linetype=2) +
  guides(colour = guide_legend(order = 1, ncol = 4, byrow = FALSE) ) +
  theme_bw() + 
  theme(legend.position="bottom", 
        legend.box = "vertical",
        legend.text=element_text(size=rel(1.2)),
        axis.text=element_text(size=rel(1.3)),
        axis.title=element_text(size=rel(1.8)),
        legend.key.width=unit(2, "cm"),
        axis.text.x = element_text(angle = 0), 
        legend.title=element_text(size=rel(1.2)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.2, color = "grey", linetype = 2),
        aspect.ratio = 0.75)
ggp_2

# save plot and table to make joint plot (Supplementary Figures and Tables)
Gene_FDR_RES_ggplot_points = RES_ggplot_points
save(Gene_FDR_RES_ggplot_points, ggp_2, file = "Gene_FDR_RES_ggplot_points.RData")

####################################################################################################################################
library(ggplot2); library(wesanderson); library(iCOBRA)
ggp_2 + 
  guides(colour = guide_legend(order = 1, ncol = 3, byrow = FALSE),
         linetype = guide_legend(order = 1))
# saved with size: 8 * 8 cm
# name FPsGeneFDR, alias Figure 4.
