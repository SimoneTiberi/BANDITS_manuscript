library(iCOBRA); library(data.table)

####################################################################################################################################
# I add a column of p.vals per method to the truth matrix.
name = "DTU_Sim_6vs6/simulation_details.txt"

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

# stratify for gene abundance:
load("tot_counts_per_gene.RData")
tot_counts_per_gene = tot_counts_per_gene[names(tot_counts_per_gene) %in% RES$gene_id]

splits = quantile(tot_counts_per_gene, probs = c(1/3,2/3))
splits

low_genes = names( tot_counts_per_gene )[ tot_counts_per_gene < splits[1]  ]
high_genes = names( tot_counts_per_gene )[ tot_counts_per_gene >= splits[2]  ]
mid_genes = names( tot_counts_per_gene )[ (tot_counts_per_gene >= splits[1]) &  (tot_counts_per_gene < splits[2])  ]

length(low_genes); 
length(mid_genes); 
length(high_genes)
[1] 4471
[1] 4471
[1] 4472

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
# I set to 1 the NA's, Inf and -1
RES[, -c(1,2)][ is.na(RES[,-c(1,2)]) ] = 1
RES[, -c(1,2)][ RES[,-c(1,2)] == Inf ] = 1
RES[, -c(1,2)][ RES[,-c(1,2)] == -1 ] = 1

# sort RES columns by name:
res = RES[, c(1,2, 2+order(colnames(RES)[-c(1,2)]) )]
RES = res

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

pow_trans <- function(n = 2){
  scales::trans_new("pow",
                    transform = function(x) sign(x) * x^n,
                    inverse = function(x) sign(x) * abs(x)^(1/n))
}

library(ggplot2)
# a function to scale the axes with any power transform.

# make 3 plots, 1 for low_genes, 1 for mid_genes and 1 for high_genes
gg_fdr = list()
for(expression in 1:3){
  keep_genes = list(low_genes, mid_genes,
                    high_genes)[[expression]]
  keep_rows = RES$gene_id %in% keep_genes
  
  padj_icobra = data.frame(RES[keep_rows,-c(1,2)])
  rownames(padj_icobra) = RES$gene_id[keep_rows]
  
  truth_cobra = data.frame(status = RES$truth[keep_rows] )
  rownames(truth_cobra) = RES$gene_id[keep_rows]
  
  # DRIMSeq and My p.vals here:
  cobra = COBRAData(padj = padj_icobra,
                    truth = truth_cobra,
                    object_to_extend = NULL)
  
  cobraperf <- calculate_performance(cobra, binary_truth = "status")
  slotNames(cobraperf)
  
  
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = marker$color, 
                                     facetted = TRUE, incloverall = FALSE,
                                     conditionalfill = FALSE)
  
  gg_fdr[[expression]] =  plot_fdrtprcurve(cobraplot, linewidth = 1, pointsize = 5, plottype = "points") +
    scale_colour_manual(values=colours[cols], aesthetics = "colour", name  = "method") +
    theme( strip.background = element_blank(),
           strip.text = element_blank(),
           legend.position="bottom", 
           legend.box = "vertical",
           legend.text=element_text(size=rel(1.5)),
           axis.text=element_text(size=rel(1.3)),
           axis.title=element_text(size=rel(2)),
           legend.key.width=unit(3, "cm"),
           axis.text.x = element_text(angle = 0, hjust=0.5),
           legend.title=element_text(size=rel(1.2)),
           aspect.ratio = 0.5 ) +
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(trans = pow_trans(2))
}

# Make plot with all 3 together:
library(ggpubr)
ggarrange(gg_fdr[[1]], gg_fdr[[2]], gg_fdr[[3]],
          labels = c("A", "B", "C"),
          common.legend = TRUE, legend = "right",
          nrow = 3, ncol = 1)
# saved with size: 12 * 18 cm
# name FDRGene6vs63in1Vertical1218, alias Figure S5.
