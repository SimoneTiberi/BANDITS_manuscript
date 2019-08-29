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
library(BANDITS)
res = top_genes(results)
fdr = res$adj.p.values
fdr_inv = res$adj.p.values_inverted
gene_id = res$Gene_id

match = match(RES$gene_id, gene_id)

RES$BANDITS = fdr[match]
RES$BANDITS_inv = fdr_inv[match]

rm(fdr); rm(fdr_inv); rm(gene_id); rm(match); rm(results)

####################################################################################################################################
# BANDITS with transcript pre-filtering:
load("BANDITS results Filtered.RData") # for 6 vs 6 filtered transcriptome
library(BANDITS)
res = top_genes(results)
fdr = res$adj.p.values
fdr_inv = res$adj.p.values_inverted
gene_id = res$Gene_id

match = match(RES$gene_id, gene_id)

RES$BANDITS_Filtered = fdr[match]
RES$BANDITS_inv_Filtered = fdr_inv[match]

rm(fdr); rm(fdr_inv); rm(gene_id); rm(match); rm(results)

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

# DRIMSeq and My p.vals here:
cobra = COBRAData(padj = padj_icobra,
                  truth = truth_cobra,
                  object_to_extend = NULL)

####################################################################################################################################
# icobra plots:
## -------------------------------------------------------------------------------------------------
cobraperf <- calculate_performance(cobra, binary_truth = "status")
slotNames(cobraperf)

library(RColorBrewer); 
colours = c(brewer.pal(3, "BuPu")[2:3],
            brewer.pal(3, "OrRd")[2:3], "white")
marker = list(color = colours)

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

library(ggplot2);
library(wesanderson)

# points only:
ggp_gene_BANDITS = plot_fdrtprcurve(cobraplot, linewidth = 1, pointsize = 5, plottype = "points",
                                    xaxisrange = c(0, 0.1),
                                    yaxisrange = c(0.82, 0.87)) +
  theme( strip.background = element_blank(),
         strip.text = element_blank(),
         legend.text=element_text(size=rel(1.5)),
         axis.text=element_text(size=rel(1.3)),
         axis.title=element_text(size=rel(2)),
         legend.key.width=unit(3, "cm"),
         legend.title=element_text(size=rel(1.2)),
         aspect.ratio = 0.75,
         axis.text.x = element_text(angle = 0, hjust=0.5) ) +
  scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1))

ggp_gene_BANDITS_3vs3 = ggp_gene_BANDITS      # for 3 vs 3 simulation;
ggp_gene_BANDITS_6vs6 = ggp_gene_BANDITS      # for 6 vs 6 simulation;

####################################################################################################################################
# repeat the above code for
# 1) 3 vs 3 simulation (ggp_gene_BANDITS_3vs3)
# 2) 6 vs 6 simulation; (ggp_gene_BANDITS_6vs6)
# Then merge plots into 1 Figure as follows:

# GENE PLOT:
ggarrange(ggp_gene_BANDITS_3vs3, ggp_gene_BANDITS_6vs6,
          labels = c("A", "B"),
          common.legend = TRUE, legend = "right",
          nrow = 2, ncol = 1)
# saved with size: 12 * 18 cm
# name 2in1BANDITSGene1218, alias Figure S3.
