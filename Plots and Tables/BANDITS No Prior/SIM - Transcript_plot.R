library(iCOBRA); library(data.table)

####################################################################################################################################
# I add a column of p.vals per method to the truth matrix.
name = "DTU_Sim_3vs3/simulation_details.txt" # for the 3 vs 3 simulation
name = "DTU_Sim_6vs6/simulation_details.txt" # for the 6 vs 6 simulation (filtered and unfiltered)

gene_to_transcript = fread(name, header = TRUE, sep="\t")

RES = data.frame( gene_id = gene_to_transcript$gene_id,
                  transcript_id = gene_to_transcript$transcript_id,
                  truth = gene_to_transcript$transcript_ds_status )

table_gene_id = table(gene_to_transcript$gene_id)
genes_atLeast2Transcripts = names( table_gene_id[table_gene_id > 1] )

# keep genes with at least 2 transcripts:
sel_atLeast2Transcripts =  RES$gene_id %in% genes_atLeast2Transcripts
RES = RES[sel_atLeast2Transcripts ,]

# remove genes with < 20 counts:
load("Min_20_gene_counts.RData")
RES = RES[RES$gene_id %in% Min_20_gene_counts, ]

# remove transcripts with < 10 counts:
load("Min_10_transcript_counts.RData")
RES = RES[RES$transcript_id %in% Min_10_transcript_counts, ]

####################################################################################################################################
# BANDITS:
load("BANDITS results.RData")
# load("BANDITS results Filtered.RData") # for 6 vs 6 filtered transcriptome
library(BANDITS)
res = top_transcripts(results)
fdr = res$adj.p.values
fdr_maxGene = res$Max_Gene_Tr.Adj.p.val
tr_id = res$Transcript_id

match = match(RES$transcript_id, tr_id)

RES$BANDITS = fdr[match]
RES$BANDITS_maxGene = fdr_maxGene[match]

rm(fdr); rm(fdr_maxGene); rm(tr_id); rm(match); 

####################################################################################################################################
# BANDITS No prior:
load("BANDITS results NoPrior.RData")
# load("BANDITS results Filtered NoPrior.RData") # for 6 vs 6 filtered transcriptome
res = top_transcripts(results)
fdr = res$adj.p.values
fdr_maxGene = res$Max_Gene_Tr.Adj.p.val
tr_id = res$Transcript_id

match = match(RES$transcript_id, tr_id)
RES$transcript_id[1:5]
tr_id[match][1:5]

RES$BANDITS_NoPrior = fdr[match]
RES$BANDITS_NoPrior_maxGene = fdr_maxGene[match]

rm(res); rm(fdr); rm(tr_id); rm(match); rm(results)

####################################################################################################################################
# I set to 1 the NA's, Inf and -1
RES[, -c(1,2,3)][ is.na(RES[,-c(1,2,3)]) ] = 1
RES[, -c(1,2,3)][ RES[,-c(1,2,3)] == Inf] = 1
RES[, -c(1,2,3)][ RES[,-c(1,2,3)] == -1 ] = 1

# sort RES columns by name:
res <- RES[, c(1,2, 3, 3 +order(colnames(RES)[-c(1,2,3)]) )]
RES = res

padj_icobra = data.frame(RES[,-c(1,2,3)])
rownames(padj_icobra) = RES$transcript_id

truth_cobra = data.frame(status = RES$truth )
rownames(truth_cobra) = RES$transcript_id

# icobra data:
cobra = COBRAData(padj = padj_icobra,
                  truth = truth_cobra,
                  object_to_extend = NULL)
cobraperf <- calculate_performance(cobra, binary_truth = "status")
slotNames(cobraperf)

## Define colour scheme:
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

library(ggplot2); library(wesanderson)

# TPR vs FDR plot, points only:
ggp_tr = plot_fdrtprcurve(cobraplot, linewidth = 1, pointsize = 5, plottype = "points") +
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

ggp_tr_3vs3 = ggp_tr      # for 3 vs 3 simulation;
ggp_tr_6vs6 = ggp_tr      # for 6 vs 6 simulation;
ggp_tr_6vs6_filt = ggp_tr # for 3 vs 3 simulation with transcript pre-filtering;

####################################################################################################################################
# repeat the above code for
# 1) 3 vs 3 simulation (ggp_gene_3vs3)
# 2) 6 vs 6 simulation; 3) 6 vs 6 simulation (filtering transcripts).
# Then merge plots into 1 Figure as follows:
library(ggpubr)

ggarrange(ggp_tr_3vs3, ggp_tr_6vs6, ggp_tr_6vs6_filt,
          labels = c("A", "B", "C"),
          common.legend = TRUE, legend = "bottom",
          nrow = 3, ncol = 1)
# saved with size: 12 * 18 cm
# name NoPrior_Tr3in1218LegendBottom, alias Figure S19.
