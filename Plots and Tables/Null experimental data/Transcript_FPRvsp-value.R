library(iCOBRA); library(data.table)

####################################################################################################################################
# I add a column of p.vals per method to the truth matrix.
load("gene_transcript_gtf.RData")

gene_to_transcript = X

# create the truth table:
RES = data.frame( gene_id = as.character(gene_to_transcript$gene_id),
                  transcript_id = as.character(gene_to_transcript$transcript_id),
                  truth = 0 )
# truth = 0 for all genes/transcripts (this is a null data analysis)

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
res = top_transcripts(results)
fdr = res$p.values
fdr_maxGene = res$Max_Gene_Tr.p.val
tr_id = res$Transcript_id

match = match(RES$transcript_id, tr_id)

RES$BANDITS = fdr[match]
RES$BANDITS_maxGene = fdr_maxGene[match]

rm(fdr); rm(fdr_maxGene); rm(gene_id); rm(match); 

####################################################################################################################################
# cjBitSeq
myRes = read.table("transcriptLevelEstimates.txt", header = T)
fdr = 1-myRes$ProbDE
tr_id = myRes$trName

match = match(RES$transcript_id, tr_id)

RES$cjBitSeq = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# DEXSeq transcripts
load("DEXSeq_TECs results.RData")
fdr = res$pvalue
tr_id = res$featureID

match = match(RES$transcript_id, tr_id)

RES$DEXSeq_TECs = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# DRIMSeq
load("DRIMSeq results.RData")
fdr = results_trancript$pvalue
tr_id = results_trancript$feature_id

match = match(RES$transcript_id, tr_id)

RES$DRIMSeq = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# rats (with bootstrap)
load("rats results.RData")
fdr = res_boot_transcripts$pval
tr_id = res_boot_transcripts$target_id

match = match(RES$transcript_id, tr_id)

RES$rats = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# SUPPA
res_SUPPA = read.table("DTU_OUTPUT_rawPval.dpsi")
tr_id = substring(rownames(res_SUPPA), 17, 31)
fdr = res_SUPPA[, 2]

match = match(RES$transcript_id, tr_id)

RES$SUPPA2 = fdr[match]

rm(fdr); rm(gene_id); rm(match)

####################################################################################################################################
# icobra data
# I set to 1 the NA's
RES[, -c(1,2,3)][ is.na(RES[,-c(1,2,3)]) ] = 1
RES[, -c(1,2,3)][ RES[,-c(1,2,3)] == Inf] = 1
RES[, -c(1,2,3)][ RES[,-c(1,2,3)] == -1 ] = 1

res <- RES[, c(1,2, 3, 3 +order(colnames(RES)[-c(1,2,3)]) )]
RES = res

grid = seq(0,1, 0.01)
fps = sapply(4:ncol(RES), function(i){
  sapply(grid, function(x){
    mean(RES[,i] < x )
  })
})
length(fps); length(rep(grid, ncol(RES) - 3)); length(rep(colnames(RES)[-c(1:3)], each = length(grid)) )

RES_ggplot = data.frame( FPs = c(fps), threshold = rep(grid, ncol(RES) - 3),
                         method = rep(colnames(RES)[-c(1:3)], each = length(grid)) )
head(RES_ggplot); tail(RES_ggplot)

# make the data.frame for the 1, 5 and 10 % points.
RES_ggplot_points = RES_ggplot[RES_ggplot$threshold %in% c(0.01, 0.05, 0.1),]

# define colour scale for plots:

## Define colour scheme:
library(RColorBrewer); #library(plotly)

colours = c(brewer.pal(3, "Greens")[2:3],
            brewer.pal(6, "BuPu")[3:6],
            brewer.pal(4, "Greys")[2:4],
            brewer.pal(6, "OrRd")[3:6],
            "white")

all_methods = c("BANDITS", "BANDITS_maxGene", "DRIMSeq", "SUPPA2", 
                "cjBitSeq", "cjBitSeq_inv", "BayesDRIMSeq", "BayesDRIMSeq_inv", 
                "DEXSeq", "DEXSeq_TECs", "DEXSeq_ECCs", "limma", "rats")
all_methods = sort(all_methods)
cols = match( c(colnames(RES)[-c(1,2,3)]), c(all_methods) )

# make a ggplot of it:
library(ggplot2); library(wesanderson); library(iCOBRA)

ggp_3 = ggplot(RES_ggplot, aes(threshold, FPs, colour = method)) +
  scale_shape_discrete(name = "FDR") +
  labs(x = "p-value") + 
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
ggp_3

# save plot and table to make joint plot (Supplementary Figures and Tables)
Tr_RES_ggplot_points = RES_ggplot_points
save(Tr_RES_ggplot_points, ggp_3, file = "Tr_RES_ggplot_points.RData")
