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
# BANDITS No prior:
load("BANDITS results NoPrior.RData")
res = top_genes(results)
fdr = res$adj.p.values
fdr_inv = res$adj.p.values_inverted
gene_id = res$Gene_id

match = match(RES$gene_id, gene_id)
RES$gene_id[1:5]
gene_id[match][1:5]

RES$BANDITS_NoPrior = fdr[match]
RES$BANDITS_NoPrior_inv = fdr_inv[match]

rm(res); rm(fdr); rm(gene_id); rm(match); rm(results)

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

library(RColorBrewer); 
colours = c(brewer.pal(3, "BuPu")[2:3],
            brewer.pal(3, "OrRd")[2:3])

# make a ggplot of it:
library(ggplot2); library(wesanderson); library(iCOBRA)

ggp_2 = ggplot(RES_ggplot, aes(threshold, FPs, colour = method)) +
  scale_shape_discrete(name = "FDR") +
  labs(x = "FDR") + 
  labs(y = "FPR") + 
  scale_colour_manual(values=colours, aesthetics = "colour") +
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
