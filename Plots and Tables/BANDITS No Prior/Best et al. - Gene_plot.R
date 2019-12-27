library(iCOBRA); library(data.table)

####################################################################################################################################
# I add a column of p.vals per method to the truth matrix.
load("Best_data/gene_transcript_gtf.RData")
gene_to_transcript = data.frame(gene_id = X$gene_id, transcript_id = X$transcript_id, stringsAsFactors = FALSE)

table_gene_id = table(gene_to_transcript$gene_id)
genes_atLeast2Transcripts = names( table_gene_id[table_gene_id > 1] )
length(genes_atLeast2Transcripts)

# load validated genes names:
genes = unique(genes_atLeast2Transcripts)
load("validated_genes.RData")
truth = ifelse(genes %in% validated, 1, 0)
sum(truth)
# 82
RES = data.frame( gene_id = genes, truth = truth )

# keep genes with at least 2 transcripts:
sel_atLeast2Transcripts =  RES$gene_id %in% genes_atLeast2Transcripts
RES = RES[sel_atLeast2Transcripts ,]

# load genes rakings:
load("gene_rank.RData")
sum(validated %in% gene_rank[1:10000])
# position of validated genes in ranking of top expressed genes (10,000 needed to include them all):
which(gene_rank%in% validated)

# filter results to only keep the top 10,000 expressed genes:
RES = RES[RES$gene_id %in% c(gene_rank[1:10000], validated), ]

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
# I set to 1 the NA's, Inf and -1
RES[, -c(1,2)][ is.na(RES[,-c(1,2)]) ] = 1
RES[, -c(1,2)][ RES[,-c(1,2)] == Inf ] = 1
RES[, -c(1,2,3)][ RES[,-c(1,2,3)] == -1 ] = 1

# sort RES columns by name:
res = RES[, c(1,2, 2+order(colnames(RES)[-c(1,2)]) )]
RES = res

padj_icobra = data.frame(RES[,-c(1,2)])
rownames(padj_icobra) = RES$gene_id

truth_cobra = data.frame(status = RES$truth )
rownames(truth_cobra) = RES$gene_id

## Define colour scheme:
library(RColorBrewer); 
colours = c(brewer.pal(3, "BuPu")[2:3],
            brewer.pal(3, "OrRd")[2:3], "white")
marker = list(color = colours)

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = marker$color, 
                                   facetted = TRUE, incloverall = FALSE,
                                   conditionalfill = FALSE)

library(ggplot2); library(wesanderson)

plot_roc(cobraplot, linewidth = 1) +
  scale_colour_manual(values = marker$color, name  = "method") +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 0.75,
        legend.position="bottom",
        legend.text=element_text(size=rel(1.2)),
        legend.key.width=unit(2, "cm"),
        legend.title=element_text(size=rel(1.2)),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.2, color = "grey", linetype = 2 ),
        axis.text.x = element_text(angle = 0, hjust=0.5) ) +
  scale_x_sqrt(  breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1) ) +
  guides(linetype = guide_legend(ncol = 3, byrow = FALSE))

# saved with size: 8 * 8 cm
# name NoPrior_ROC_Best_xsqrt, alias Figure S20.

####################################################################################################################################
# AUC and Median ranking:
####################################################################################################################################
library(ROCR)
par(mfrow = c(1,1))
preds = apply(cobra@padj, 2, function(u) 
  prediction(-u, cobra@truth))

# Study the (average) position of the 16 validated genes in the ranking of ds genes:
ranking = apply(RES[,-c(1:2)], 2, rank)

####################################################################################################################################
# Top 100-200 ranked genes:
####################################################################################################################################
# identified genes in the top 100 by each method:
RES_ranking_01 = apply(RES[,-c(1,2)], 2, function(x){
  ifelse( validated %in% RES$gene_id[ order(x)[1:100] ], 1, 0)
})
rownames(RES_ranking_01) = validated
top_100 = colSums(RES_ranking_01)
top_100

# identified genes in the top 100 by each method:
RES_ranking_01 = apply(RES[,-c(1,2)], 2, function(x){
  ifelse( validated %in% RES$gene_id[ order(x)[1:200] ], 1, 0)
})
rownames(RES_ranking_01) = validated
top_200 = colSums(RES_ranking_01)
top_200

####################################################################################################################################
# Gene ontology analysis:
####################################################################################################################################
# need to convert ensambl ids into Entrez ids:
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=RES$gene_id,
  mart=mart)

RES$Entrez = genes$entrezgene_id[ match(RES$gene_id, genes$ensembl_gene_id) ]

library(limma)
# gene set/gene onthology:
x = 1-RES$truth
names(x) = RES$Entrez
x = x[x < 0.01]
go.de = goana(names(x), coef = 1, geneid = names(x), species="Hs")
top = topGO(go.de, number = Inf)
# GO terms with pval < 0.01
top_GO = rownames(top)[top$P.DE < 0.01]

GO_res = apply(RES[,-c(1,2, ncol(RES))], 2, function(FDR){
  x = RES$Entrez[ FDR < 0.01 ]
  go.de.BAN = goana(x, coef = 1, geneid = x, FDR = 0.01, species="Hs")
  top.BAN = topGO(go.de.BAN, number = Inf)
  
  top.BAN_GO     = rownames(top.BAN)[ order(top.BAN$P.DE)[1:length(top_GO)] ]
  top.BAN_GO_FDR = rownames(top.BAN)[ order(top.BAN$P.DE)[1:length(top_GO_FDR)] ]
  
  c(sum(top_GO %in% top.BAN_GO) )
})

GO_res_05 = apply(RES[,-c(1,2, ncol(RES))], 2, function(FDR){
  x = RES$Entrez[ FDR < 0.05 ]
  go.de.BAN = goana(x, coef = 1, geneid = x, FDR = 0.01, species="Hs")
  top.BAN = topGO(go.de.BAN, number = Inf)
  
  top.BAN_GO     = rownames(top.BAN)[ order(top.BAN$P.DE)[1:length(top_GO)] ]
  top.BAN_GO_FDR = rownames(top.BAN)[ order(top.BAN$P.DE)[1:length(top_GO_FDR)] ]
  
  c(sum(top_GO %in% top.BAN_GO) )
})

res = cbind( round(apply(ranking[RES$truth ==1,], 2, median)), # 
             #apply(ranking[RES$truth ==1,], 2, mean),
             sapply( preds, function(x){ performance(x, measure="auc")@y.values[[1]] }), # AUC
             sapply( preds, function(x){ performance(x, measure="auc", fpr.stop=0.1)@y.values[[1]] }), # pAUC 0.1
             sapply( preds, function(x){ performance(x, measure="auc", fpr.stop=0.2)@y.values[[1]] }), # pAUC 0.2
             top_100, top_200,
             top_GO_01 = GO_res/length(top_GO), top_GO_05 = GO_res_05/length(top_GO)
)
ord = 1:4 # order(res[,1])
ord
library(xtable)
xtable(res[ord,], digits = 2, colnames = FALSE)
round(res[ord,], 2)
# SABLE S11