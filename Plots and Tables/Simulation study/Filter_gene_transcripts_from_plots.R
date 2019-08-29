# load estimated transcript-level counts:
library(data.table); library(tximport)
files = paste0("/6vs6_sim/aligned_reads/sample",c(1:12),"/quant.sf")

txi <- tximport(files = files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
counts <- txi$counts

# Load gene-transcript matching:
gene_tr_id = data.table::fread("DTU_Sim_6vs6/simulation_details.txt", header = TRUE)
gene_tr_id = data.frame(gene_id = gene_tr_id$gene_id, transcript_id = gene_tr_id$transcript_id )

# Match gene and transcript ids:
matches = match( rownames(counts), gene_tr_id$transcript_id)
gene_id = gene_tr_id$gene_id[matches]

# Select genes with at least 20 counts:
tot_counts_per_tr = rowSums(counts)
counts_per_gene = split(tot_counts_per_tr, gene_id)
tot_counts_per_gene = sapply(counts_per_gene, sum)
Min_20_gene_counts = tot_counts_per_gene >= 20
mean(Min_20_gene_counts)
# 0.7184223
Min_20_gene_counts = names(Min_20_gene_counts)[Min_20_gene_counts]

save(Min_20_gene_counts, file = "Min_20_gene_counts.RData")

counts_Sel = counts[gene_id %in% Min_20_gene_counts, ]
# Select transcripts with at least 10 counts:
Min_10_transcript_counts = rownames(counts_Sel)[ rowSums(counts_Sel) >= 10]
save(Min_10_transcript_counts, file = "Min_10_transcript_counts.RData")

# Select genes with at least 1200 counts (100 per sample on average):
Min_1200_gene_counts = tot_counts_per_gene >= 1200
Min_1200_gene_counts = names(Min_1200_gene_counts)[Min_1200_gene_counts]

save(Min_1200_gene_counts, file = "Min_1200_gene_counts.RData")

# store gene-level counts for each gene (needed for the stratification plot):
save(tot_counts_per_gene, file = "tot_counts_per_gene.RData")
