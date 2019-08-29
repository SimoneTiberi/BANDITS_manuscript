# Load salmon estmimated counts
library(data.table)
library(tximport)
data_dir = "null_data/aligned_reads/"
sample_names = c("P_3", "P_6", "P_8",  "P_1",  "P_4",  "P_5")
files = quant_files = file.path(data_dir, sample_names, "quant.sf")

txi <- tximport(files = files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
counts <- txi$counts

# Load truth table:
load("null_data/gene_transcript_gtf.RData")
gene_tr_id = data.frame(gene_id = X$gene_id, transcript_id = X$transcript_id, stringsAsFactors = FALSE)

# Match gene and transcript ids:
matches = match( rownames(counts), gene_tr_id$transcript_id)
gene_id = gene_tr_id$gene_id[matches]

# Select genes with at least 20 counts:
tot_counts_per_tr = rowSums(counts)
counts_per_gene = split(tot_counts_per_tr, gene_id)
tot_counts_per_gene = sapply(counts_per_gene, sum)
Min_20_gene_counts = tot_counts_per_gene >= 20
mean(Min_20_gene_counts)
Min_20_gene_counts = names(Min_20_gene_counts)[Min_20_gene_counts]

save(Min_20_gene_counts, file = "Min_20_gene_counts.RData")

counts_Sel = counts[gene_id %in% Min_20_gene_counts, ]
# Select transcripts with at least 10 counts:
Min_10_transcript_counts = rownames(counts_Sel)[ rowSums(counts_Sel) >= 10]
save(Min_10_transcript_counts, file = "Min_10_transcript_counts.RData")

# store gene-level counts for each gene (needed for the stratification plot):
save(tot_counts_per_gene, file = "tot_counts_per_gene.RData")
