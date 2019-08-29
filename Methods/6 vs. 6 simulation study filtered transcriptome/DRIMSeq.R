# Load salmon estmimated counts
library(data.table)
library(tximport)
files = paste0("/6vs6_sim/aligned_reads/sample",c(1:12),"/quant.sf")

txi <- tximport(files = files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
counts <- txi$counts

### TRANSCRIPTS FILTERING:
load("/home/stiber/SIM5_Simone/FILTERED/BANDITS/transcripts_to_keep.RData")
counts = counts[rownames(counts) %in% transcripts_to_keep, ]
# ~1/3 of the counts kept
###

# Load truth table:
truth_dir="DTU_Sim_6vs6/simulation_details.txt"
truth = fread(truth_dir, header = TRUE, sep="\t")
# or use fread
#head(truth)

matches = match( rownames(counts), truth$transcript_id)
#head( rownames(counts) )
#head( truth$transcript_id[matches] )

colnames(counts) = paste("sample", 1:12)

gene_id = as.character( truth$gene_id[matches] )

N_1 = N_2 = 6

library(DRIMSeq)
#library(edgeR)

samples <- data.frame(sample_id = colnames(counts),
                      group = c( rep("A", N_1), rep("B", N_2) ))

counts_df = data.frame(counts, gene_id = gene_id, feature_id = rownames(counts))
#head(counts_df)

samples <- data.frame(sample_id = colnames(counts_df)[1:{N_1 + N_1}],
                      group = c( rep("A", N_1), rep("B", N_2) ))
#samples

# Create a dmDSdata object
d <- dmDSdata(counts = counts_df, samples = samples)
#d

#table(samples(d)$group)
d <- dmFilter(d, min_samps_gene_expr = 1, min_samps_feature_expr = 1, min_samps_feature_prop= 0)
#plotData(d)

design_full <- model.matrix(~ group, data = samples(d))
#design_full

library(BiocParallel)
BPPARAM = MulticoreParam(workers=12)

# infer the precision parameters:
d <- dmPrecision(d, genewise_precision = TRUE, 
                 design = design_full, BPPARAM = BPPARAM)

# We fit the model
# ?dmFit
d <- dmFit(d, design = design_full, verbose = 1, BPPARAM = BPPARAM)
#d

# We test the genes
# ?dmTest
d <- dmTest(d, coef = "groupB", verbose = 1, BPPARAM = BPPARAM)
#d

results_gene	    = results(d, level = "gene")
results_trancript = results(d, level = "feature")
# Proportions has the transcript proportions for the two conditions:
props = proportions(d)

# Save results:
save(results_gene, results_trancript, props, file = "DRIMSeq results.RData")
