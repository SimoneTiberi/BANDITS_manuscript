set.seed(61217)

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
truth = data.frame( gene_id = as.character(X$gene_id), transcript_id = as.character(X$transcript_id) )

matches = match( rownames(counts), truth$transcript_id)

colnames(counts) = paste("sample", 1:6)

gene_id = as.character( truth$gene_id[matches] )

library(DEXSeq)

sampleTable = data.frame( row.names = paste0("sample", 1:6),
                          condition = c("A", "A", "A",
                                        "B", "B", "B") )

dxd = DEXSeqDataSet(countData = round( counts ),
                    sampleData=sampleTable,
                    design= ~ sample + exon + condition:exon,
                    featureID = rownames(counts),
                    groupID = gene_id, 
                    transcripts = rownames(counts))

# normalization:
dxd = estimateSizeFactors( dxd )

library(BiocParallel)
BPPARAM = MulticoreParam(workers=12)

# dispersion estimate:
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)

# NB test:
dxd = testForDEU( dxd, reduced = ~ sample + exon, BPPARAM=BPPARAM)
# Error in designAndArgChecker(object, betaPrior = FALSE) : 
# full model matrix is less than full rank

res = DEXSeqResults( dxd, independentFiltering = FALSE )

qval = perGeneQValue(res)

res_gene = data.frame(gene = names(qval), qval)

save(res, res_gene, file="DEXSeq_TECs results.RData")
