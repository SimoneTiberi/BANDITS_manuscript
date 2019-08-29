set.seed(61217)

# Load salmon estmimated counts
library(data.table)
library(tximport)
files = paste0("/6vs6_sim/aligned_reads/sample",c(1:12),"/quant.sf")

txi <- tximport(files = files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
counts <- txi$counts

# Load truth table:
truth_dir="DTU_Sim_6vs6/simulation_details.txt"
truth = fread(truth_dir, header = TRUE, sep="\t")

matches = match( rownames(counts), truth$transcript_id)

colnames(counts) = paste("sample", 1:12)

gene_id = as.character( truth$gene_id[matches] )

library(DEXSeq)

sampleTable = data.frame( row.names = paste0("sample", 1:12),
                          condition = c("A", "A", "A", "A", "A", "A",
                                        "B", "B", "B", "B", "B", "B") )

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
