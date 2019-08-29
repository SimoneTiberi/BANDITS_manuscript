source('BayesDRIMSeq_fun.R')

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

library('LaplacesDemon')
library('foreach')
library('doMC')

set.seed(61217)

nGenes <- length(unique(gene_id))
genes <- as.factor(gene_id)

nTranscripts <- length(genes)				# total number of transcripts
nSamplesA <- 3						# number of replicates for 1st condition 
nSamplesB <- 3						# number of replicates for 2nd condition
nSamples <- nSamplesA + nSamplesB			# total number of replicates
# 	Simulate data.frame with counts (no DTU evidence at all)
countDataFrame <- data.frame(counts)

# 	Call the BayesDRIMSEQ function as:
myRes <- laplaceDM(
  count_data = countDataFrame, 		# data.frame of counts with dimension: nTranscripts x nSamples
  gene_data = genes, 			# factor with `nGenes` levels with length: nTranscripts
  grouping = as.factor(
    c(rep('A',nSamplesA),
      rep('B',nSamplesB))), 	# factor with 2 levels and length nSamples
  min_reads_filter = 10, 			# threshold used to filter our low expressed transcripts
  nCores = 12, 				# number of paraller workers
  lambdaRate = 0.5			# positive prior parameter \lambda
)

save(myRes, file = "BayeDRIMSeq_results.RData")
