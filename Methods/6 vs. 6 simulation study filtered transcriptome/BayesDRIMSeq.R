source('BayesDRIMSeq_fun.R')

# Load salmon estmimated counts
library(data.table)
library(tximport)
files = paste0("/6vs6_sim/aligned_reads/sample",c(1:12),"/quant.sf")

txi <- tximport(files = files, type = "salmon", txOut = TRUE, dropInfReps = TRUE)
counts <- txi$counts

### TRANSCRIPTS FILTERING:
load("transcripts_to_keep.RData")
counts = counts[rownames(counts) %in% transcripts_to_keep, ]
# ~1/3 of the counts kept
###

# Load truth table:
truth_dir="DTU_Sim_6vs6/simulation_details.txt"
truth = fread(truth_dir, header = TRUE, sep="\t")

matches = match( rownames(counts), truth$transcript_id)

colnames(counts) = paste("sample", 1:12)

gene_id = as.character( truth$gene_id[matches] )

library('LaplacesDemon')
library('foreach')
library('doMC')

set.seed(61217)

nGenes <- length(unique(gene_id))
genes <- as.factor(gene_id)

nTranscripts <- length(genes)				# total number of transcripts
nSamplesA <- 6						# number of replicates for 1st condition 
nSamplesB <- 6						# number of replicates for 2nd condition
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
