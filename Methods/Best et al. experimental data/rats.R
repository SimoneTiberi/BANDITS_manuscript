# installing wasabi from binary:
# R CMD INSTALL -l /home/stiber/R/x86_64-pc-linux-gnu-library/3.6 wasabi-master/

# installing rats from binary:
# R CMD INSTALL -l /home/stiber/R/x86_64-pc-linux-gnu-library/3.6 RATS-master/

# in R:
R

library(wasabi);
library(rats)

# Load salmon estmimated counts
library(data.table)
library(tximport)
data_dir = "Best_data/aligned_reads/"
samples = c("control_1", "control_2", "control_3", "treated_1", "treated_2", "treated_3")
files = quant_files = file.path(data_dir, sample_names, "quant.sf")

txi <- tximport(files = files, type = "salmon", txOut = TRUE, dropInfReps = FALSE)
counts <- txi$counts

# Load truth table:
load("Best_data/gene_transcript_gtf.RData")
truth = data.frame( gene_id = as.character(X$gene_id), transcript_id = as.character(X$transcript_id) )

matches = match( rownames(counts), truth$transcript_id)

colnames(counts) = paste("sample", 1:6)

gene_id = as.character( truth$gene_id[matches] )

myannot  <- data.frame( target_id = rownames(counts), parent_id = gene_id )        # Transcript and gene IDs for the above data.

mycond_A = mycond_B = list()
for(i in 1:3){
  mycond_A[[i]] <- data.table( target_id = myannot$target_id, txi$infReps[[i]]  )
}
for(i in 4:6){
  mycond_B[[i-3]] <- data.table( target_id = myannot$target_id, txi$infReps[[i]]  )
}

# WITH bootstrap replicates:
mydtu_bootstrap <- call_DTU(annot = myannot, 
                            boot_data_A = mycond_A,
                            boot_data_B = mycond_B,
                            threads = 12, qboot  = TRUE)

#mydtu_bootstrap
res_boot = mydtu_bootstrap$Genes
res_boot_transcripts = mydtu_bootstrap$Transcripts

save(res_boot, res_boot_transcripts, file = "rats results.RData")
