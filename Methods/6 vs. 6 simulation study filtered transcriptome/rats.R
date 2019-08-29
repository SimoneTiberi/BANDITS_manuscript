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
files = paste0("/home/stiber/SIM5_Simone/FILTERED/STAR_MultiMapping/sample",c(1:12),"/quant.sf")
#file.exists(files)

txi <- tximport(files = files, type = "salmon", txOut = TRUE, dropInfReps = FALSE)
counts <- txi$counts

### TRANSCRIPTS FILTERING:
load("transcripts_to_keep.RData")
SEL_rows = rownames(counts) %in% transcripts_to_keep
counts = counts[SEL_rows, ]
# ~1/3 of the counts kept
###

# Load truth table:
truth_dir="/home/stiber/data_SIM5_Simone/simulation_details.txt"
truth = fread(truth_dir, header = TRUE, sep="\t")

matches = match( rownames(counts), truth$transcript_id)

colnames(counts) = paste("sample", 1:12)

gene_id = as.character( truth$gene_id[matches] )

library(rats)

myannot  <- data.frame( target_id = rownames(counts), parent_id = gene_id )        # Transcript and gene IDs for the above data.

mycond_A = mycond_B = list()
for(i in 1:6){
  mycond_A[[i]] <- data.table( target_id = myannot$target_id, txi$infReps[[i]][SEL_rows,]  )
}
for(i in 7:12){
  mycond_B[[i-6]] <- data.table( target_id = myannot$target_id, txi$infReps[[i]][SEL_rows,]  )
}
# SEL_rows selects transcripts to keep.

# WITH bootstrap replicates:
mydtu_bootstrap <- call_DTU(annot = myannot, 
                            boot_data_A = mycond_A,
                            boot_data_B = mycond_B,
                            threads = 12, qboot  = TRUE)

#mydtu_bootstrap
res_boot = mydtu_bootstrap$Genes
res_boot_transcripts = mydtu_bootstrap$Transcripts

save(res_boot, res_boot_transcripts, file = "rats results.RData")
