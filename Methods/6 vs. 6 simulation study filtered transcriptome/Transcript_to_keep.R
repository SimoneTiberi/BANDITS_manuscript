# Generate the transcripts to keep (i.e., not filtered) from all methods.

library(BANDITS)
data_dir = "6vs6_sim/aligned_reads/"

gene_tr_id = data.table::fread("DTU_Sim_6vs6/simulation_details.txt", header = TRUE)
gene_tr_id = cbind( as.character(gene_tr_id$gene_id), as.character(gene_tr_id$transcript_id) )

sample_names = paste0("sample", seq_len(12))
quant_files = file.path(data_dir, sample_names, "quant.sf")

library(tximport)
txi = tximport(files = quant_files, type = "salmon", txOut = TRUE, dropInfReps=TRUE)
counts = txi$counts

transcripts_to_keep = filter_transcripts(gene_to_transcript = gene_tr_id,
                                         transcript_counts = counts, 
                                         min_transcript_proportion = 0.01,
                                         min_transcript_counts = 10, 
                                         min_gene_counts = 20)

save(transcripts_to_keep, file = "transcripts_to_keep.RData")