library(BANDITS)

data_dir = "6vs6_sim/aligned_reads/kallisto"

gene_tr_id = data.table::fread("DTU_Sim_6vs6/simulation_details.txt", header = TRUE)
gene_tr_id = cbind( as.character(gene_tr_id$gene_id), as.character(gene_tr_id$transcript_id) )

sample_names = paste0("sample", seq_len(12))
quant_files = file.path(data_dir, sample_names, "quant.sf")

library(tximport)
txi = tximport(files = quant_files, type = "kallisto", txOut = TRUE, dropInfReps=TRUE)
counts = txi$counts

samples_design = data.frame(sample_id = sample_names,
                            group = rep( c("A", "B"), each = 6))

eff_len = eff_len_compute(x_eff_len = txi$length)

kallisto_equiv_classes = file.path(data_dir, "ec", sample_names, "pseudoalignments.ec")
kallisto_equiv_counts  = file.path(data_dir, "ec", sample_names, "pseudoalignments.tsv")

input_data = create_data(salmon_or_kallisto = "kallisto",
                         gene_to_transcript = gene_tr_id,
                         kallisto_equiv_classes = kallisto_equiv_classes,
                         kallisto_equiv_counts = kallisto_equiv_counts,
                         kallisto_counts = counts,
                         eff_len = eff_len, n_cores = 12)

input_data = filter_genes(input_data, min_counts_per_gene = 10)

set.seed(61217)
precision = prior_precision(gene_to_transcript = gene_tr_id,
                            transcript_counts = counts, n_cores = 12)


set.seed(61217)
results = test_DTU(BANDITS_data = input_data,
                   precision = precision$prior,
                   samples_design = samples_design,
                   group_col_name = "group",
                   R = 10^4, burn_in = 2*10^3, n_cores = 12,
                   gene_to_transcript = gene_tr_id)

save(results, file = "BANDITS kallisto results.RData")
