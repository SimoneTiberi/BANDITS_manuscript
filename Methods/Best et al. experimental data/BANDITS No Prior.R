library(BANDITS)

data_dir = "Best_data/aligned_reads/"

load("Best_data/gene_transcript_gtf.RData")
gene_tr_id = cbind( as.character(X$gene_id), as.character(X$transcript_id) )

sample_names = c("control_1", "control_2", "control_3",
                 "treated_1", "treated_2", "treated_3")
quant_files = file.path(data_dir, sample_names, "quant.sf")

library(tximport)
txi = tximport(files = quant_files, type = "salmon", txOut = TRUE, dropInfReps=TRUE)
counts = txi$counts

samples_design = data.frame(sample_id = sample_names,
                            group = rep( c("A", "B"), each = 3))

eff_len = eff_len_compute(x_eff_len = txi$length)

equiv_classes_files = file.path(data_dir, sample_names, 
                                "aux_info", "eq_classes.txt")

input_data = create_data(salmon_or_kallisto = "salmon",
                         gene_to_transcript = gene_tr_id,
                         salmon_path_to_eq_classes = equiv_classes_files,
                         eff_len = eff_len, n_cores = 12)

input_data = filter_genes(input_data, min_counts_per_gene = 10)

set.seed(61217)
results = test_DTU(BANDITS_data = input_data,
                   samples_design = samples_design,
                   group_col_name = "group",
                   R = 10^4, burn_in = 2*10^3, n_cores = 12,
                   gene_to_transcript = gene_tr_id)

save(results, file = "BANDITS results NoPrior.RData")
