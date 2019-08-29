gtf="Homo_sapiens.GRCh38.92.chr.gtf"

suppa="/software/SUPPA-2.3"

export LC_ALL=C

########## To generate the events from the GTF file one has to run the following command:
python3.4 $suppa/suppa.py generateEvents -i $gtf -o events -f ioi

########## need to compute an expression_file: TPM with 1 col per sample, 1 row per transcript
# in R:
R

library(data.table); library(tximport)
data_dir = "Best_data/aligned_reads/"
samples = c("control_1", "control_2", "control_3", "treated_1", "treated_2", "treated_3")
files = quant_files = file.path(data_dir, sample_names, "quant.sf")

txi <- tximport(files = files, type = "salmon", txOut = TRUE, abundanceCol = "TPM")
TPM <- txi$abundance
colnames(TPM) = paste0("sample", 1:6)
head(TPM)

write.table(TPM[,1:3], file = "expression_file_Cond_A.csv", row.names=TRUE, col.names=TRUE, sep="\t")
write.table(TPM[,4:6], file = "expression_file_Cond_B.csv", row.names=TRUE, col.names=TRUE, sep="\t")
# \t for tab separated.

# OPEN WITH TEXTEDIT AND then manually remove quotes.

########## At the moment the PSI per transcript isoform is calculated in the following way:
python3.4 $suppa/suppa.py psiPerIsoform -g $gtf -e expression_file_Cond_A.csv -o psi_A
python3.4 $suppa/suppa.py psiPerIsoform -g $gtf -e expression_file_Cond_B.csv -o psi_B

########## Differential splicing:
python3.4 $suppa/suppa.py diffSplice -m empirical -gc -i events.ioi -p psi_A_isoform.psi psi_B_isoform.psi \
-e expression_file_Cond_A.csv expression_file_Cond_B.csv -o DTU_OUTPUT

# SUPPA2 includes an option to correct for multiple testing using the Benjamini-Hochberg method across all events
# from the same gene, as they cannot be considered to be entirely independent of each other, for which the false discovery rate (FDR) cut-off can be given as input.

# -gc, --gene-correction
# Boolean. If True, SUPPA correct the p-values by gene. (Default: False).

########## Sort results in R:
# DTU with transcripts:
# in R:
R

res_SUPPA = read.table("DTU_OUTPUT.dpsi", header = TRUE)
head(res_SUPPA); dim(res_SUPPA)

# get the gene names:              
gene_names = rownames(res_SUPPA)
gene_names = substring(gene_names, 1, 15)

res_by_gene = split(res_SUPPA[,2], gene_names)
min_p_val = lapply(res_by_gene, min, na.rm = T)
min_p_val = unlist(min_p_val)

res_SUPPA = min_p_val

save(res_SUPPA, file = "SUPPA2 results.RData")
