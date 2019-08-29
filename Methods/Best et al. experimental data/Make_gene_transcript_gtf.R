# use the gtf to make a correspondance between gene id and transcript id:
library(rtracklayer)
gtf_import_rtracklayer_original = import("Homo_sapiens.GRCh38.92.chr.gtf", format = "gtf")

transcript = gtf_import_rtracklayer_original$transcript_id
gene = gtf_import_rtracklayer_original$gene_id

gene_sel       = gene[is.na(transcript) ==FALSE]
transcript_sel = transcript[is.na(transcript) ==FALSE]

X = data.frame(gene_id = gene_sel, transcript_id = transcript_sel)

save(X, file = "Best_data/gene_transcript_gtf.RData")
