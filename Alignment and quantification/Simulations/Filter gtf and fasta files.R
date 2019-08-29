######################################################################################################################################################
# FILTER the transcript fasta file:
######################################################################################################################################################

library(Biostrings)
fasta = readDNAStringSet("Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa")
fasta

fasta_names = names(fasta)
head(fasta_names)

fasta_tr_id = substring(fasta_names, 1, 15) # the first 15 char specify the tr id
head(fasta_tr_id)

# load the transcripts to keep, filtering was done via:
# BANDITS::transcripts_to_keep = filter_transcripts(gene_to_transcript = gene_tr_id,
#                                         transcript_counts = counts, min_transcript_proportion = 0.01,
#                                         min_transcript_counts = 10, min_gene_counts = 20)

load("transcripts_to_keep.RData")
SEL_tr = fasta_tr_id %in% transcripts_to_keep
mean(SEL_tr)
# ~1/3 of the transcripts are kept

fasta_SEL = fasta[SEL_tr]

writeXStringSet(fasta_SEL,
                file = "Homo_sapiens.GRCh37.71.SELECTED_FILT01_STARMultiMap.cdna.all.fa")


######################################################################################################################################################
# FILTER the gtf file to match the transcripts from the filtered transcript fasta file (Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa):
######################################################################################################################################################
rm(list = ls())

# load the gtf file:
library(rtracklayer, lib.loc = "/home/stiber/R/x86_64-pc-linux-gnu-library/3.6")
gtf_import_rtracklayer = import("Homo_sapiens.GRCh37.71.sorted_matchFilteredFasta&Filtered.gtf", format = "gtf")

# load the FILTERED transcript fasta file:
library(Biostrings)
fasta = readDNAStringSet("Homo_sapiens.GRCh37.71.SELECTED_FILT01_STARMultiMap.cdna.all.fa")

fasta_names = names(fasta)
head(fasta_names)

fasta_tr_id = substring(fasta_names, 1, 15) # the first 15 char specify the tr id
head(fasta_tr_id)
mean(fasta_tr_id %in% gtf_import_rtracklayer$transcript_id)
# 1
# All fasta tr ids are in the gtf tr ids

# filter 1: transcripts should also happear in the fasta file:
FINAL_SEL = gtf_import_rtracklayer$transcript_id %in% fasta_tr_id
mean(FINAL_SEL)
# ~58 % of transcripts are filtered from the gtf.

# FILTER gtf to keep only transcripts in 'FINAL_SEL'
gtf_import_rtracklayer@ranges@start = gtf_import_rtracklayer@ranges@start[FINAL_SEL]
gtf_import_rtracklayer@ranges@width = gtf_import_rtracklayer@ranges@width[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$source = gtf_import_rtracklayer@elementMetadata@listData$source[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$type = gtf_import_rtracklayer@elementMetadata@listData$type[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$score = gtf_import_rtracklayer@elementMetadata@listData$score[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$phase = gtf_import_rtracklayer@elementMetadata@listData$phase[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$gene_id = gtf_import_rtracklayer@elementMetadata@listData$gene_id[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$transcript_id = gtf_import_rtracklayer@elementMetadata@listData$transcript_id[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$exon_number = gtf_import_rtracklayer@elementMetadata@listData$exon_number[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$gene_name = gtf_import_rtracklayer@elementMetadata@listData$gene_name[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$gene_biotype = gtf_import_rtracklayer@elementMetadata@listData$gene_biotype[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$transcript_name = gtf_import_rtracklayer@elementMetadata@listData$transcript_name[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$exon_id = gtf_import_rtracklayer@elementMetadata@listData$exon_id[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@listData$protein_id = gtf_import_rtracklayer@elementMetadata@listData$protein_id[FINAL_SEL]
gtf_import_rtracklayer@seqnames = gtf_import_rtracklayer@seqnames[FINAL_SEL]
gtf_import_rtracklayer@strand = gtf_import_rtracklayer@strand[FINAL_SEL]
gtf_import_rtracklayer@elementMetadata@nrows = sum(FINAL_SEL)

gtf_import_rtracklayer

# Check that all TR ids in the gtf file also happear in the FASTA file:
mean(gtf_import_rtracklayer$transcript_id %in% fasta_tr_id)

# Check that all TR ids in the fasta file also happear in the gtf file:
mean(fasta_tr_id %in% gtf_import_rtracklayer$transcript_id )

# Export the filtered gtf file:
setwd("/home/stiber/SIM5_Simone/FILTERED/BANDITS")
export(gtf_import_rtracklayer, 'Homo_sapiens.GRCh37.71.sorted_matchFilteredFasta&Filtered_FILT01.gtf', 
       format = 'gtf') # It may sort the columns differently...LOOK for a different method to save the file!
