######################################################################################################################################################
# FILTER the transcript fasta file:
######################################################################################################################################################

# load the transcript fast file:
library(Biostrings)
fasta = readDNAStringSet("Homo_sapiens.GRCh37.71.cdna.all.fa")

x = as.character(fasta)
dups = duplicated(x)
mean(dups)
# ~7% are duplicated
head(names(fasta))

# split the fasta file names:
ss <- strsplit(names(fasta), " ")

# get the gene biotype:
biotype = sapply(ss, .subset, 5) # gene biotype
table(biotype, dups)

# get the chromosome type
coord = sapply(ss, .subset, 3) # Chromosome coordinates
head(coord)
coordss <- strsplit(coord,":") # I split the coordinates according
head(coordss)
chr <- sapply(coordss, .subset, 3)
table(chr)

# FILTER:
# keep Chromosomes: 1:22, "X" and "Y"
# keep gene biotype = protein_coding
sel_chr_geneCoding <- chr %in% c(1:22,"X","Y") & biotype=="gene_biotype:protein_coding"
mean(sel_chr_geneCoding)
# ~ 25 % of the transcripts are removed from the fasta file.

# save the filtered fasta file, also excluding duplicated transcripts (less than 1%):
writeXStringSet(unique(fasta[sel_chr_geneCoding]),
                file = "Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa")

######################################################################################################################################################
# FILTER the gtf file to match the transcripts from the filtered transcript fasta file (Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa):
######################################################################################################################################################
rm(list = ls())

# load the gtf file:
library(rtracklayer)
gtf_import_rtracklayer = import("Homo_sapiens.GRCh37.71.sorted.gtf", format = "gtf")

# load the FILTERED transcript fasta file:
library(Biostrings)
fasta = readDNAStringSet("Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa")

fasta_names = names(fasta)
head(fasta_names)

fasta_tr_id = substring(fasta_names, 1, 15) # the first 15 char specify the tr id
head(fasta_tr_id)
mean(fasta_tr_id %in% gtf_import_rtracklayer$transcript_id)
# 1
# All fasta tr ids are in the gtf tr ids

# filter 1: transcripts should also happear in the fasta file.
gtf_sel = gtf_import_rtracklayer$transcript_id %in% fasta_tr_id
mean(gtf_sel)

# filter 2: I only keep exon files, not start-end codons or CDS.
exon_sel = gtf_import_rtracklayer$type == "exon"
mean(exon_sel)

# FILTER including all 3 filters above
FINAL_SEL = gtf_sel & exon_sel
mean(FINAL_SEL)

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
export(gtf_import_rtracklayer, 'Homo_sapiens.GRCh37.71.sorted_matchFilteredFasta&Filtered.gtf', 
       format = 'gtf') # It may sort the columns differently...LOOK for a different method to save the file!
