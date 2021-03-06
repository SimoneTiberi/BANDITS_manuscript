# in R
R

system.file("python_scripts", package = "DEXSeq")

# /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts
#################################################################################################################################################################
# in bash:

# install HTSeq:
# pip install HTSeq --user

# Note: Filtering transcripts was done before aligning reads:
# reads were aligned to the filtered transcriptome/genome via STAR

# dexseq_prepare_annotation.py:
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
Homo_sapiens.GRCh37.71.sorted_matchFilteredFasta\&Filtered_FILT01.gtf DEXSeq.gff

# dexseq_count.py:
bam_dir="/6vs6_sim/aligned_reads_filt01"
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample1Aligned.sortedByCoord.out.bam sample1.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample2Aligned.sortedByCoord.out.bam sample2.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample3Aligned.sortedByCoord.out.bam sample3.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample4Aligned.sortedByCoord.out.bam sample4.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample5Aligned.sortedByCoord.out.bam sample5.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample6Aligned.sortedByCoord.out.bam sample6.txt -p yes -f bam -r pos

python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample7Aligned.sortedByCoord.out.bam sample7.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample8Aligned.sortedByCoord.out.bam sample8.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample9Aligned.sortedByCoord.out.bam sample9.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample10Aligned.sortedByCoord.out.bam sample10.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample11Aligned.sortedByCoord.out.bam sample11.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/sample12Aligned.sortedByCoord.out.bam sample12.txt -p yes -f bam -r pos

#################################################################################################################################################################
# DEXSeq in R
R

set.seed(61217)

library(DEXSeq)

sampleTable = data.frame( row.names = paste0("sample", 1:12),
                          condition = c("A", "A", "A", "A", "A", "A",
                                        "B", "B", "B", "B", "B", "B" ) )
sampleTable

names = paste0("sample",1:12,".txt")
dxd = DEXSeqDataSetFromHTSeq(names,
                             sampleData=sampleTable,
                             flattenedfile="DEXSeq.gff" )
dxd
head( counts(dxd), 5 )
head( featureCounts(dxd), 5 )

# normalization:
dxd = estimateSizeFactors( dxd )

library(BiocParallel)
BPPARAM = MulticoreParam(workers=12)

# dispersion estimate:
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)

# NB test:
dxd = testForDEU( dxd, reduced = ~ sample + exon, BPPARAM=BPPARAM)

res = DEXSeqResults( dxd, independentFiltering = FALSE )

qval = perGeneQValue(res)

res_gene = data.frame(gene = names(qval), qval)

save(res, res_gene, file="DEXSeq results.RData")

#################################################################################################################################################################
# limma in R
R

library(DEXSeq); library(limma); library(edgeR)

names = paste0("sample",1:12,".txt")
data = lapply(names, data.table::fread, header =F)

# I create the counts table.
counts = ( cbind(data[[1]][,2], data[[2]][,2], data[[3]][,2], data[[4]][,2], data[[5]][,2], data[[6]][,2], 
                 data[[7]][,2], data[[8]][,2], data[[9]][,2], data[[10]][,2], data[[11]][,2], data[[12]][,2]) )
rownames(counts) = data[[1]][,1][[1]]
colnames(counts) = paste0("sample",1:12)

dge <- DGEList(counts=counts)

GeneID = substring(rownames(counts), 1, 15)
dge$genes <- data.frame(GeneID=GeneID)

# I do the filtering before.
A <- rowSums(dge$counts)
dge <- dge[A>=10,, keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)
design = cbind(rep(1,12), c(rep(0, 6), rep(1, 6)))

v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)

ex <- diffSplice(fit, geneid="GeneID")

res = topSplice(ex, coef=2, test="simes", number = 10^6) # I select all genes

save(res, file = "limma results.RData")
