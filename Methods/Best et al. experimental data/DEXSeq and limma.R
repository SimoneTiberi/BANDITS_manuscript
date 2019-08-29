# in R
R

system.file("python_scripts", package = "DEXSeq")

# /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts
#################################################################################################################################################################
# in bash:

# install HTSeq:
# pip install HTSeq --user

# dexseq_prepare_annotation.py:
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
Homo_sapiens.GRCh38.92.chr.gtf DEXSeq.gff

# dexseq_count.py:
bam_dir="/Best_data/aligned_reads"
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/control_1Aligned.sortedByCoord.out.bam sample1.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/control_2Aligned.sortedByCoord.out.bam sample2.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/control_3Aligned.sortedByCoord.out.bam sample3.txt -p yes -f bam -r pos

python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/treated_1Aligned.sortedByCoord.out.bam sample4.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/treated_2Aligned.sortedByCoord.out.bam sample5.txt -p yes -f bam -r pos
python3 /home/stiber/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py DEXSeq.gff $bam_dir/treated_3Aligned.sortedByCoord.out.bam sample6.txt -p yes -f bam -r pos

#################################################################################################################################################################
# DEXSeq in R
R

set.seed(61217)

library(DEXSeq)

sampleTable = data.frame( row.names = paste0("sample", 1:6),
                          condition = c("A", "A", "A",
                                        "B", "B", "B" ) )
sampleTable

names = paste0("sample",1:6,".txt")
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

names = paste0("sample",1:6,".txt")
data = lapply(names, data.table::fread, header =F)

# I create the counts table.
counts = ( cbind(data[[1]][,2], data[[2]][,2], data[[3]][,2], 
                 data[[4]][,2], data[[5]][,2], data[[6]][,2]) )
rownames(counts) = data[[1]][,1][[1]]
colnames(counts) = paste0("sample",1:6)

dge <- DGEList(counts=counts)

GeneID = substring(rownames(counts), 1, 15)
dge$genes <- data.frame(GeneID=GeneID)

# I do the filtering before.
A <- rowSums(dge$counts)
dge <- dge[A>=10,, keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)
design = cbind(rep(1,6), c(rep(0, 3), rep(1, 3)))

v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)

ex <- diffSplice(fit, geneid="GeneID")

res = topSplice(ex, coef=2, test="simes", number = 10^6) # I select all genes

save(res, file = "limma results.RData")
