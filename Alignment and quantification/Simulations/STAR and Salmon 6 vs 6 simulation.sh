########################################################################################################################################################################
# STAR
########################################################################################################################################################################
# GRCh37.71
# genome index directory:
GDIR="/6vs6_sim/STAR/sjdbOverhang100"

# path to STAR
star="path_to_STAR"

# Generate Genome index
$star --runMode genomeGenerate --runThreadN 12 --genomeDir $GDIR  \
	--genomeFastaFiles ensembl_Homo_sapiens.GRCh37.71.fa \
	--sjdbGTFfile      Homo_sapiens.GRCh37.71.sorted_matchFilteredFasta\&Filtered.gtf \
	--sjdbOverhang 100

# Homo_sapiens.GRCh37.71.sorted_matchFilteredFasta\&Filtered.gtf is the gtf file created via 'Clean gtf and fasta files.R' script

# input files
fastqDir="/6vs6_sim/reads"

# output files directory:
outDir="6vs6_sim/aligned_reads"

cd $outDir

# align reads for our 12 samples:
$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_1_1.fq.gz) <(zcat $fastqDir/sample_1_2.fq.gz) \
--outFileNamePrefix sample1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_2_1.fq.gz) <(zcat $fastqDir/sample_2_2.fq.gz) \
--outFileNamePrefix sample2 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_3_1.fq.gz) <(zcat $fastqDir/sample_3_2.fq.gz) \
--outFileNamePrefix sample3 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_4_1.fq.gz) <(zcat $fastqDir/sample_4_2.fq.gz) \
--outFileNamePrefix sample4 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_5_1.fq.gz) <(zcat $fastqDir/sample_5_2.fq.gz) \
--outFileNamePrefix sample5 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_6_1.fq.gz) <(zcat $fastqDir/sample_6_2.fq.gz) \
--outFileNamePrefix sample6 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_7_1.fq.gz) <(zcat $fastqDir/sample_7_2.fq.gz) \
--outFileNamePrefix sample7 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_8_1.fq.gz) <(zcat $fastqDir/sample_8_2.fq.gz) \
--outFileNamePrefix sample8 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_9_1.fq.gz) <(zcat $fastqDir/sample_9_2.fq.gz) \
--outFileNamePrefix sample9 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_10_1.fq.gz) <(zcat $fastqDir/sample_10_2.fq.gz) \
--outFileNamePrefix sample10 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_11_1.fq.gz) <(zcat $fastqDir/sample_11_2.fq.gz) \
--outFileNamePrefix sample11 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 12 --genomeDir $GDIR --readFilesIn <(zcat $fastqDir/sample_12_1.fq.gz) <(zcat $fastqDir/sample_12_2.fq.gz) \
--outFileNamePrefix sample12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

########################################################################################################################################################################
# Salmon
########################################################################################################################################################################
# path to salmon
salmon="path_to_salmon"

# the reference transcriptome:
fasta="Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa"

# use the alignment output of STAR to compute the equivalence classes and quantify transcript abundance estimates with Salmon:
$salmon quant -t $fasta -l A -a sample1Aligned.toTranscriptome.out.bam -o sample1 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample2Aligned.toTranscriptome.out.bam -o sample2 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample3Aligned.toTranscriptome.out.bam -o sample3 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample4Aligned.toTranscriptome.out.bam -o sample4 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample5Aligned.toTranscriptome.out.bam -o sample5 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample6Aligned.toTranscriptome.out.bam -o sample6 -p 12 --dumpEq --numBootstraps 100

$salmon quant -t $fasta -l A -a sample7Aligned.toTranscriptome.out.bam -o sample7 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample8Aligned.toTranscriptome.out.bam -o sample8 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample9Aligned.toTranscriptome.out.bam -o sample9 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample10Aligned.toTranscriptome.out.bam -o sample10 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample11Aligned.toTranscriptome.out.bam -o sample11 -p 12 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a sample12Aligned.toTranscriptome.out.bam -o sample12 -p 12 --dumpEq --numBootstraps 100
