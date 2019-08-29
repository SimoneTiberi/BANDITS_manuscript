########################################################################################################################################################################
# STAR
########################################################################################################################################################################
# GRCh37.71
# genome index directory:
GDIR="/3vs3_sim/STAR/sjdbOverhang100"

# path to STAR
star="path_to_STAR"

# Generate Genome index
$star --runMode genomeGenerate --runThreadN 12 --genomeDir $GDIR  \
	--genomeFastaFiles ensembl_Homo_sapiens.GRCh37.71.fa \
	--sjdbGTFfile      Homo_sapiens.GRCh37.71.sorted_matchFilteredFasta\&Filtered.gtf \
	--sjdbOverhang 100

# Homo_sapiens.GRCh37.71.sorted_matchFilteredFasta\&Filtered.gtf is the gtf file created via 'Clean gtf and fasta files.R' script

# input files
fastqDir="/3vs3_sim/reads"

# output files directory:
outDir="3vs3_sim/aligned_reads"

cd $outDir

# align reads for our 6 samples:
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
