########################################################################################################################################################################
# STAR
########################################################################################################################################################################
# GRCh37.71
# genome index directory:
GDIR="/STAR/sjdbOverhang100"

# path to STAR
star="path_to_STAR"

# Generate Genome index
$star --runMode genomeGenerate --runThreadN 12 --genomeDir $GDIR  \
	--genomeFastaFiles Homo_sapiens.GRCh38.dna.toplevel.fa \
	--sjdbGTFfile      Homo_sapiens.GRCh38.92.chr.gtf \
	--sjdbOverhang 100 --limitGenomeGenerateRAM 310000000000

# input files
fastqDir="/Best_data/reads"

# output files directory:
outDir="Best_data/aligned_reads"

cd $outDir

# align reads for our 6 samples:
$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR1513329_1.fastq $fastqDir/SRR1513329_2.fastq \
--outFileNamePrefix control_1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR1513330_1.fastq $fastqDir/SRR1513330_2.fastq \
--outFileNamePrefix control_2 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR1513331_1.fastq $fastqDir/SRR1513331_2.fastq \
--outFileNamePrefix control_3 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR1513332_1.fastq $fastqDir/SRR1513332_2.fastq \
--outFileNamePrefix treated_1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR1513333_1.fastq $fastqDir/SRR1513333_2.fastq \
--outFileNamePrefix treated_2 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR1513334_1.fastq $fastqDir/SRR1513334_2.fastq \
--outFileNamePrefix treated_3 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM
########################################################################################################################################################################
# Salmon
########################################################################################################################################################################
# use gffread to create a (transcript) fasta file from the gtf and (genome) fasta files used to run STAR:
gffread -w STAR_fasta_tr.fa -g Homo_sapiens.GRCh38.dna.toplevel.fa Homo_sapiens.GRCh38.92.chr.gtf

fasta="STAR_fasta_tr.fa"

# path to salmon
salmon="path_to_salmon"

# use the alignment output of STAR to compute the equivalence classes and quantify transcript abundance estimates with Salmon:
$salmon quant -t $fasta -l A -a control_1Aligned.toTranscriptome.out.bam -o control_1 -p 20 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a control_2Aligned.toTranscriptome.out.bam -o control_2 -p 20 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a control_3Aligned.toTranscriptome.out.bam -o control_3 -p 20 --dumpEq --numBootstraps 100

$salmon quant -t $fasta -l A -a treated_1Aligned.toTranscriptome.out.bam -o treated_1 -p 20 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a treated_2Aligned.toTranscriptome.out.bam -o treated_2 -p 20 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a treated_3Aligned.toTranscriptome.out.bam -o treated_3 -p 20 --dumpEq --numBootstraps 100
