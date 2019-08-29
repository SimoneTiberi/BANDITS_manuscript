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
	--sjdbOverhang 77 --limitGenomeGenerateRAM 310000000000

# input files
fastqDir="/null_data/reads"

# output files directory:
outDir="null_data/aligned_reads"

cd $outDir

# align reads for our 6 samples:
$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR493937_1.fastq $fastqDir/SRR493937_2.fastq \
--outFileNamePrefix P_1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR493941_1.fastq $fastqDir/SRR493941_2.fastq \
--outFileNamePrefix P_3 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR493945_1.fastq $fastqDir/SRR493945_2.fastq \
--outFileNamePrefix P_4 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR493949_1.fastq $fastqDir/SRR493949_2.fastq \
--outFileNamePrefix P_5 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR493953_1.fastq $fastqDir/SRR493953_2.fastq \
--outFileNamePrefix P_6 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM

$star --runMode alignReads --runThreadN 20 --genomeDir $GDIR --readFilesIn $fastqDir/SRR493957_1.fastq $fastqDir/SRR493957_2.fastq \
--outFileNamePrefix P_8 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM
########################################################################################################################################################################
# Salmon
########################################################################################################################################################################
# use gffread to create a (transcript) fasta file from the gtf and (genome) fasta files used to run STAR:
gffread -w STAR_fasta_tr.fa -g Homo_sapiens.GRCh38.dna.toplevel.fa Homo_sapiens.GRCh38.92.chr.gtf

fasta="STAR_fasta_tr.fa"

# path to salmon
salmon="path_to_salmon"

# use the alignment output of STAR to compute the equivalence classes and quantify transcript abundance estimates with Salmon:
$salmon quant -t $fasta -l A -a P_1Aligned.toTranscriptome.out.bam -o P_1 -p 20 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a P_3Aligned.toTranscriptome.out.bam -o P_3 -p 20 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a P_4Aligned.toTranscriptome.out.bam -o P_4 -p 20 --dumpEq --numBootstraps 100

$salmon quant -t $fasta -l A -a P_5Aligned.toTranscriptome.out.bam -o P_5 -p 20 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a P_6Aligned.toTranscriptome.out.bam -o P_6 -p 20 --dumpEq --numBootstraps 100
$salmon quant -t $fasta -l A -a P_8Aligned.toTranscriptome.out.bam -o P_8 -p 20 --dumpEq --numBootstraps 100
