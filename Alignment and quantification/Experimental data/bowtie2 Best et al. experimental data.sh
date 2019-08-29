########################################################################################################################
# input files
fastqDir="/Best_data/reads"

# output files for the cjBitSeq output
outDir="Best_data/aligned_reads_bowtie2"

# the reference transcriptome:
# For consistency with the other methods, we use the same transcriptome used in the other tools:
# we use gffread to create a (transcript) fasta file from the gtf and (genome) fasta files used to run STAR:
gffread -w STAR_fasta_tr.fa -g Homo_sapiens.GRCh38.dna.toplevel.fa Homo_sapiens.GRCh38.92.chr.gtf

fasta="STAR_fasta_tr.fa"

bowtie2="path_to_bowtie2"

########################################################################################################################
# Build reference:
$bowtie2/bowtie2-build -f $fasta reference_bowtie2

# Align reads
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR1513329_1.fastq $fastqDir/SRR1513329_2.fastq -S sample1_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR1513330_1.fastq $fastqDir/SRR1513330_2.fastq -S sample2_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR1513331_1.fastq $fastqDir/SRR1513331_2.fastq -S sample3_BT2.sam 

bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR1513332_1.fastq $fastqDir/SRR1513332_2.fastq -S sample4_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR1513333_1.fastq $fastqDir/SRR1513333_2.fastq -S sample5_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR1513334_1.fastq $fastqDir/SRR1513334_2.fastq -S sample6_BT2.sam 
