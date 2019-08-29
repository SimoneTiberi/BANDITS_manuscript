########################################################################################################################
# input files
fastqDir="/null_data/reads"

# output files for the cjBitSeq output
outDir="null_data/aligned_reads_bowtie2"

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
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR493937_1.fastq $fastqDir/SRR493937_2.fastq -S sample1_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR493941_1.fastq $fastqDir/SRR493941_2.fastq -S sample2_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR493945_1.fastq $fastqDir/SRR493945_2.fastq -S sample3_BT2.sam 

bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR493949_1.fastq $fastqDir/SRR493949_2.fastq -S sample4_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR493953_1.fastq $fastqDir/SRR493953_2.fastq -S sample5_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/SRR493957_1.fastq $fastqDir/SRR493957_2.fastq -S sample6_BT2.sam 
