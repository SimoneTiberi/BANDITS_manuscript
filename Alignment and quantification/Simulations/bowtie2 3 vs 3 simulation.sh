########################################################################################################################
# input files
fastqDir="/3vs3_sim/reads"

# output files for the cjBitSeq output
outDir="3vs3_sim/aligned_reads_bowtie2"

# the reference transcriptome:
fasta="Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa"
# Homo_sapiens.GRCh37.71.sorted_matchFilteredFasta\&Filtered.gtf is the gtf file created via 'Clean gtf and fasta files.R' script

bowtie2="path_to_bowtie2"

########################################################################################################################
# Build reference:
$bowtie2/bowtie2-build -f $fasta reference_bowtie2

# Align reads
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/sample_1_1.fq.gz -2 $fastqDir/sample_1_2.fq.gz -S sample1_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/sample_2_1.fq.gz -2 $fastqDir/sample_2_2.fq.gz -S sample2_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/sample_3_1.fq.gz -2 $fastqDir/sample_3_2.fq.gz -S sample3_BT2.sam 

bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/sample_4_1.fq.gz -2 $fastqDir/sample_4_2.fq.gz -S sample4_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/sample_5_1.fq.gz -2 $fastqDir/sample_5_2.fq.gz -S sample5_BT2.sam 
bowtie2 -q -k 100 --verbose --no-mixed --no-discordant --threads 12 -x reference_bowtie2 -1 $fastqDir/sample_6_1.fq.gz -2 $fastqDir/sample_6_2.fq.gz -S sample6_BT2.sam 
