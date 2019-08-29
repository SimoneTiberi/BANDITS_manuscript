# Fist align reads with Bowtie2, then run the following script.

# setup:
# git clone https://github.com/mqbssppe/cjBitSeq.git
# cd cjBitSeq
# make

# install the following R packages:
# parallel, doParallel, Matrix, foreach, doMC, fields.

cjbs="/software/cjBitSeq"
PATH=/software/cjBitSeq/:$PATH
# CHANGE R PATH to the lastest version (with packages above installed):
PATH=/software/R/R-3.6.0/bin/:$PATH

# export R libraries:
alias R='R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/3.6/ /usr/local/R/R-3.6.0/bin/R'
export R_LIBS=/home/stiber/R/x86_64-pc-linux-gnu-library/3.6/

PATH=/software/BitSeq/:$PATH
########################################################################################################################
# the reference transcriptome:
fasta="STAR_fasta_tr.fa"

# the reference genome:
gtf="Homo_sapiens.GRCh38.92.chr.gtf"

# input aligned reads folder:
input_dir="null_data/aligned_reads_bowtie2"
########################################################################################################################
# Run parseAlignment from Bitseq:
parseAlignment $input_dir/sample1_BT2.sam -o conditionA_1.prob --trSeqFile $fasta --uniform
parseAlignment $input_dir/sample2_BT2.sam -o conditionA_2.prob --trSeqFile $fasta --uniform
parseAlignment $input_dir/sample3_BT2.sam -o conditionA_3.prob --trSeqFile $fasta --uniform
parseAlignment $input_dir/sample4_BT2.sam -o conditionB_1.prob --trSeqFile $fasta --uniform
parseAlignment $input_dir/sample5_BT2.sam -o conditionB_2.prob --trSeqFile $fasta --uniform
parseAlignment $input_dir/sample6_BT2.sam -o conditionB_3.prob --trSeqFile $fasta --uniform
# parseAlignment takes about 15 min per sample.

# note: when computing non-uniform read distribution (without option --uniform) option --procN allows parallelisation.
# We advice not using more than default 3 CPUs as this tends to decrease the overall performance.

########################################################################################################################
# Annotation pre-processing:

R CMD BATCH \
'--args STAR_fasta_tr.fa Homo_sapiens.GRCh38.92.chr.gtf out_ex' \
$cjbs/prepareAnnotationForcjBitSeq.R

# Alternatively (if no output is produced), run: 'cjBitSeq_PrepareAnnotation.R' script

########################################################################################################################
# Run DTU cjBitSeq method:
# Note: you have to be in the folder where the .prob files are!
# Note: cjBitSeq looks for "conditionA" and "conditionB" names: the .prob files need to follow that structure.

# Limit cjBitSeq usage to 12 cores:
NSLOTS=12
#!/bin/bash
export PARALLEL=-j$NSLOTS
echo $PARALLEL

$cjbs/cjBitSeqWithinGene output_DTU out_ex 12 conditionA_1.prob conditionA_2.prob conditionA_3.prob C conditionB_1.prob conditionB_2.prob conditionB_3.prob
