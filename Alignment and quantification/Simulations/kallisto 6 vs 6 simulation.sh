# the reference transcriptome
fasta="Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa"

# input files
fastqDir="/6vs6_sim/reads"

# output files
kallistoDir="/6vs6_sim/aligned_reads/kallisto"

cd $kallistoDir

# index name
idx="Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa.idx_kallisto"

# Build index
kallisto index -i $idx $fasta

# quant (to quantify transcript abundance)
kallisto quant -i  $idx -o $kallistoDir/sample1 --threads 12 \
$fastqDir/sample_1_1.fq.gz $fastqDir/sample_1_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample2 --threads 12 \
$fastqDir/sample_2_1.fq.gz $fastqDir/sample_2_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample3 --threads 12 \
$fastqDir/sample_3_1.fq.gz $fastqDir/sample_3_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample4 --threads 12 \
$fastqDir/sample_4_1.fq.gz $fastqDir/sample_4_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample5 --threads 12 \
$fastqDir/sample_5_1.fq.gz $fastqDir/sample_5_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample6 --threads 12 \
$fastqDir/sample_6_1.fq.gz $fastqDir/sample_6_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample7 --threads 12 \
$fastqDir/sample_7_1.fq.gz $fastqDir/sample_7_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample8 --threads 12 \
$fastqDir/sample_8_1.fq.gz $fastqDir/sample_8_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample9 --threads 12 \
$fastqDir/sample_9_1.fq.gz $fastqDir/sample_9_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample10 --threads 12 \
$fastqDir/sample_10_1.fq.gz $fastqDir/sample_10_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample11 --threads 12 \
$fastqDir/sample_11_1.fq.gz $fastqDir/sample_11_2.fq.gz 

kallisto quant -i  $idx -o $kallistoDir/sample12 --threads 12 \
$fastqDir/sample_12_1.fq.gz $fastqDir/sample_12_2.fq.gz 

# output files, equivalence classes:
mkdir ec
kallistoDir_ec="/6vs6_sim/aligned_reads/kallisto/ec"

# pseudo (to obtain the equivalence classes):
kallisto pseudo -i  $idx -o $kallistoDir_ec/sample1 \
$fastqDir/sample_1_1.fq.gz $fastqDir/sample_1_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample2 \
$fastqDir/sample_2_1.fq.gz $fastqDir/sample_2_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample3 \
$fastqDir/sample_3_1.fq.gz $fastqDir/sample_3_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample4 \
$fastqDir/sample_4_1.fq.gz $fastqDir/sample_4_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample5 \
$fastqDir/sample_5_1.fq.gz $fastqDir/sample_5_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample6 \
$fastqDir/sample_6_1.fq.gz $fastqDir/sample_6_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample7 \
$fastqDir/sample_7_1.fq.gz $fastqDir/sample_7_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample8 \
$fastqDir/sample_8_1.fq.gz $fastqDir/sample_8_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample9 \
$fastqDir/sample_9_1.fq.gz $fastqDir/sample_9_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample10 \
$fastqDir/sample_10_1.fq.gz $fastqDir/sample_10_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample11 \
$fastqDir/sample_11_1.fq.gz $fastqDir/sample_11_2.fq.gz 

kallisto pseudo -i  $idx -o $kallistoDir_ec/sample12 \
$fastqDir/sample_12_1.fq.gz $fastqDir/sample_12_2.fq.gz 
