# the reference transcriptome:
fasta="Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa"

# input files
fastqDir="/6vs6_sim/reads"

# output files
$salmonDir="/6vs6_sim/aligned_reads/Salmon"

# input files
fastqDir="/home/stiber/data_SIM5_Simone/reads"

# output files, in Taupo
salmonDir="/home/stiber/SIM5_Simone/FILTERED/Salmon"

cd $salmonDir

# index name
idx="Homo_sapiens.GRCh37.71.SELECTED.cdna.all.fa.idx_Salmon0.8"

# Build index
salmon index -i $idx -t $fasta -p 12 --type quasi -k 31

# quant
salmon quant -i  $idx -l A -1 $fastqDir/sample_1_1.fq.gz -2 $fastqDir/sample_1_2.fq.gz \
-p 12 -o $salmonDir/sample1 --dumpEq

salmon quant -i  $idx -l A -1 $fastqDir/sample_2_1.fq.gz -2 $fastqDir/sample_2_2.fq.gz \
  -p 12 -o $salmonDir/sample2 --dumpEq

salmon quant -i  $idx -l A -1 $fastqDir/sample_3_1.fq.gz -2 $fastqDir/sample_3_2.fq.gz \
  -p 12 -o $salmonDir/sample3 --dumpEq

salmon quant -i  $idx -l A -1 $fastqDir/sample_4_1.fq.gz -2 $fastqDir/sample_4_2.fq.gz \
  -p 12 -o $salmonDir/sample4 --dumpEq
  
salmon quant -i  $idx -l A -1 $fastqDir/sample_5_1.fq.gz -2 $fastqDir/sample_5_2.fq.gz \
  -p 12 -o $salmonDir/sample5 --dumpEq
  
salmon quant -i  $idx -l A -1 $fastqDir/sample_6_1.fq.gz -2 $fastqDir/sample_6_2.fq.gz \
  -p 12 -o $salmonDir/sample6 --dumpEq
  
salmon quant -i  $idx -l A -1 $fastqDir/sample_7_1.fq.gz -2 $fastqDir/sample_7_2.fq.gz \
-p 12 -o $salmonDir/sample7 --dumpEq

salmon quant -i  $idx -l A -1 $fastqDir/sample_8_1.fq.gz -2 $fastqDir/sample_8_2.fq.gz \
  -p 12 -o $salmonDir/sample8 --dumpEq
  
salmon quant -i  $idx -l A -1 $fastqDir/sample_9_1.fq.gz -2 $fastqDir/sample_9_2.fq.gz \
  -p 12 -o $salmonDir/sample9 --dumpEq
  
salmon quant -i  $idx -l A -1 $fastqDir/sample_10_1.fq.gz -2 $fastqDir/sample_10_2.fq.gz \
  -p 12 -o $salmonDir/sample10 --dumpEq
  
salmon quant -i  $idx -l A -1 $fastqDir/sample_11_1.fq.gz -2 $fastqDir/sample_11_2.fq.gz \
  -p 12 -o $salmonDir/sample11 --dumpEq
  
salmon quant -i  $idx -l A -1 $fastqDir/sample_12_1.fq.gz -2 $fastqDir/sample_12_2.fq.gz \
  -p 12 -o $salmonDir/sample12 --dumpEq
