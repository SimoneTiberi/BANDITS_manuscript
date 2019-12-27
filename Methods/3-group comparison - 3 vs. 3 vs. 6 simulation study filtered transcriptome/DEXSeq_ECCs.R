export LC_ALL=C

# upgrade pip locally (user):
# pip install --upgrade pip --user

# update python pandas:
# python3 -m pip install --upgrade pandas --user

# git clone git@github.com:Oshlack/ec-dtu-pipe.git
ec="/software/ec-dtu-paper-master"

python3 $ec/create_salmon_ec_count_matrix.py \
/6vs6_sim/aligned_reads/sample1/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample2/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample3/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample4/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample5/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample6/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample7/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample8/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample9/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample10/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample11/aux_info/eq_classes.txt \
/6vs6_sim/aligned_reads/sample12/aux_info/eq_classes.txt \
sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10,sample11,sample12 \
ec_matrix.txt


# in R:
R

set.seed(61217)

# Load truth table:
truth_dir="DTU_Sim_6vs6/simulation_details.txt"
truth = data.table::fread(truth_dir, header = TRUE, sep="\t" )

ensg <- data.frame(gene_name = truth$gene_id, gene_id = truth$gene_id, transcript = truth$transcript_id,
                   stringsAsFactors = FALSE)
# by default, data.frame converts strings into factors!!!

library(dplyr)
library(data.table)
library(DEXSeq)
library(DRIMSeq)

# We are now ready to load our data:
ecm <- read.delim('ec_matrix.txt', sep='\t')
ecm <- inner_join(ecm, ensg, by='transcript')

### TRANSCRIPTS FILTERING:
load("transcripts_to_keep.RData")
ecm = ecm[ecm$transcript %in% transcripts_to_keep, ]
# ~51% of the ECs are kept
###

# Next we remove equivalence classes that map to multiple different genes:
multi_ecs <- data.table(ecm)[, length(unique(gene_id)), keyby=ec_names]
multi_ecs <- multi_ecs[multi_ecs$V1 > 1]
multi_ecs <- multi_ecs$ec_names
df <- ecm[!ecm$ec_names %in% multi_ecs,]

# We now define our sample groups:
group <- c( rep(0, 3), rep(1, 3), rep(2, 6) )
samples <- paste0("sample", 1:12)
sampleTable <- data.frame(sample_id = samples,
                          condition = group)

# Our EC data frame still potentially contains repeat entries as multiple transcripts mapping to the same EC will be divided into multiple rows. Therefore, we will remove that transcript ID column, and consider all distinct equivalence classes per gene. As we will use DRIMSeq to prepare and filter our data, we will also change create a 'feature_id' column which holds the equivalence class IDs:
df$feature_id <- df$ec_names
df <- distinct(df[,c(samples, 'gene_id', 'feature_id')])

# Our data is now ready to load into a DRIMSeq object.
d <- dmDSdata(counts = df, samples = sampleTable)

# Optionally, we may filter our data further by using some standard transcript and gene-level filters:
d <- dmFilter(d, min_samps_feature_expr=0,
              min_feature_expr=0,
              min_samps_feature_prop=0,
              min_samps_gene_expr=0,
              min_gene_expr=10)

# Our data is now ready to load into DEXSeq to perform differential transcript usage analysis:
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))

dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)
# We will use MulticoreParam to allow DEXSeq to take advantage of multiple for cores when performing analysis. This greatly speeds up analysis. First, we estimate the library sizes, then estimate dispersions, which we can visualise:
BPPARAM=MulticoreParam(workers=12)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM)

dxd <- testForDEU(dxd, BPPARAM=BPPARAM)
res <- DEXSeqResults(dxd)

qval = perGeneQValue(res)
res_gene = data.frame(gene = names(qval), qval)

save(res, res_gene, file="DEXSeq_ECs results.RData")
