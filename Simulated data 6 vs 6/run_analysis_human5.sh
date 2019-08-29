#!/bin/bash

# The pipieline for our simulation was modified from: https://github.com/markrobinsonuzh/diff_splice_paper
# In partcular, we considered the human simulation with DGE and introduced two changes:
# 1) we increased the sample size of each group from 3 to 6 (we now have a 6 vs 6 group comparison);
# 2) we changed the definition of DTU: we randomly permute the relative abundance of the top four expressed transcripts,
# while Soneson et al. invert the relative abundance of the top two expressed transcripts.

## Define paths to software and reference files
BASEDIR=/DTU_Sim_6vs6
SRATOOLKIT=$BASEDIR/software/sratoolkit.2.5.0-1-ubuntu64
REFERENCEDIR=$BASEDIR/hsapiens/reference_files
RSEM=$BASEDIR/software/rsem-1.2.21
ASTALAVISTA=$BASEDIR/software/astalavista-3.2/bin
FIGDIR=$BASEDIR/hsapiens/figures
RCODEGEN=$BASEDIR/software/Rcode
ROUT=$BASEDIR/hsapiens/Rout

# output of RSEM files and RNA-seq reads:
RSEM_FILES=$BASEDIR/rsem_files
READS=$BASEDIR/reads

## ------------------------- INPUT PREPARATION ----------------------------- ##

## The basis for the simulation is generated from a sample downloaded from the SRA. 
## This sample was sequenced with an Illumina HiSeq, with a paired-end protocol,
## and with a read length of 101 bp.
## Extract the fastq files from the sra archive
$SRATOOLKIT/fastq-dump --split-files \
-Q 33 -O $REFERENCEDIR $REFERENCEDIR/SRR493366.sra

## Extract only lines corresponding to primary assembly chromosomes from gtf file
grep -v "^H" $REFERENCEDIR/Homo_sapiens.GRCh37.71.gtf \
> $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.gtf

## Extract only lines corresponding to protein coding genes from gtf file
grep "protein_coding" $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.gtf \
> $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf

## Prepare the reference files for RSEM
$RSEM/rsem-prepare-reference \
--gtf $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf --bowtie2 \
$REFERENCEDIR/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa \
$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71

## Estimate the model file with RSEM
$RSEM/rsem-calculate-expression \
--paired-end --bowtie2 --seed 123 \
$REFERENCEDIR/SRR493366_1.fastq $REFERENCEDIR/SRR493366_2.fastq \
$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71 \
$REFERENCEDIR/rsem_model/SRR493366

## Plot some characteristics of the estimated model
$RSEM/rsem-plot-model $REFERENCEDIR/rsem_model/SRR493366 $FIGDIR/RSEM_model.pdf

## Modify the quality score distribution in the RSEM model file so that the probability 
## of quality score 2 is 0. Otherwise, the quality scores of the simulated data may be very low.
## Also make sure that the transition probabilities into this state are 0.
## ->> $REFERENCEDIR/rsem_model/SRR493366.stat/SRR493366.highQ.model

## Estimate and plot the isoform percentage distributions from the RSEM results
R CMD BATCH --no-restore --no-save "--args referencefile='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' outdir='$FIGDIR'" $RCODEGEN/isopct_distribution.R $ROUT/isopct_distribution_human.Rout

## Run ASTALAVISTA to classify splicing events
$ASTALAVISTA/astalavista -t asta \
-i $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf -e [ASE,ASI,DSP,VST]
gunzip $REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_sorted.gtf_astalavista.gtf.gz

## -------------------------- DATA SIMULATION ------------------------------ ##

## Estimate the simulation parameters for the individual samples
R CMD BATCH --no-restore --no-save "--args path_to_generate_rsem_files='$RCODEGEN/generate_rsem_files_function.R' seed=123 isoform_results_file='$REFERENCEDIR/rsem_model/SRR493366.isoforms.results' nbr_per_group=6 meandisp.file='$REFERENCEDIR/Pickrell.Cheung.Mu.Phi.Estimates.rds' outdirbase='$BASEDIR' librarysize=40000000 keepchr=NULL nbr_diff_spliced=1000 nbr_diff_expr=1000 fold_changes='expon'" $RCODEGEN/generate_rsem_files_human_run.R $ROUT/generate_rsem_files_human_run_de.Rout

## Generate truth files
R CMD BATCH --no-restore --no-save "--args path_to_generate_truth_file='$RCODEGEN/generate_truth_table_function.R' path_to_final_summary='$BASEDIR/simulation_details.txt' out.file='$BASEDIR/truth_human_non_null.txt' astalavista.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_sorted.gtf_astalavista.gtf' gtf.file='$REFERENCEDIR/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf' flattened.gtf.file=NULL missing.annot.file=NULL" $RCODEGEN/generate_truth_table_run.R $ROUT/generate_truth_table_human_nonnull_de.Rout

## Simulate reads for 12 samples (40 million pairs/sample),
## non-null DTU situation, with differential expression
for n in 1 2 3 4 5 6 7 8 9 10 11 12
do
$RSEM/rsem-simulate-reads \
$REFERENCEDIR/rsem_reference/Homo_sapiens.GRCh37.71 \
$REFERENCEDIR/rsem_model/SRR493366.stat/SRR493366.highQ.model \
$BASEDIR/rsem_files/sample${n}.txt \
0.05 40000000 $READS/sample_${n} \
--seed 123
done
