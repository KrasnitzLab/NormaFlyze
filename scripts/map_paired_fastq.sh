#!/bin/bash

#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -pe threads 2
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/output_paired_real.out -j y

export BOWTIE_PATH=/mnt/wigclust1/data/safe/kostic/bowtie-1.2.1.1

export DATA_PATH=/mnt/wigclust1/data/safe/kostic/cut_sites

export RESULTS_PATH=/mnt/wigclust1/data/safe/kostic/bin_mapping

export BOWTIE_INDEXES=$BOWTIE_PATH/indexes

cd $DATA_PATH

# see manual for more: http://bowtie-bio.sourceforge.net/manual.shtml#the--n-alignment-mode
# -q for fastq formatted read file
# -S causes the alignments to be printed in SAM format
# -t outputs time for each phase
# -p <int> launches parallel search threads
# -n <int1> -e <int2> allows int1 mismatches in the seed and a sum of int2 for quality scores of mismatches at any other location
# -m 1 for unique alignments

#paired-end read options
# -1 and -2 for paired end read files
# -X <int> gives the maximum insert size for paired end read alignment (default 250)
# -I <int> gives the minimum insert size for paired end read alignment (default 0)
# --allow-contain allows mapping of paired end reads that overlap

#bowtie [options] <index base name> <-1 fastq1 -2 fastq2> <output file>
$BOWTIE_PATH/bowtie -q -S -t -p 2 -n 2 -e 200 --chunkmbs 256 -m 1 -X 800 --best --strata --allow-contain hybrid_index -1 $DATA_PATH/CATG_paired_150bp_07_13_fastq_1.fq -2 $DATA_PATH/CATG_paired_150bp_07_13_fastq_2.fq $RESULTS_PATH/paired_7_13.sam
