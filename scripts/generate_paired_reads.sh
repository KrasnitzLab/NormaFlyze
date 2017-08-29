#!/bin/bash

#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/output_paired_reads.out -j y

export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts

export DATA_PATH=/mnt/wigclust1/data/safe/kostic/hybrid_genome

export RESULTS_PATH=/mnt/wigclust1/data/safe/kostic/cut_sites

cd $RESULTS_PATH

#arg1 = the path containing the hybrid genome .fa files
#arg2 = the desired legnth of the reads
#arg3 = the restriction enzyme (Nla3) cutsite
python $PYTHON_PATH/generate_paired_reads.py $DATA_PATH 150 "CATG"
#output in RESULTS_PATH = 2 fastq files representing paired end reads of length 150 starting and ending at CATG sites