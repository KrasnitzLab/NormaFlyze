#!/bin/bash

#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/outputfile.out -j y

export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts

export DATA_PATH=/mnt/wigclust1/data/safe/kostic/hybrid_genome

export RESULTS_PATH=/mnt/wigclust1/data/safe/kostic/cut_sites

cd $RESULTS_PATH

# locator outputs just the seqs, new one on each line
# python $PYTHON_PATH/cut_site_locator.py $DATA_PATH 120 "CATG"

# finder outputs the FASTQ formatted sequences
python $PYTHON_PATH/cut_site_finder.py $DATA_PATH 120 "CATG"