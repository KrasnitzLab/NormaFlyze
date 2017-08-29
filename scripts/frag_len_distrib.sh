#!/bin/bash

#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/len_distrib.out -j y

export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts

export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data_2

cd $DATA_PATH

python $PYTHON_PATH/frag_len_distrib.py $DATA_PATH
# sort -n $DATA_PATH/data_frag_len_distrib.txt -o $DATA_PATH/data_frag_len_distrib.txt