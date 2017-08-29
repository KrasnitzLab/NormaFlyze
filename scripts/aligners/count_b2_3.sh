#!/bin/bash 
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" 
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/count_b2_3.out -j y 
export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data
export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts
export SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19
cd $DATA_PATH 
python $PYTHON_PATH/count_and_dedup.py <($SAM_PATH/samtools view b2_3_tagged.bam) b2_3_varbin_count.txt b2_3_varbin_stats.txt
