#!/bin/bash 
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" 
#$ -pe threads 2 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/b7_9.out -j y 
export BOWTIE_PATH=/mnt/wigclust1/data/safe/kostic/bowtie-1.2.1.1
export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data/b7
export RESULTS_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data/
export BOWTIE_INDEXES=$BOWTIE_PATH/indexes
export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts
cd $DATA_PATH 
gunzip -c s_1_1_sequence.b9.txt.gz | python $PYTHON_PATH/remove_seq_adapter.py GATCGGAAGAGCGG > s_1_1_sequence.b9.txt
gunzip -c s_1_2_sequence.b9.txt.gz | python $PYTHON_PATH/remove_seq_adapter.py GATCGGAAGAGCGG > s_1_2_sequence.b9.txt
$BOWTIE_PATH/bowtie -q -S -t -n 2 -e 200 --chunkmbs 256 -3 0 -5 16 -m 1 -X 800 --allow-contain --best --strata hybrid_index -1 s_1_1_sequence.b9.txt -2 s_1_2_sequence.b9.txt $RESULTS_PATH/aligned_b7_9.sam
rm s_1_1_sequence.b9.txt
rm s_1_2_sequence.b9.txt
