#!/bin/bash 
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" 
#$ -pe threads 2 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/proc_b2_1.out -j y 
export BOWTIE_PATH=/mnt/wigclust1/data/safe/kostic/bowtie-1.2.1.1
export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data
export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts
export SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19
cd $DATA_PATH 
$SAM_PATH/samtools view -bhS aligned_b2_1.sam | $SAM_PATH/samtools view -bh -F 4 - | $SAM_PATH/samtools sort -n - b2_1_sorted  
python $PYTHON_PATH/add_varietal_tag_paired.py <( gunzip -c $DATA_PATH/b2/s_1_1_sequence.b1.txt.gz ) <( gunzip -c $DATA_PATH/b2/s_1_2_sequence.b1.txt.gz ) <( $SAM_PATH/samtools view -h b2_1_sorted.bam ) | $SAM_PATH/samtools view -Sbh - > b2_1_tagged.bam
