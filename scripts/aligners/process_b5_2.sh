#!/bin/bash 
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" 
#$ -pe threads 2 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/proc_b5_2.out -j y 
export BOWTIE_PATH=/mnt/wigclust1/data/safe/kostic/bowtie-1.2.1.1
export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data
export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts
export SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19
cd $DATA_PATH 
$SAM_PATH/samtools view -bhS aligned_b5_2.sam | $SAM_PATH/samtools view -bh -F 4 - | $SAM_PATH/samtools sort -n - b5_2_sorted  
python $PYTHON_PATH/add_varietal_tag_paired.py <( gunzip -c $DATA_PATH/b5/s_1_1_sequence.b2.txt.gz ) <( gunzip -c $DATA_PATH/b5/s_1_2_sequence.b2.txt.gz ) <( $SAM_PATH/samtools view -h b5_2_sorted.bam ) | $SAM_PATH/samtools view -Sbh - > b5_2_tagged.bam
