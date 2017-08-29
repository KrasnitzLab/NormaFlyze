#!/bin/bash 
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" 
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/count_b26_10GC.out -j y 
export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data_2
export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts
export SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19
export GENOME_PATH=/mnt/wigclust1/data/safe/kostic/hybrid_genome
cd $DATA_PATH 
$SAM_PATH/samtools sort b26_10_tagged.bam b26_10_sorted
$PYTHON_PATH/compute_frag_gc_content.py <( $SAM_PATH/samtools view b26_10_sorted.bam ) $GENOME_PATH b26_10_GCadded.sam
python $PYTHON_PATH/count_and_dedup_GC.py <( $SAM_PATH/samtools view b26_10_tagged.bam )  b26_10_varbinGC_count.txt b26_10_varbinGC_stats.txt
