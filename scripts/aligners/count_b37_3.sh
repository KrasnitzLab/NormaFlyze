#!/bin/bash 
#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" 
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/count_b37_3.out -j y 
export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data_2
export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts
export SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19
cd $DATA_PATH 
awk '{ if($9 > 0) print $3,"	",$4,"	",$4+$9 }' <($SAM_PATH/samtools view b37_3_tagged.bam) | python $PYTHON_PATH/compute_frag_gc_content.py - $GENOME_PATH b37_3_GCcontent.txtpython $PYTHON_PATH/count_and_dedup_GC.py <($SAM_PATH/samtools view b37_3_tagged.bam)b37_3_GCcontent.txt b37_3_varbinGC_count.txt b37_3_varbinGC_stats.txt
