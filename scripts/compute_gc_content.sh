#!/bin/bash

#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/gc.out -j y

export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts

export DATA_PATH=/mnt/wigclust1/data/safe/kostic/bin_mapping

export GENOME_PATH=/mnt/wigclust1/data/safe/kostic/hybrid_genome

cd $DATA_PATH

python $PYTHON_PATH/compute_gc_content.py  $GENOME_PATH $DATA_PATH/hg19_bin_boundaries_sorted_125_600.txt
#hybrid_bin_boundaries_sorted_125_600.txt

python $PYTHON_PATH/compute_len_distrib.py $DATA_PATH/hg19_bin_boundaries_sorted_125_600.txt $DATA_PATH/sorted_mappers_hg19_125_600.bed
#sorted_mappers_125_600.bed
#hybrid_bin_boundaries_sorted_125_600.txt

mv hybrid_genome_GCcontent10_01.txt hg19_GC.txt
mv hybrid_bin_fragment_stats10_01.txt hg19_LEN.txt

#range125_600

# original names of output files:
# hybrid_genome_GCcontent" + time.strftime("%m_%d") + ".txt"
# "hybrid_bin_fragment_stats" + time.strftime("%m_%d") + ".txt"