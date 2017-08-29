#!/bin/bash

#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/outputfile_zoning.out -j y

export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts

export DATA_PATH=/mnt/wigclust1/data/safe/kostic/bin_mapping

export GENOME_PATH=/mnt/wigclust1/data/safe/kostic/hybrid_genome

export SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19

cd $DATA_PATH

#get chromosome sizes
#input = path to the hybrid genome; the desired name of the output file
python $PYTHON_PATH/chrom_sizes.py $GENOME_PATH chrom_sizes_hg_dm_combined_consecutive.txt

#convert to bam and sort the alignments
$SAM_PATH/samtools view -bhS paired_7_13.sam | $SAM_PATH/samtools sort -n - paired_sorted_7_13

#only keep alignments that fall within your data range as determined by ___
$SAM_PATH/samtools view paired_sorted_7_13.bam | awk '($9 >= 150 && $9 <= 500) || ($9 <= -150 && $9 >= -500)' - > ps_data_lim_150_500.sam

#output = .bed file containing all alignments on each chromosome; .txt file containing the number of these alignments for each chromosome
python $PYTHON_PATH/get_goodzones.py ps_data_lim_150_500.sam

#sort in lexicographic order
sort -k 1,1 -k 2,2n hybrid_unique_mappers_07_28.bed > sorted_mappers_150_500.bed

#arg 1 = hybrid_num_mappable.txt ; arg 2 = sorted_mappers.bed ; arg 3 = chromosome sizes file ; arg 4 = the location for the output file
#arg 2 must be in lex order (chr1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,2l,2r,3,3l,3r,4,4_dm,etc)
python $PYTHON_PATH/sep_bin_boundaries.py $DATA_PATH/hybrid_num_mappable_07_28.txt $DATA_PATH/sorted_mappers_150_500.bed $DATA_PATH/chrom_sizes_hg_dm_combined_consecutive.txt $DATA_PATH

#sort by absolute position - first all hg then all dm
sort -k 3,3n hybrid_bin_boundaries_07_28.txt > hybrid_bin_boundaries_sorted_150_500.txt
