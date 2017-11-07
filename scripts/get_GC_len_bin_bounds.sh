#!/bin/bash

#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/outputfile_GCzoning.out -j y

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

#only keep alignments that fall within your data range as determined by fragment length distribution
$SAM_PATH/samtools view paired_sorted_7_13.bam | awk '($9 >= 125 && $9 <= 600) || ($9 <= -125 && $9 >= -600)' - > ps_data_lim_125_600.sam

#output = .bed file containing all alignments on each chromosome; .txt file containing the number of these alignments for each chromosome
python $PYTHON_PATH/get_goodzones.py ps_data_lim_125_600.sam

#get GC content of the simulated fragments
python $PYTHON_PATH/compute_quantile_gc_content.py hybrid_unique_mappers_08_08.bed $GENOME_PATH fragment_GC_125_600.txt 

#arg 1 = chromosome sizes file ; arg2 = the fragment GC file; arg 3 = the output file
python $PYTHON_PATH/sep_GC_len_bin_bounds.py $DATA_PATH/chrom_sizes_hg_dm_combined_consecutive.txt $DATA_PATH/fragment_GC_125_600.txt GC_bin_bounds_125_600.txt
