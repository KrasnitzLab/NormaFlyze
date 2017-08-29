#!/bin/bash

#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -pe threads 1 
#$ -l vf=16G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/outputfile.out -j y

#path to Bowtie software
export BOWTIE_PATH=/mnt/wigclust1/data/safe/kostic/bowtie-1.2.1.1
#path to folder where python scripts are held
export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts
#path to the hybrid genome .fa files
export DATA_PATH=/mnt/wigclust1/data/safe/kostic/hybrid_genome_2

cd $DATA_PATH

#extract all the human chromosomes
tar -zxvf chromFa.tar.gz

#mask the pseudoautosomal region on the hg19 Y chromosome
#arg1 = infile, arg2 = outfile
python $PYTHON_PATH/hg19_mask_Y.py chrY.fa chrY_masked.fa
mv chrY_masked.fa chrY.fa

#remove "hap" chromosomes from the hg19 genome
rm *hap*

#add _dm to all fly chromosomes and separate into files
python $PYTHON_PATH/sep_drosoph_genome.py <(gunzip -c dm6.fa.gz) $DATA_PATH

cat *.fa > whole_hybrid_genome.fa

#build the bowtie index
$BOWTIE_PATH/bowtie-build -f whole_hybrid_genome.fa $BOWTIE_PATH/indexes/hybrid_indices_2