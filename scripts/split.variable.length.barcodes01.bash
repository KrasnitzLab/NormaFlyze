#!/bin/bash

#$ -q "all.q@wigclust2[3]"
#$ -pe threads 1 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/splitBarcodes.out -j y



#####   Run this in some directory.  It will create 4 subdirs b15, b26, b37, b48.  Also, a couple other subdirs and files.
#####   Within the 4 subdirs run the other barcode split program.  

cd /mnt/wigclust1/data/safe/kostic/SNS_data_SKBR3/

/mnt/wigclust15/data/safe/leey/ToWhom/zihua/code/sfbc.differnet.barcodes.length/trimBarcodeFragments -nThreads=8 -inputFiles=/mnt/wigclust1/data/safe/kostic/SNS_data_SKBR3/1_S1_L001_R1_001.fastq.gz,/mnt/wigclust1/data/safe/kostic/SNS_data_SKBR3/1_S1_L001_R2_001.fastq.gz "b15:GTCAGCT;TCTCACCT;GGTATTCGT;AATTGATGCT;AGGACCT;TAACAGCT;TGAGTTGGT;CAGTGAACCT" "b26:CACTCGT;CTAGACGT;ACGACAGCT;TGGCATACGT;GAGTGCT;GTTATCCT;ACAGCAGCT;TGAACATGCT" "b37:TCAAGCT;ATTGACGT;CAGATTCCT;TGTCGATAGT;AGGAGGT;AGTCTGGT;TCTGTACCT;CGACATACGT" "b48:CGTTCGT;GCATAGGT;AAGGATGCT;GTCGTAAGCT;CTATGGT;GACTTCCT;GGCACTGCT;GACTGTTCCT" &





# /mnt/wigclust15/data/safe/leey/ToWhom/zihua/code/sfbc.differnet.barcodes.length/trimBarcodeFragments -nThreads=8 
# -inputFiles=/mnt/wigclust5/data/safe/kendall/nla3_88/b1/r1.fastq.gz,/mnt/wigclust5/data/safe/kendall/nla3_88/b1/r2.fastq.gz 
# "b1:GTCAGCT;TCTCACCT;GGTATTCGT;AATTGATGCT" "b2:CACTCGT;CTAGACGT;ACGACAGCT;TGGCATACGT" "b3:TCAAGCT;ATTGACGT;CAGATTCCT;TGTCGATAGT" "b4:CGTTCGT;GCATAGGT;AAGGATGCT;GTCGTAAGCT" "b5:AGGACCT;TAACAGCT;TGAGTTGGT;CAGTGAACCT" "b6:GAGTGCT;GTTATCCT;ACAGCAGCT;TGAACATGCT" "b7:AGGAGGT;AGTCTGGT;TCTGTACCT;CGACATACGT" "b8:CTATGGT;GACTTCCT;GGCACTGCT;GACTGTTCCT" &

