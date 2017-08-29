#!/bin/bash

#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]"
#$ -pe threads 2 
#$ -l vf=8G 
#$ -o /mnt/wigclust1/data/safe/kostic/output/output_paired.out -j y


export PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts

export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data

export SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19

cd $DATA_PATH

$SAM_PATH/samtools view -bS uniquename.sam > uniquename.bam


79)
Convert the output to .bam file format:
/filepath/samtools-0.1.16/samtools view -Sb -o /filepath/SRR054616.bam /filepath/SRR054616.sam


80) (//sort so that read /1 and /2 are one after another)
Sort the .bam file:
/filepath/samtools-0.1.16/samtools sort /filepath/SRR054616.bam /filepath/SRR054616.sorted



add varietal tag

python /mnt/wigclust5/data/safe/kendall/nla3_88/add.varietal.tag01.bt1.bypair.py <(gunzip -c /mnt/wigclust5/data/safe/kendall/nla3_88/b2/s_1_1_sequence.b2.txt.gz) <(gunzip -c /mnt/wigclust5/data/safe/kendall/nla3_88/b2/s_1_2_sequence.b2.txt.gz) <(/mnt/wigclust5/data/safe/kendall/samtools-0.1.19/samtools view -h /mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014296.r1.ucsc.hg19dm6.bt1.sorted.bam) 1 | /mnt/wigclust5/data/safe/kendall/samtools-0.1.19/samtools view -Sbh - > /mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014296.r1.ucsc.hg19dm6.bt1.sorted.vt.bypair.bam

(create an index)
/mnt/wigclust5/data/safe/kendall/samtools-0.1.19/samtools index /mnt/wigclust5/data/safe/kendall/nla3_88/b2/BEI014296.r1.ucsc.hg19dm6.bt1.sorted.vt.bypair.bam





81)
Remove reads likely to be PCR duplicates:
/filepath/samtools-0.1.16/samtools rmdup -s /filepath/SRR054616.sorted.bam /filepath/SRR054616.rmdup.bam




82)
Create a .bam file index:
/filepath/samtools-0.1.16/samtools index /filepath/SRR054616.rmdup.bam
83)
Create a .sam file from the sorted .bam file with duplicates removed:
/filepath/samtools-0.1.16/samtools view -o /filepath/SRR054616.rmdup.sam /filepath/SRR054616.rmdup.bam



/mnt/wigclust5/data/safe/kendall/bowtie-0.12.7/bowtie -S -t -n 2 -e 170 --chunkmbs 256 -3 0 -5 16 -m 1 --best --strata hg19dm6 <(gunzip -c /mnt/wigclust5/data/safe/kendall/nla3_88/b5/s_1_1_sequence.b10.txt.gz | python /mnt/wigclust5/data/safe/kendall/pythonpgms/adapter.clip02.min58.py GATCGGAAGAGCGG) 2> /mnt/wigclust5/data/safe/kendall/nla3_88/b5/BEI014304.r1.bt1.report.txt | /mnt/wigclust5/data/safe/kendall/samtools-0.1.16/samtools view -Sbho - - | /mnt/wigclust5/data/safe/kendall/samtools-0.1.16/samtools sort - /mnt/wigclust5/data/safe/kendall/nla3_88/b5/BEI014304.r1.ucsc.hg19dm6.bt1.sorted

/mnt/wigclust5/data/safe/kendall/samtools-0.1.16/samtools index /mnt/wigclust5/data/safe/kendall/nla3_88/b5/BEI014304.r1.ucsc.hg19dm6.bt1.sorted.bam

/mnt/wigclust5/data/safe/kendall/bowtie-0.12.7/bowtie -S -t -n 2 -e 170 --chunkmbs 256 -3 0 -5 16 -m 1 --best --strata hg19dm6 <(gunzip -c /mnt/wigclust5/data/safe/kendall/nla3_88/b5/s_1_2_sequence.b10.txt.gz | python /mnt/wigclust5/data/safe/kendall/pythonpgms/adapter.clip02.min58.py GATCGGAAGAGCGG) 2> /mnt/wigclust5/data/safe/kendall/nla3_88/b5/BEI014304.r2.bt1.report.txt | /mnt/wigclust5/data/safe/kendall/samtools-0.1.16/samtools view -Sbho - - | /mnt/wigclust5/data/safe/kendall/samtools-0.1.16/samtools sort - /mnt/wigclust5/data/safe/kendall/nla3_88/b5/BEI014304.r2.ucsc.hg19dm6.bt1.sorted

/mnt/wigclust5/data/safe/kendall/samtools-0.1.16/samtools index /mnt/wigclust5/data/safe/kendall/nla3_88/b5/BEI014304.r2.ucsc.hg19dm6.bt1.sorted.bam



