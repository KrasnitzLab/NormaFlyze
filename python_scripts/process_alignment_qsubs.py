#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import os

def main():

	project_dir = "/mnt/wigclust1/data/safe/kostic/SNS_data_2"
	bash_dir = "/mnt/wigclust1/data/safe/kostic/scripts/aligners"

	for file in os.listdir(project_dir):
		if os.path.isfile(os.path.join(project_dir,file)) and ("aligned" in file) and (".sam" in file):

			subdir = file.split("_")[1]
			fastq = (file.split("_")[2]).split(".")[0]
			fname = subdir + "_" + fastq
			qsubFile = bash_dir + "/process_" + fname + ".sh"
			QSUB = open(qsubFile, "w")

			header = ('#!/bin/bash \n#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" '
							'\n#$ -pe threads 2 '
							'\n#$ -l vf=8G '
							'\n#$ -o /mnt/wigclust1/data/safe/kostic/output/proc_''' + fname + '.out -j y \n')
			QSUB.write(header)

			paths = ('export BOWTIE_PATH=/mnt/wigclust1/data/safe/kostic/bowtie-1.2.1.1' 
							'\nexport DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data_2' 
							'\nexport PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts'
							'\nexport SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19'
							'\ncd $DATA_PATH \n')
			QSUB.write(paths)

			#extract only the correctly aligned pairs, sort them in a bam file, add varietal tags
			command = ('$SAM_PATH/samtools view -bhS ' + file + ' | $SAM_PATH/samtools view -bh -F 4 - | $SAM_PATH/samtools sort -n - ' + fname + '_sorted  '
						'\npython $PYTHON_PATH/add_varietal_tag_paired.py ' 
						'<( gunzip -c $DATA_PATH/' + subdir + '/s_1_1_sequence.b' + fastq + '.txt.gz ) '
						'<( gunzip -c $DATA_PATH/' + subdir + '/s_1_2_sequence.b' + fastq + '.txt.gz ) '
						'<( $SAM_PATH/samtools view -h ' + fname + '_sorted.bam ) | $SAM_PATH/samtools view -Sbh - > ' + fname + '_tagged.bam\n' 
						)
			QSUB.write(command)

			# command = ('$SAM_PATH/samtools sort ' + fname + '_tagged.bam tagged_out ' 
			# 			'\n$SAM_PATH/samtools index tagged_out.bam')
			# QSUB.write(command)

			QSUB.close()

			time.sleep(1)
			thisCommand = "chmod 755 " + qsubFile
			os.system(thisCommand)
			thisCommand = "qsub " + qsubFile
			os.system(thisCommand)


if __name__ == "__main__":
	main()


