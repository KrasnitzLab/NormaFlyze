#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import os

def main():

	project_dir = "/mnt/wigclust1/data/safe/kostic/SNS_data_2"
	bash_dir = "/mnt/wigclust1/data/safe/kostic/scripts/aligners"

	for f in os.listdir(project_dir):
		if os.path.isfile(os.path.join(project_dir,f)) and ("tagged.bam" in f):

			fname = "_".join(f.split("_")[0:2])
			qsubFile = bash_dir + "/count_" + fname + ".sh"
			QSUB = open(qsubFile, "w")

			header = ('#!/bin/bash \n#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" '
							'\n#$ -pe threads 1 '
							'\n#$ -l vf=8G '
							'\n#$ -o /mnt/wigclust1/data/safe/kostic/output/count_''' + fname + '.out -j y \n')
			QSUB.write(header)

			paths = ('export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data_2' 
							'\nexport PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts'
							'\nexport SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19'
							'\ncd $DATA_PATH \n')
			QSUB.write(paths)

			#in count_and_dedup.py, change the input bin_boundaries file and the chromosome length file to the appropriate ones
			command = ('python $PYTHON_PATH/count_and_dedup.py <($SAM_PATH/samtools view ' + f + ')' 
						' ' + fname + '_varbin_count.txt ' + fname + '_varbin_stats.txt\n')
			QSUB.write(command)

			QSUB.close()

			time.sleep(1)
			thisCommand = "chmod 755 " + qsubFile
			os.system(thisCommand)
			thisCommand = "qsub " + qsubFile
			os.system(thisCommand)


if __name__ == "__main__":
	main()