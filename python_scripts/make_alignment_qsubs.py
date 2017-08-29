#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import os

def main():

	project_dir = "/mnt/wigclust1/data/safe/kostic/SNS_data_SKBR3"
	bash_dir = "/mnt/wigclust1/data/safe/kostic/scripts/aligners"

	for subdir, dirs, files in os.walk(project_dir):
		for file in files:
			if "s_1_1" in file:
				pair = file.split("_")
				pair = "_".join(pair[0:2]) + "_2_" + str(pair[3])

				sub =  subdir.split("/")[-1]
				tag = file.split(".")[1][1:]
				qsubFile = bash_dir + "/align_" + sub + "_" + tag + ".sh"
				QSUB = open(qsubFile, "w")

				samfile = "aligned_" + sub + "_" + tag + ".sam"

				header = ('#!/bin/bash \n#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" '
									'\n#$ -pe threads 2 '
									'\n#$ -l vf=8G '
									'\n#$ -o /mnt/wigclust1/data/safe/kostic/output/''' + sub + '_' + tag + '.out -j y \n')
				QSUB.write(header)

				paths = ('export BOWTIE_PATH=/mnt/wigclust1/data/safe/kostic/bowtie-1.2.1.1' 
							'\nexport DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data_SKBR3/' + sub + ''
							'\nexport RESULTS_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data_SKBR3/' 
							'\nexport BOWTIE_INDEXES=$BOWTIE_PATH/indexes'
							'\nexport PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts'
							'\ncd $DATA_PATH \n')
				QSUB.write(paths)

				unzipped1 = ".".join(file.split(".")[:-1])
				unzipped2 = ".".join(pair.split(".")[:-1])
				command = ('gunzip -c ' + file  +  ' | python $PYTHON_PATH/remove_seq_adapter.py GATCGGAAGAGCGG' + ' > ' + unzipped1 + '\n'
							'gunzip -c ' + pair + ' | python $PYTHON_PATH/remove_seq_adapter.py GATCGGAAGAGCGG' + ' > ' + unzipped2 + '\n')
				QSUB.write(command)

				#command = '$BOWTIE_PATH/bowtie -q -S -t -n 2 -e 200 --chunkmbs 256 -3 0 -5 16 -m 1 --best --strata hybrid_index ' + unzipped1 + ' $RESULTS_PATH/' + samfile + '\n'
				command = ('$BOWTIE_PATH/bowtie -q -S -t -n 2 -e 200 --chunkmbs 256 -3 0 -5 16' 
								' -m 1 -X 800 --allow-contain --best --strata hybrid_index'
				 				' -1 ' + unzipped1 + ' -2 ' + unzipped2 +  ' $RESULTS_PATH/' + samfile + '\n')
				QSUB.write(command)

				command = 'rm ' + unzipped1 + '\nrm ' + unzipped2 + '\n'
				QSUB.write(command)

				QSUB.close()

				time.sleep(1)
				thisCommand = "chmod 755 " + qsubFile
				os.system(thisCommand)
				thisCommand = "qsub " + qsubFile
				os.system(thisCommand)






if __name__ == "__main__":
	main()


