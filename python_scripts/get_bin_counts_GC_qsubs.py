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
			qsubFile = bash_dir + "/count_GC_" + fname + ".sh"
			print(qsubFile)
			QSUB = open(qsubFile, "w")

			header = ('#!/bin/bash \n#$ -q "all.q@wigclust1[7-9]","all.q@wigclust2[0-4]" '
							'\n#$ -pe threads 1 '
							'\n#$ -l vf=8G '
							'\n#$ -o /mnt/wigclust1/data/safe/kostic/output/count_''' + fname + 'GC.out -j y \n')
			QSUB.write(header)

			paths = ('export DATA_PATH=/mnt/wigclust1/data/safe/kostic/SNS_data_2' 
							'\nexport PYTHON_PATH=/mnt/wigclust1/data/safe/kostic/python_scripts'
							'\nexport SAM_PATH=/mnt/wigclust1/data/software/samtools/samtools-0.1.19'
							'\nexport GENOME_PATH=/mnt/wigclust1/data/safe/kostic/hybrid_genome'
							'\ncd $DATA_PATH \n')
			QSUB.write(paths)

			command = ('$SAM_PATH/samtools sort ' + f + ' ' + fname + '_sorted\n'
					 + '$PYTHON_PATH/compute_frag_gc_content.py <( $SAM_PATH/samtools view ' + fname +'_sorted.bam ) $GENOME_PATH ' + fname +'_GCadded.sam\n')

			#QSUB.write(command)

			#calc GC content of each fragment
			keep = '$3,"\\t",$4,"\\t",$4+$9'
			#command = ('$SAM_PATH/samtools view -f 0x23 ' + f + " | awk '{print " + keep + " }' > bedFile_" + fname + ".bed\n")
							#+ 'python $PYTHON_PATH/compute_frag_gc_content.py bedFile_' + fname +'.bed $GENOME_PATH ' + fname + '_GCcontent.txt\n')
			#QSUB.write(command)

			#in count_and_dedup.py, change the input bin_boundaries file and the chromosome length file to the appropriate ones
			command = ('python $PYTHON_PATH/count_and_dedup_GC.py <( $SAM_PATH/samtools view ' + f + ' ) ' + fname + '_GCcontent.txt'
						' ' + fname + '_varbinGC_count.txt ' + fname + '_varbinGC_stats.txt\n'
						+ 'rm bedFile_' + fname + '.bed\n' )
			#QSUB.write(command)

			QSUB.close()

			time.sleep(1)
			thisCommand = "chmod 755 " + qsubFile
			os.system(thisCommand)
			thisCommand = "dos2unix " + qsubFile
			os.system(thisCommand)
			thisCommand = "qsub " + qsubFile
			os.system(thisCommand)


if __name__ == "__main__":
	main()