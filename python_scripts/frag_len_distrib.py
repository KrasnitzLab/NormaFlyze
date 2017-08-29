#!/usr/bin/env python

# this script finds the length of the chromosomes in the input files
# output: a file where each line contains "chrom name:lenght of chrom:asbolute position"

import re
import sys
import random
import string
import os

def main():

	directory = sys.argv[1]
	#outfilename = "data_frag_len_distrib_hg_dm.txt"

	len_density = {}
	dm_density = {}

#
	hg_lens = []
	dm_lens = []
#

	for f in os.listdir(directory):
		if os.path.isfile(os.path.join(directory,f)) and ("aligned" in f):
			print f
			#
			outfilename = "len_distrib_" + f
			#
			IN = open(directory + "/" + f, "r")
			for line in IN:
				if "@" in line:
					continue

				info = line.split("\t")
				chrom = info[2].strip()
				length = info[8].strip()
				
				if (not "-" in length) and (not int(length) == 0):
					if "_dm" in chrom:
						dm_lens.append(length)
						# try:
						# 	dm_density[length] += 1
						# except KeyError:
						# 	dm_density[length] = 1
					else:
						hg_lens.append(length)
						# try:
						# 	len_density[length] += 1
						# except KeyError:
						# 	len_density[length] = 1
				else:
					continue

#####
			OUT = open(outfilename, "w")
			for i in range(len(hg_lens)):
				OUT.write("hg\t" + hg_lens[i] + "\n")
			for i in range(len(dm_lens)):
				OUT.write("dm\t" + dm_lens[i] + "\n")
			OUT.close
#####		

	# OUT = open(outfilename, "w")
	# for key,val in len_density.items():
	# 	OUT.write("hg\t" + key + "\t" + str(val) + "\n")
	# for key,val in dm_density.items():
	# 	OUT.write("dm\t" + key + "\t" + str(val) + "\n")
	# OUT.close()



if __name__ == "__main__":
	main()
