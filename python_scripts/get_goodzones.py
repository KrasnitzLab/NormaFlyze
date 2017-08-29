#!/usr/bin/env python

import re
import sys
import random
import string
import time

def main():

	file = sys.argv[1]
	outfile = "hybrid_unique_mappers_" + time.strftime("%m_%d") + ".bed"
	mapfile = "hybrid_num_mappable_" + time.strftime("%m_%d") + ".txt"
	OUTFILE = open(outfile, "w")
	MAPFILE = open(mapfile, "w")
	
	inGood = False
	prevChrom = ""
	thisStart = -1
	thisEnd = -1
	count = 0

	#the .sam alignment file
	INFILE = open(file, "r")
	for line in INFILE:

		#skip past the header section of the SAM file
		if line[0] == "@":
			continue

		pair_line = INFILE.next()

		if "random" in line or "Un" in line or "chrM" in line:
			continue
		
		read1_data = line.rstrip().split()
		read2_data = pair_line.rstrip().split()
		######
		read1_id = read1_data[0]
		read2_id = read2_data[0]

		read1_chrom = read1_id.split(".")[0]
		read2_chrom = read2_id.split(".")[0]

		seg_start1 = int(read1_id.split(".")[1])
		seg_start2 = int(read2_id.split(".")[1])

		seg_end1 = int(read1_id.split(".")[2])
		seg_end2 = int(read2_id.split(".")[2])
		######
		code1 = int(read1_data[1])
		code2 = int(read2_data[1])
		######
		align_chrom1 = read1_data[2]
		align_chrom2 = read2_data[2]
		######
		align_posit1 = int(read1_data[3])
		align_posit2 = int(read2_data[3])
		######
		rev_read_len = len(read2_data[9])


		#simplify this by just spliitting at / and checking if same TODO
		if (not read1_id == read2_id) or not (align_chrom1 == align_chrom2):
			print("pair lines not in correct order")
			print("prob with: " + read1_id + " and " + read2_id)
			break

		if prevChrom == "":
			prevChrom = align_chrom1


		# code 99 means forward alignment of pair 1 successful
		# code 147 means reverse alignment of pair 2 successful
		if (code1 == 99 and code2 == 147
			and (align_posit1 - 1) == seg_start1
			and (align_posit2 - 1) == seg_end2-rev_read_len):

			OUTFILE.write(align_chrom1)
			OUTFILE.write("\t")
			OUTFILE.write(str(seg_start1))
			OUTFILE.write("\t")
			OUTFILE.write(str(seg_end2))
			OUTFILE.write("\n")

			if not (prevChrom == align_chrom1) and not(prevChrom == "*"):
				MAPFILE.write(prevChrom + "\t" + str(count) + "\n")	
				count = 0

			prevChrom = align_chrom1
			count = count+1

	MAPFILE.write(align_chrom1 + "\t" + str(count) + "\n")

	INFILE.close()
	MAPFILE.close()
	OUTFILE.close()


if __name__ == "__main__":
	main()
