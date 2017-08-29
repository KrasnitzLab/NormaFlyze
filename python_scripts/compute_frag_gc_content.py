#!/usr/bin/env python

import re
import sys
import random
import string
import time


def main():

	print("computing GC content")
	start_time = time.time()
	print(str(start_time))

	samfile = sys.argv[1]
	SAM = open(samfile, "r")

	directory = sys.argv[2]

	outfile = sys.argv[3]
	GC_content = open(outfile, "w")

	chrom_file = False
	prevChrom = ""
	inHeader = True

	for line in SAM:

		if inHeader:
			if line[0] != "@":
				inHeader = False
			else:
				GC_content.write(line.rstrip())
				continue

		arow = line.rstrip().split("\t")
		pair_line = SAM.next()
		pairrow = pair_line.rstrip().split("\t")

		if not arow[0] == pairrow[0]:
			print("not ordered correctly")

		currChrom = arow[2]
		fragStart = min(int(arow[3]), int(pairrow[3]))
		fragEnd = fragStart + abs(int(arow[8]))

		if not currChrom == prevChrom:
			if chrom_file:
				chrom_file.close()
			chrom_file = open(directory + "/" + currChrom + ".fa", "r")
			chr = []
			x = ""
			y = ""
			chrom_file.readline()
			for x in chrom_file:
				chr.append(x.rstrip())
			x = "".join(chr)
			y = x.upper()
			prevChrom = currChrom

		print("curr chrom: " + currChrom)
		print(str(fragStart) +"\t"+ str(fragEnd+1))
		gcContent = float(len(re.findall("[CG]", y[(fragStart):(fragEnd+1)]))) / float(len(re.findall("[ACGT]", y[(fragStart):(fragEnd+1)])))

		arow.append("GC:" + str(gcContent))
		pairrow.append("GC:" + str(gcContent))

		GC_content.write("\t".join(arow) + "\n")
		GC_content.write("\t".join(pairrow) + "\n")
		GC_content.flush()

	SAM.close()
	GC_content.close()

	print("done computing GC content, time elapsed = " + str(time.time() - start_time))


if __name__ == "__main__":
	main()
