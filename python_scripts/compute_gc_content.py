#!/usr/bin/env python

import re
import sys
import random
import string
import time


def main():

	directory = sys.argv[1]
	bounds = sys.argv[2]

	bins = fileToArray(bounds, False)
	INFILE = False
	OUTFILE = open("hybrid_genome_GCcontent" + time.strftime("%m_%d") + ".txt","w")
	OUTFILE.write("bin.chrom\tbin.start.chrompos\tbin.start.abspos\tbin.end.chrompos\tbin.length\tmappable.positions\tgc.content\n")

	prevChrom = ""

	for line in bins:
		currChrom = line[0].strip()
		binStart = int(line[1].strip())
		binEnd = int(line[3].strip())

		if not currChrom == prevChrom:
			if INFILE:
				INFILE.close()
			INFILE = open(directory+ "/" + currChrom + ".fa", "r")
			chr = []
			x = ""
			y = ""
			INFILE.readline()
			for x in INFILE:
				chr.append(x.rstrip())
			x = "".join(chr)
			y = x.upper()
			print "after read " + currChrom
			prevChrom = currChrom

		gcContent = float(len(re.findall("[CG]", y[(binStart):(binEnd+1)]))) / float(len(re.findall("[ACGT]", y[(binStart):(binEnd+1)])))
		OUTFILE.write("\t".join(line))
		OUTFILE.write("\t")
		OUTFILE.write(str(gcContent))
		OUTFILE.write("\n")
		OUTFILE.flush()

	INFILE.close()
	OUTFILE.close()


def fileToArray(inputFile, skipFirst):
	input = open(inputFile, "r")

	ra = []
	if skipFirst:
		input.readline()
	for x in input:
		arow = x.rstrip().split("\t")
		ra.append(arow)
	input.close()
	return(ra)


if __name__ == "__main__":
	main()
