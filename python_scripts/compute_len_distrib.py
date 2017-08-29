#!/usr/bin/env python

import re
import sys
import random
import string
import numpy
import time



def main():

	bounds = sys.argv[1]
	frags = sys.argv[2]

	bins = fileToArray(bounds, False)
	f_dict = fileToDict(frags, False)
	OUTFILE = open("hybrid_bin_fragment_stats" + time.strftime("%m_%d") + ".txt","w")
	OUTFILE.write("bin.chrom\tbin.start.chrompos\tbin.start.abspos\tbin.end.chrompos\tmedian.len\tmedian.len.ratio\n")

	med = {}
	mn = {}
	med_sum = 0
	mean_sum = 0

	prevpos = 0
	prevChrom = ""
	for line in bins:

		currChrom = line[0].strip()
		binStart = int(line[1].strip())
		binEnd = int(line[3].strip())
		lens = []
		if not currChrom == prevChrom:
			prevpos = 0
			prevChrom = currChrom

		fragments = f_dict[currChrom]

		while int(fragments[prevpos][0]) >= binStart and int(fragments[prevpos][1]) <= binEnd:
			lens.append( int(fragments[prevpos][1]) - int(fragments[prevpos][0]) )
			prevpos += 1
			print "appending"
			if prevpos >= len(fragments):
				print("broke at " + str(prevpos) + " listlen " + str(len(fragments)))
				break

		key = "\t".join(line[0:2])
		med[key] = numpy.median(lens)
		mn[key] = numpy.mean(lens)
		med_sum += med[key] 
		mean_sum += mn[key]
		print("median " + str(med[key]))

	for line in bins:

		key = "\t".join(line[0:2])

		OUTFILE.write("\t".join(line[:-2]))
		OUTFILE.write("\t")
		OUTFILE.write(str( med[key] ))
		OUTFILE.write("\t")
		OUTFILE.write(str( float( med[key] / med_sum ) ))
		OUTFILE.write("\n")
		OUTFILE.flush()

	OUTFILE.close()


def fileToDict(inputFile, skipFirst):
	infile = open(inputFile, "r")
	outDict = {}

	if skipFirst:
		infile.readline()
	for line in infile:
		arow = line.rstrip().split("\t")
		chrom = arow[0]
		try:
			outDict[chrom].append((arow[1],arow[2]))
		except KeyError:
			outDict[chrom] = []
			outDict[chrom].append((arow[1],arow[2]))

	infile.close()
	return(outDict)



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
