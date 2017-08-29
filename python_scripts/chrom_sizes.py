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
	outfilename = sys.argv[2]

	hg19_chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
	dm_chroms = ["2L", "2R", "3L", "3R" , "4", "X", "Y"]

	INFILE = 0
	OUTFILE = open(outfilename, "w")

	abspos = 0
	for x in hg19_chroms:
		infilename = directory + "/chr" + x + ".fa"
		print infilename
		if INFILE:
			INFILE.close()
		INFILE = open(infilename, "r")
		INFILE.readline()
		chromlength = 0
		for y in INFILE:
			aline = y.rstrip()
			chromlength += len(aline)
		OUTFILE.write("chr" + x)
		OUTFILE.write("\t")
		OUTFILE.write(str(chromlength))
		OUTFILE.write("\t")
		OUTFILE.write(str(abspos))
		abspos += chromlength
		OUTFILE.write("\n")
		
		
	INFILE.close()

	INFILE1 = 0
	for x in dm_chroms:
		infilename = directory + "/chr" + x + "_dm.fa"
		print infilename
		if INFILE1:
			INFILE1.close()
		INFILE1 = open(infilename, "r")
		INFILE1.readline()
		chromlength = 0
		for y in INFILE1:
			aline = y.rstrip()
			chromlength += len(aline)
		OUTFILE.write("chr" + x + "_dm")
		OUTFILE.write("\t")
		OUTFILE.write(str(chromlength))
		OUTFILE.write("\t")
		OUTFILE.write(str(abspos))
		OUTFILE.write("\n")
		abspos += chromlength

	OUTFILE.close()



if __name__ == "__main__":
	main()
