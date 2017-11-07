#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import bisect
import numpy

def main():

	GCbincount = 10
	LENbincount = 10

	#the chrom lenghts and relative positions
	CHROMLEN = open(sys.argv[1], "r")

	#GC content of each fragment
	GC = open(sys.argv[2], "r")

	OUTFILE = open(sys.argv[3], "w")

	f_list = makeList(GC)
	num_f = len(f_list)

	gcB = int(num_f / GCbincount)
	rem = num_f % GCbincount

	#sort the list by GC and separate into 10 bins
	gc_list = sorted(f_list, key = itemgetter(1))
	print(gc_list)
	gc_bins = numpy.array_split(gc_list, GCbincount)
	#print("gc_bin number: " + str(len(gc_bins)))

	#sort GC bins by fragment length and separate into 10 sub-bins
	lenB = int(gcB / LENbincount)
	rem = gcB % LENbincount

	len_bins = []
	for i in range(0, len(gc_bins)):
		this_sub = gc_bins[i].tolist()
		gc_min = this_sub[0][1]
		gc_max = this_sub[-1][1]
		this_sub.sort(key = itemgetter(2))
		split_this = numpy.array_split(this_sub, LENbincount)
		len_bins.append(split_this)
		print("len list length: " + str(len(len_bins)))

		for j in len_bins[i]:
			lst = j.tolist()
			firstEntry = lst[0]
			lastEntry = lst[-1]
			OUTFILE.write(gc_min + "\t")
			OUTFILE.write(gc_max + "\t")
			OUTFILE.write(firstEntry[2] + "\t")
			OUTFILE.write(lastEntry[2] + "\t")
			OUTFILE.write(str(len(lst)) + "\n")

	OUTFILE.close()
	CHROMLEN.close()
	GC.close()

def makeList(GC_fragments):
	d = []
	GC_fragments.readline()
	for line in GC_fragments:
		sep = line.split("\t")
		name = "\t".join(sep[0:3])
		gc_content = float(sep[3])
		length = int(sep[2]) - int(sep[1])
		vals = (name, gc_content, length)
		d.append(vals)
	return(d)


if __name__ == "__main__":
	main()
