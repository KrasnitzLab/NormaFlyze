#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import bisect


def main():

	bincount = 5000
	#numner of mappable positions on each chromosome
	MAP = open(sys.argv[1], "r")
	#all the mappable positions
	GOOD = open(sys.argv[2], "r")
	#the chrom lenghts and relative positions
	CHROMLEN = open(sys.argv[3], "r")
	directory = sys.argv[4]

	outfile = "/hybrid_bin_boundaries_" + time.strftime("%m_%d") + ".txt"
	OUTFILE = open(directory + outfile, "w")

	chromlen = dict()
	for line in CHROMLEN:
		arow = line.rstrip().split("\t")
		thisChrom = arow[0]
		thisChromlen = int(arow[1])
		thisAbspos = long(arow[2])
		chromlen[thisChrom] = [thisChromlen, thisAbspos]

	print("chromlen dict")
	print chromlen

	chromarray = []
	chroms = dict()
	totalUnique = long(0)
	#put total num mappable posits per chrom into a dictionary
	for line in MAP:
		arow = line.rstrip().split("\t")
		thisChrom = arow[0]
		numUniquePos = long(arow[1])
		if "chrM" in thisChrom or "random" in thisChrom or "Un" in thisChrom:
			continue
		totalUnique += numUniquePos
		chroms[thisChrom] = numUniquePos

	bincountUsed = 0
	for k, v in chroms.items():
		chromBincount = float(bincount) * (float(v) / float(totalUnique))
		i = int(chromBincount)
		# print i
		bincountUsed += i 
		r = chromBincount - i
		chromarray.append([k, i, r])

	# Sort chromarray from highest to lowest r value
	a = []
	for i in chromarray:
		bisect.insort(a, (-i[2], i))

	chromarray = []
	for j in a:
		chromarray.append(j[1]) 
		# This j[1] index has to match the index of i in the bisect.insort 2nd parameter.

	a = []
	
	# print chromarray
	# done sorting chromarray

	# Add the remaining reads to bincount of chromosomes in that range
	remain = bincount - bincountUsed
	for i in range(remain):
		chromarray[i][1] += 1

	chroms2 = dict()

	for i in range(len(chromarray)):
		numUniquePos = chroms[chromarray[i][0]]
		chrombins = chromarray[i][1]
		binlength = numUniquePos / chrombins
		remainder = numUniquePos % chrombins
		chroms2[chromarray[i][0]] = [chrombins, binlength, remainder]

	
	chromorder = sorted(chroms2.keys())
	print("chromorder")
	print chromorder

	print 
	print "Starting to get bin boundaries"
	print

	for chrom in chromorder:

		currChrom = chrom
		binStart = 0
		binEnd = 0
		numBins = 0
		numreads = 0

		# print("current chrom: " + currChrom)
		# print("num bins needed: " + str(chroms2[currChrom][0]))
		# print("num reads per bin: " + str(chroms2[currChrom][1]))

		while numBins < chroms2[currChrom][0]:
			binReads = chroms2[currChrom][1]
			numReads = 0

			if chroms2[currChrom][2] > 0:
					binReads += 1
					chroms2[currChrom][2] -= 1

			while numReads < binReads:

				currEntry = GOOD.readline()
				if currEntry == '':
					break
				currChrom = currEntry.split()[0]
				binEnd = int(currEntry.split()[2])
				if not currChrom == chrom:
					print("something wrong with counting")
				numReads += 1

			#print("done reading for this bin, num reads: " + str(numReads))

			binStartAbspos = chromlen[currChrom][1] + binStart

			#for the last bin: ends at end of chromosome
			if numBins+1 == chroms2[currChrom][0]:
				binEnd = chromlen[currChrom][0]
				# print("binEnd for last: "+ str(binEnd))
				# print("last segment " + currEntry)

			OUTFILE.write(currChrom)
			OUTFILE.write("\t")
			OUTFILE.write(str(binStart))
			OUTFILE.write("\t")
			OUTFILE.write(str(binStartAbspos))
			OUTFILE.write("\t")
			OUTFILE.write(str(binEnd))
			OUTFILE.write("\t")
			#actual length of the bin
			OUTFILE.write(str(binEnd - binStart))
			OUTFILE.write("\t")
			#number of reads in the bin
			OUTFILE.write(str(numReads))
			OUTFILE.write("\n")

			binStart = binEnd
			numBins += 1


	OUTFILE.close()
	MAP.close()
	GOOD.close()
	CHROMLEN.close()


if __name__ == "__main__":
	main()
