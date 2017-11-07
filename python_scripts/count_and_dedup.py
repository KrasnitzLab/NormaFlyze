#!/usr/bin/env python

import sys
import bisect


def main():

	infilename = sys.argv[1]
	outfilename = sys.argv[2]
	statfilename = sys.argv[3]

	chrominfo = fileToDictionary("/mnt/wigclust1/data/safe/kostic/bin_mapping/chrom_sizes_hg_dm_combined_consecutive.txt", 0)
	bins = fileToArray("/mnt/wigclust1/data/safe/kostic/bin_mapping/hg19_bin_boundaries_sorted_125_600.txt", 0)
		#hybrid_bin_boundaries_sorted_125_600.txt", 0)

	INFILE = open(infilename, "r")
	OUTFILE = open(outfilename, "w")
	STATFILE = open(statfilename, "w")

	binCounts = []
	for i in range(len(bins)):
		binCounts.append(0)

	#the absolute start positions of bins on the chromosomes
	binStarts = []
	for i in range(len(bins)):
		binStarts.append(long(bins[i][2]))

	counter = 0
	dups = 0
	totalReads = 0
	prevChromposTags = dict()
	prevChromPosits = dict()
	
	for x in INFILE:

		pair = INFILE.next()
		arow = x.rstrip().split("\t")
		pairrow = pair.rstrip().split("\t")

		if not pairrow[0] == arow[0]:
			print "pairs not sorted right"
			break

		thisChrom = arow[2]
		thisChrompos = min(arow[3],pairrow[3])
		
		#the flag
		thisCode = int(arow[1])
		code2 = int(pairrow[1])
		
		thisTag = ""
		if len(arow) > 11:
			for i in range(11, len(arow)):
				if arow[i][0:5] == "XV:Z:":
					thisTag = arow[i][5:]
					break
		
		if thisChrom.find("random") > -1 or thisChrom.find("Un") > -1:
			continue
		if thisChrom == "chrM":
			continue
		if thisChrom == "":
			continue
		if chrominfo.has_key(thisChrom):
			pass
		else:
			continue

		totalReads += 1

		#if this read maps to same place as a previous one, check if tags are the same
		#if tags are same, this is a PCR duplicate, do not count in the bins
		try:
			if thisChrompos in prevChromPosits[thisChrom] and prevChromposTags.has_key(thisTag):
				prevChromposTags[thisTag] += 1
				dups += 1
				continue
		except KeyError:
			prevChromPosits[thisChrom] = []

		counter += 1
			
		prevChromPosits[thisChrom].append(thisChrompos)
		prevChromposTags[thisTag] = 1
		
		thisChrominfo = chrominfo[thisChrom]
		thisAbspos = long(thisChrompos) + long(thisChrominfo[2])
		
		#get index where thisAbspos would fit in the binStarts array
		#increment the bin count of the bin where this read would fall
		indexDown = bisect.bisect(binStarts, thisAbspos)
		binCounts[indexDown-1] += 1


		#if the paired other side falls out of the bins, keep track - make a counter and put in stats

	###### end for line in the sam file

	for i in range(len(binCounts)):
		thisRatio = float(binCounts[i]) / ( float(counter) / float(len(bins)) )
		#print chromosome name, start pos, and absolute start pos
		OUTFILE.write("\t".join(bins[i][0:3]))
		OUTFILE.write("\t")
		OUTFILE.write(str(binCounts[i]))
		OUTFILE.write("\t")
		OUTFILE.write(str(thisRatio))
		OUTFILE.write("\n")

	binCounts.sort()
	
	STATFILE.write("TotalReads\tDupsRemoved\tReadsKept\tMedianBinCount\n")
	STATFILE.write(str(totalReads))
	STATFILE.write("\t")
	STATFILE.write(str(dups))
	STATFILE.write("\t")
	STATFILE.write(str(counter))
	STATFILE.write("\t")
	STATFILE.write(str(binCounts[len(bins)/2]))
	STATFILE.write("\n")

	INFILE.close()
	OUTFILE.close()
	STATFILE.close()


def fileToDictionary(inputFile, indexColumn):
	input = open(inputFile, "r")

	rd = dict()
	for x in input:
		arow = x.rstrip().split("\t")
		id = arow[indexColumn]
		if rd.has_key(id):
			print "duplicate knowngene id = " + id
			print "arow =   " + str(arow)
			print "rd[id] = " + str(rd[id])
		else:
			rd[id] = arow
		
	input.close()
	return(rd)


def fileToArray(inputFile, skipFirst):
	input = open(inputFile, "r")

	ra = []

	for i in range(skipFirst):
		input.readline()

	for x in input:
		arow = x.rstrip().split("\t")
		ra.append(arow)
		
	input.close()
	return(ra)


if __name__ == "__main__":
	main()
