#!/usr/bin/env python

import sys
import bisect
import numpy


def main():

	print("in count_and_dedupGC")
	infilename = sys.argv[1]
	outfilename = sys.argv[2]
	statfilename = sys.argv[3]

	chrominfo = fileToDictionary("/mnt/wigclust1/data/safe/kostic/bin_mapping/chrom_sizes_hg_dm_combined_consecutive.txt", 0)
	bins = fileToArray("/mnt/wigclust1/data/safe/kostic/bin_mapping/hybrid_bin_boundaries_sorted_125_600.txt", 0)
	GC_bins = fileToArray("/mnt/wigclust1/data/safe/kostic/bin_mapping/GC_bin_bounds_125_600.txt", 0)

	INFILE = open(infilename, "r")
	OUTFILE = open(outfilename, "w")
	STATFILE = open(statfilename, "w")

	print infilename

	#make 100 GC/len bins for each genomic bin
	binCounts = []
	gcbinCounts = []


	binCounts = numpy.zeros((len(bins), len(GC_bins)))

	#the absolute start positions of bins on the chromosomes
	binStarts = []
	for i in range(len(bins)):
		binStarts.append(long(bins[i][2]))

	lenStarts = []
	for i in range(len(GC_bins)):
		lenStarts.append( float(GC_bins[i][0]) + float(GC_bins[i][2]) )

	counter = 0
	dups = 0
	totalReads = 0
	flycount = 0
	flyNLA = 0
	humancount = 0
	humanNLA = 0
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
		#the leftmost alignment position of the fragment
		thisChrompos = min(arow[3],pairrow[3])
		
		if thisChrom.find("random") > -1 or thisChrom.find("Un") > -1:
			continue
		if thisChrom == "chrM":
			continue
		if thisChrom == "":
			continue
		if arow[8] == 0:
			continue
		if chrominfo.has_key(thisChrom):
			pass
		else:
			continue

		#the flag
		thisCode = int(arow[1])
		code2 = int(pairrow[1])

		thisTag = ""
		thisGC = 0.0
		if len(arow) > 11:
			for i in range(11, len(arow)):
				if arow[i][0:5] == "XV:Z:":
					thisTag = arow[15][5:]
				if arow[i][0:3] == "GC:":
					thisGC = float(arow[i][3:])
				if arow[i][0:4] == "NLA:":
					thisNLA = int(arow[i][4:])

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

		if (thisChrom.find("_dm") > -1):
			flycount += 1
			flyNLA += thisNLA
		else:
			humancount += 1
			humanNLA += thisNLA
			
		prevChromPosits[thisChrom].append(thisChrompos)
		prevChromposTags[thisTag] = 1
		
		thisChrominfo = chrominfo[thisChrom]
		thisAbspos = long(thisChrompos) + long(thisChrominfo[2])
		
		#get index where thisAbspos would fit in the binStarts array
		indexGenome = bisect.bisect(binStarts, thisAbspos)
		#for each alignment, get lenght of frag from SAM format
		fl = abs(int(arow[8]))
		lenIndex = bisect.bisect(lenStarts, float(thisGC) + float(fl))
		binCounts[indexGenome-1][lenIndex-1] += 1


	###### end for line in the sam file
	print("starting to add to output file")
	for p in range(len(binCounts)):
		for q in range(len(binCounts[p])):
			#chrom name, start position of genomic bin, absolute start position of genome, start GC, end GC, start len, end len, bin count for GC/len, total count for genomic bin
			OUTFILE.write("\t".join(bins[p][0:3]))
			OUTFILE.write("\t")
			#gc content min, gc max, len min, len max
			OUTFILE.write("\t".join(GC_bins[q][0:4]))
			OUTFILE.write("\t")
			OUTFILE.write(str(binCounts[p][q]))
			OUTFILE.write("\t")
			OUTFILE.write(str(sum(binCounts[p])))
			OUTFILE.write("\n")

	
	STATFILE.write("TotalReads\tDupsRemoved\tReadsKept\tFlyNLAEff\tMeanHumanNLAEff\n")
	STATFILE.write(str(totalReads))
	STATFILE.write("\t")
	STATFILE.write(str(dups))
	STATFILE.write("\t")
	STATFILE.write(str(counter))
	STATFILE.write("\t")
	#print("final fly: " + str(float(flyNLA)/float(flycount)))
	STATFILE.write(str(float(flycount-1)/float(flyNLA + flycount - 1)))
	STATFILE.write("\t")
	#print("final hum: " + str(float(humanNLA/humancount)))
	STATFILE.write(str(float(humancount-1)/float(humanNLA + humancount - 1)))
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
