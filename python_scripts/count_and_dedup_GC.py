#!/usr/bin/env python

import sys
import bisect


def main():

	print("in count_and_dedupGC")
	infilename = sys.argv[1]
	GCfilename =sys.argv[2]
	outfilename = sys.argv[3]
	statfilename = sys.argv[4]

	chrominfo = fileToDictionary("/mnt/wigclust1/data/safe/kostic/bin_mapping/chrom_sizes_hg_dm_combined_consecutive.txt", 0)
	bins = fileToArray("/mnt/wigclust1/data/safe/kostic/bin_mapping/hybrid_bin_boundaries_sorted_125_600.txt", 0)
	GC_bins = fileToArray("/mnt/wigclust1/data/safe/kostic/bin_mapping/GC_bin_bounds_125_600.txt", 0)

	INFILE = open(infilename, "r")
	GCFILE = open(GCfilename, "r")
	OUTFILE = open(outfilename, "w")
	STATFILE = open(statfilename, "w")

	#make 100 GC/len bins for each genomic bin
	binCounts = []
	gcbinCounts = []
	print("GC_bins length: " + str(len(GC_bins)))

	for j in range(len(GC_bins)):
		gcbinCounts.append(0)
	for i in range(len(bins)):
		binCounts.append(gcbinCounts)

	#the absolute start positions of bins on the chromosomes
	binStarts = []
	for i in range(len(bins)):
		binStarts.append(long(bins[i][2]))

	# gcStarts = []
	lenStarts = []
	for i in range(len(GC_bins)):
		lenStarts.append( float(GC_bins[i][0]) + float(GC_bins[i][2]) )

	counter = 0
	dups = 0
	totalReads = 0
	prevChromposTags = dict()
	prevChromPosits = dict()

	#skip first line of file - labels
	GCFILE.readline()
	
	for x in INFILE:

		pair = INFILE.next()
		arow = x.rstrip().split("\t")
		pairrow = pair.rstrip().split("\t")

		if (int(arow[8]) <= 0)
		####somethign here TODO

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
		if arow[8] = 0:
			continue
		if chrominfo.has_key(thisChrom):
			pass
		else:
			continue

		#the flag
		thisCode = int(arow[1])
		code2 = int(pairrow[1])

		thisTag = ""
		if len(arow) > 11:
			for i in range(11, len(arow)):
				if arow[i][0:5] == "XV:Z:":
					thisTag = arow[i][5:]
					break
		thisGC = 0.0
		if len(arow) > 12:
			if arow[-1][0:3] == "GC:":
				thisGC = float(arow[-1][4])

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
		indexGenome = bisect.bisect(binStarts, thisAbspos)

		#for each alignment, get lenght of frag from SAM format and calc GC content
		fl = abs(int(arow[8]))
		lenIndex = bisect.bisect(lenStarts, float(thisGC) + float(fl))
		binCounts[indexGenome-1][lenIndex-1] += 1


	###### end for line in the sam file

	for i in range(len(binCounts)):
		for j in range(binCounts[i]):
			#chrom name, start position of genomic bin, absolute start positions
			OUTFILE.write("\t".join(bins[i][0:3]))
			OUTFILE.write("\t")
			#gc content min, gc max, len min, len max
			OUTFILE.write("\t".join(GC_bins[j][0:4]))
			OUTFILE.write("\t")
			OUTFILE.write(str(binCounts[i][j]))
			OUTFILE.write("\t")




		thisRatio = float(binCounts[i]) / ( float(counter) / float(len(bins)) )
		#print chromosome name, start pos, and absolute start pos
		OUTFILE.write("\t".join(bins[i][0:3]))
		OUTFILE.write("\t")
		OUTFILE.write(str(binCounts[i][j]))
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
