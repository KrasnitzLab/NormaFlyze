#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import os
import subprocess

def main():

	infilename = "/mnt/wigclust1/data/safe/kostic/SNS_data_2/nla3_hybrid_guide_02.txt"
	#hg19_only/nla3_hg19_guide_02.txt"
	#nla3_hybrid_guide_02.txt"
	output_dir = "/mnt/wigclust1/data/safe/kostic/SNS_data_2/"
	#hg19_only"
	
	outfilename = "/mnt/wigclust1/data/safe/kostic/SNS_data_2/GClen_uber_10by10_count_data.txt"
	#hg19_only/GClen_uber_hg19_count_data.txt"
	#GClen_uber_10by10_count_data.txt"

	nlaout = "/mnt/wigclust1/data/safe/kostic/SNS_data_2/GClen_nla_efficiencies.txt"
	NLA = open(nlaout, "w")

	guide = fileToGuide(infilename)
	ncells = len(guide["seq.unit.id"])

	uber = dict()
	#index of col where the bincount is found
	colnum = 7
	samples = []

	NLA3fly = []
	NLA3hum = []

	for i in range(ncells):
		if guide["process"][i] == "0":
			continue

		samples.append(guide["seq.unit.id"][i])
		uber[guide["seq.unit.id"][i]] = []
		thisFilename = output_dir + "/" + guide["barcode.group"][i] + "_" + guide["barcode"][i] + "_varbinGC_10by10_count.txt"

		print thisFilename

		statFilename = output_dir + "/" + guide["barcode.group"][i] + "_" + guide["barcode"][i] + "_varbinGC_10by10_stats.txt"

		# "_hg19_count.txt"
		# "_hg19_stats.txt"

		# thisFilename = output_dir + "/" + guide["barcode.group"][i] + "_" + guide["barcode"][i] + "_hg19_varbinGC_count.txt"

		# statFilename = output_dir + "/" + guide["barcode.group"][i] + "_" + guide["barcode"][i] + "_hg19_varbinGC_stats.txt"


		STAT = open(statFilename,"r")
		STAT.readline()
		#NLAinfo = STAT.readline().split("\t")
		#NLA3hum.append(float(NLAinfo[4]))
		#NLA3fly.append(float(NLAinfo[3]))


		IN = open(thisFilename, "r")
		rowcount = 0
		for x in IN:
			arow = x.rstrip().split("\t")
			#append the bin count
			uber[guide["seq.unit.id"][i]].append(arow[colnum])
			rowcount += 1
		IN.close()
		if rowcount < 1:
			print "ERROR: ", samples[len(samples) - 1]
		
		OUT = open(outfilename, "w")
	
	print samples

	for j in range(len(samples) - 1):
		OUT.write(samples[j])
		OUT.write("\t")
		NLA.write(samples[j])
		NLA.write("\t")

	OUT.write(samples[-1])
	OUT.write("\n")
	#NLA.write(samples[-1])
	#NLA.write("\n")

	# NLA.write("human Nla3 eff\t")
	# for l in range(len(NLA3hum)):
	# 	NLA.write(str(NLA3hum[l]) + "\t")
	# NLA.write("\nhuman Nla3 eff\t")
	# for l in range(len(NLA3hum)):
	# 	NLA.write(str(NLA3fly[l]) + "\t")
	# NLA.write("\n")

	for k in range(len(uber[samples[0]])):
		for j in range(len(samples) - 1):
			OUT.write(uber[samples[j]][k])
			OUT.write("\t")
		OUT.write(uber[samples[-1]][k])
		OUT.write("\n")

	OUT.close()
	STAT.close()
	NLA.close()
	

def fileToGuide(infilename):
	INFILE = open(infilename, "r")
	
	guide = dict()
	x = INFILE.readline()
	colnames = x.rstrip().split("\t")
	for c in colnames:
		c = c.strip()
	ncols = len(colnames)

	a = []
	for x in INFILE:
		arow = x.rstrip().split("\t")
		for b in arow:
			b = b.strip()
		if len(arow) < ncols:
			n = len(arow)
			for i in range(n, ncols):
				arow.append("")
		a.append(arow)

	z = zip(*a)
	if len(z) == ncols:
		for i in range(ncols):
			if guide.has_key(colnames[i]):
				print "ERROR: Duplicate column name."
			else:
				guide[colnames[i]] = z[i]
	else:
		print "ERROR: Data and colnames length not equal."
		
	INFILE.close()
	return guide


if __name__ == "__main__":
	main()
