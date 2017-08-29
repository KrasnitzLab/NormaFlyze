#!/usr/bin/env python

import sys
import time
import random
import math
import os
import string

def main():

	fastqfilename1 = sys.argv[1]
	fastqfilename2 = sys.argv[2]
	samfilename = sys.argv[3]

	#####  READ FASTQ AND BUILD VT INDEX
	
	FQ1 = open(fastqfilename1, "r")
	FQ2 = open(fastqfilename2, "r")
	
	counter = 0
	prevReadId = ""
	prevVt1 = ""
	prevVt2 = ""
	vt = dict()
	
	for x in FQ1:
		
		x2 = FQ2.readline()
		aline = x.rstrip()
		aline2 = x2.rstrip()
		
		if counter % 4 == 0:
			#exclude the @ at beginning and the /1 or /2 at the end for paired reads
			readId = aline[1:-2]
			prevVt = ""
		
		if counter % 4 == 1:
			vt1 = aline[3:12]
			vt2 = aline2[3:12]
			vt[readId] = [vt1, vt2]

		if counter % 4 == 2:
			pass

		if counter % 4 == 3:
			pass		
				
		counter += 1
	
	FQ1.close()
	FQ2.close()

	#####  READ SAM AND ADD VT INFO
	
	SAM = open(samfilename, "r")
	
	inHeader = True
	
	for x in SAM:
		if inHeader:
			if x[0] != "@":
				inHeader = False
			else:
				print x.rstrip()
				continue
			
		pair_line = SAM.next()
		arow = x.rstrip().split("\t")
		pairrow = pair_line.rstrip().split("\t")

		thisReadId = arow[0]  ### the id of both the forward and reverse strand
		
		if not thisReadId == pairrow[0]:
			print "pairs not in order"

		flag1 = arow[1]  		###  the flag for alignment
		flag2 = pairrow[1]

		thisVt = "NA"
		if vt.has_key(thisReadId):
			if flag1 == "99" and flag2 == "147": 				
				thisVt = vt[thisReadId][0] + "." + vt[thisReadId][1]
			elif flag1 == "83" and flag2 == "163":	
				thisVt = vt[thisReadId][1] + "." + vt[thisReadId][0]
			else:
				pass
			
		arow.append("XV:Z:" + thisVt)
		pairrow.append("XV:Z:" + thisVt)
		print "\t".join(arow)
		print "\t".join(pairrow)
	
	SAM.close()
	
	
if __name__ == "__main__":
	main()
