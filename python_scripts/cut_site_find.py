#!/usr/bin/env python

import re
import sys
import random
import string
import os
import cPickle


def main():

	#takes the name of a directory and a string representing the restriction enzyme site
	#the directory should contain the desired chromosome .fa files
	directory = sys.argv[1]
	site = sys.argv[2]
	outfile = str(site) + "_indices.pkl"

	f = open(outfile,"a+b")
	empty = {}
	cPickle.dump(empty, f)
	f.close()

	f = open(outfile, "rb")
	dicto = cPickle.load(f)
	f.close()

	for file in os.listdir(directory):
		print("in here 1")
		infile = open(os.path.abspath(directory + "/"+ file))
		seq = []
		x = ""
		y = ""
		infile.readline()

		#removes all white space and makes all uppercase
		for x in infile:
			seq.append(x.rstrip())
		x = "".join(seq)
		y = x.upper()
		
		thisStart = 0
		thisEnd = len(y)
		thisIndex = 0
		dicto[file] = []

		while thisIndex != -1:
			("in here 2")
			thisIndex = y.find(site, thisStart, thisEnd)
			if thisIndex == -1:
				break
			dicto[file].append(thisIndex)
			thisStart = thisIndex + len(site)

		f = open(outfile, "r+b")
		cPickle.dump(dicto, f)
		f.close()
			
		infile.close()


if __name__ == "__main__":
	main()
