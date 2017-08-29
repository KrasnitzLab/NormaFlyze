#!/usr/bin/env python

#arguments: directory of files to be searched; length of reads to be found on either side of the string specified; the string to be found
#output: a text file containing reads of length (argv[3]) each on a new line

import re
import sys
import random
import string
import os

def main():

	directory = sys.argv[1]
	length = int(sys.argv[2])
	site = sys.argv[3]
	outfile = str(site) + "_fragments_" + str(length) + "bp_6_28_raw.txt"
	length_file = "lengths.txt"
	f = open(outfile, "w+")
	lf = open(length_file, "w+")

	for file in os.listdir(directory):
		infile = open(os.path.abspath(directory + "/"+ file))

		header = infile.readline()
		seq = infile.read()
		infile.close()
		nospace = ''.join(seq.split())
		y = nospace.upper()
		
		thisStart = 0
		thisEnd = len(y)
		thisIndex = 0
		seg = ""

		listindex = []
		i = y.find(site, thisStart)
		while i >= 0:
			listindex.append(i)
			i = y.find(site, i + 4)

		while thisIndex < len(listindex):
			curr_loc = listindex[thisIndex]

			#if first location is near the beginning
			if thisIndex == 0:

				if curr_loc+4 - length < 0:
					seg = y[0 : curr_loc+4]
					f.write(seg+ "\n")
					lf.write(str(len(seg)) + "\n")

				else:
					seg = y[curr_loc+4 - length : curr_loc+4]
					f.write( seg + "\n")
					lf.write(str(len(seg)) + "\n")

				if len(listindex) > 1:
				
					next_loc = listindex[1]

					if next_loc - curr_loc <= length:
						seg = y[curr_loc+4 : next_loc+4]
						f.write(seg + "\n")
						lf.write(str(len(seg)) + "\n")

					else:
						seg = y[curr_loc+4 : curr_loc+4 + length]
						f.write(seg + "\n")
						lf.write(str(len(seg)) + "\n")
				else:

					if curr_loc+4+length > thisEnd:
						seg = y[curr_loc+4 :]
						f.write( seg + "\n")
						lf.write(str(len(seg)) + "\n")
					else:
						seg = y[curr_loc+4 : curr_loc+4+length]
						f.write(seg + "\n")
						lf.write(str(len(seg)) + "\n")


			#if last location is near end
			elif thisIndex+1 == len(listindex):

				prev_loc = listindex[thisIndex - 1]

				if curr_loc+4+length > thisEnd:
					seg = y[curr_loc+4 :]
					f.write( seg + "\n")
					lf.write(str(len(seg)) + "\n")
				else:
					seg = y[curr_loc+4 : curr_loc+4+length]
					f.write(seg + "\n")
					lf.write(str(len(seg)) + "\n")

				if not curr_loc - prev_loc <= length:
					seg = y[curr_loc+4-length : curr_loc+4]
					f.write(seg + "\n")
					lf.write(str(len(seg)) + "\n")

			#in between locations
			else:
				prev_loc = listindex[thisIndex - 1]
				next_loc = listindex[thisIndex + 1]

				if not curr_loc - prev_loc <= length:
					seg = y[curr_loc+4-length : curr_loc+4]
					f.write(seg + "\n")
					lf.write(str(len(seg)) + "\n")

				if next_loc - curr_loc <= length:
					seg = y[curr_loc+4 : next_loc+4]
					f.write(seg + "\n")
					lf.write(str(len(seg)) + "\n")

				else:
					seg = y[curr_loc+4 : curr_loc+4 + length]
					f.write(seg + "\n")
					lf.write(str(len(seg)) + "\n")

			thisIndex = thisIndex + 1
			
	f.close()
	lf.close()


if __name__ == "__main__":
	main()
