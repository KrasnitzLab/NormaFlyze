#!/usr/bin/env python

#arguments: directory of files to be searched; length of reads to be found on either side of the string specified; the string to be found
#output: 2 text files formatted as FASTQ paired end read files; 1 text file containing fragment lengths

import re
import sys
import random
import string
import os
import time

# args: [seq]: a string representing a DNA sequence
# return: the reverse complement of [seq]
def get_rev_comp(seq):

	base_pairs = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	rev = "".join([base_pairs[base] for base in reversed(seq)])
	return rev
	#from: http://crazyhottommy.blogspot.com/2013/10/python-code-for-getting-reverse.html

def main():

	directory = sys.argv[1]
	length = int(sys.argv[2])
	site = sys.argv[3]
	outfile1 = str(site) + "_paired_" + str(length) + "bp_" + time.strftime("%m_%d") +  "_fastq_1.fq"
	outfile2 = str(site) + "_paired_" + str(length) + "bp_" + time.strftime("%m_%d") +  "_fastq_2.fq"
	outfile3 = str(site) + "_paired_" + str(length) + "bp_" + time.strftime("%m_%d") + "_frag_lengths.txt"

	fow = open(outfile1, "w+")
	rev = open(outfile2, "w+")
	fl = open(outfile3,"w+")

	#length of string being searched
	sl = len(site)

	for file in os.listdir(directory):
		infile = open(os.path.abspath(directory + "/"+ file))

		header = infile.readline()
		header = header.replace(">", "")
		header = header.strip()

		seq = infile.read()
		infile.close()
		nospace = ''.join(seq.split())
		y = nospace.upper()
		
		thisStart = 0
		thisEnd = len(y)
		thisIndex = 0
		seg = ""
		revcomp = ""
		prev = 0

		#create list containing all locations of "site" in this file
		listindex = []
		i = y.find(site, thisStart)
		while i >= 0:
			fl.write(str(i - prev) + "\n")
			prev = i
			listindex.append(i)
			i = y.find(site, i + sl)

		#generate the reads
		while thisIndex < len(listindex):
			
			curr_loc = listindex[thisIndex] + sl
			prev_loc = 0
			next_loc = thisEnd

			if not thisIndex == 0:
				prev_loc = listindex[thisIndex - 1] + sl

			if not thisIndex+1 == len(listindex):
				next_loc = listindex[thisIndex + 1] + sl


			if curr_loc - prev_loc <= length:
				seg = y[prev_loc : curr_loc]
				fow.write("@" + header + "." + str(prev_loc) + "." + str(curr_loc) + "/1"+ "\n")
				fow.write(seg + "\n")
				fow.write("+" +"\n" + "B" * len(seg) + "\n")

				revcomp = get_rev_comp(seg)
				rev.write("@" + header + "." + str(prev_loc) + "." + str(curr_loc) + "/2"+ "\n")
				rev.write(revcomp + "\n")
				rev.write("+" +"\n" + "B" * len(revcomp) + "\n")

			else:
				seg = y[prev_loc : prev_loc+length]
				fow.write("@" + header + "." + str(prev_loc) + "." + str(curr_loc) + "/1"+ "\n")
				fow.write(seg + "\n")
				fow.write("+" +"\n" + "B" * len(seg) + "\n")

				revcomp = y[curr_loc-length : curr_loc]
				revcomp = get_rev_comp(revcomp)
				rev.write("@" + header + "." + str(prev_loc) + "." + str(curr_loc) + "/2"+ "\n")
				rev.write(revcomp + "\n")
				rev.write("+" +"\n" + "B" * len(revcomp) + "\n")

			#the last match
			if thisIndex+1 == len(listindex):

				if curr_loc == thisEnd:
					break

				if next_loc - curr_loc <= length:
					seg = y[curr_loc : next_loc]
					fow.write("@" + header + "." + str(curr_loc) + "." + str(next_loc) + "/1"+ "\n")
					fow.write(seg + "\n")
					fow.write("+" +"\n" + "B" * len(seg) + "\n")

					revcomp = get_rev_comp(seg)
					rev.write("@" + header + "." + str(curr_loc) + "." + str(next_loc) + "/2"+ "\n")
					rev.write(revcomp + "\n")
					rev.write("+" +"\n" + "B" * len(revcomp) + "\n")

				else:
					seg = y[curr_loc : curr_loc+length]
					fow.write("@" + header + "." + str(curr_loc) + "." + str(thisEnd) + "/1"+ "\n")
					fow.write(seg + "\n")
					fow.write("+" +"\n" + "B" * len(seg) + "\n")

					revcomp = y[thisEnd-length :]
					revcomp = get_rev_comp(revcomp)
					rev.write("@" + header + "." + str(curr_loc) + "." + str(thisEnd) + "/2"+ "\n")
					rev.write(revcomp + "\n")
					rev.write("+" +"\n" + "B" * len(revcomp) + "\n")

			thisIndex = thisIndex + 1
			
	fow.close()
	rev.close()
	fl.close()


if __name__ == "__main__":
	main()
