#!/usr/bin/env python

import re
import sys
import random
import string
import os

def main():
	
	infile = open(sys.argv[1], "r")
	path_out = sys.argv[2]

	opened = False

	for line in infile:
		if not line.find(">") == -1:
			if(opened):
				outfile.close()
			opened = True
			title = line[1:].strip() + "_dm"
			line = ">" + title + "\n"
			outfile = open(path_out + "/" + title + ".fa", "w+")
		outfile.write(line)

	outfile.close()


if __name__ == "__main__":
	main()