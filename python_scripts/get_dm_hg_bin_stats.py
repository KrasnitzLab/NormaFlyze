#!/usr/bin/env python

import sys
import time
from operator import itemgetter
import os
import numpy

def main():

	bins = "/mnt/wigclust1/data/safe/kostic/bin_mapping/hybrid_bin_boundaries_sorted_07_25.txt"
	uber = "/mnt/wigclust1/data/safe/kostic/SNS_data/uber_varbin_count_data.txt"

	numpy.loadtext(uber, skiprows = 1)

	UBER = open(uber, "r")
	UBER.readline()
	for line in UBER:





if __name__ == "__main__":
	main()