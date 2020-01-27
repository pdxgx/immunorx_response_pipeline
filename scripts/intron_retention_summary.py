#!/usr/bin/env python

from __future__ import print_function
from collections import defaultdict
from numpy import median
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", "-m",
                        type=str,
                        required=True,
                        help="path to manifest file with tumor-normal pair info")
    parser.add_argument("--output-file", "-f",
                        type=str,
                        required=True,
                        help="path to output directory")
    parser.add_argument("--melanocyte-outliers", "-o",
                        type=str,
                        required=True,
                        help="path to file with melanocyte outlier info")
    parser.add_argument("--tumor-outliers", "-t",
                        type=str,
                        required=True,
                        help="path to file with tumor outlier info")
    args = parser.parse_args()

	# Store patient info
	patients = defaultdict(set)
	with open(os.path.abspath(args.manifest)) as f:
		for line in f:
			tokens = line.strip().split('\t')
			if tokens[2] != 'NA' and tokens[3] != 'NA':
				patients[tokens[0]].add((tokens[2], tokens[3]))

	# Store melanocyte outliers
	melanocyte_introns = set()
	with open(os.path.abspath(args.melanocyte_outliers)) as f:
		melanocyte_samples = f.readline().strip().split('\t')[1:]
		for line in f:
			tokens = line.strip().split('\t')
			if '1' in tokens[1:]:
				melanocyte_introns.add(tokens[0])

	# Write filtered intron outliers
	with open(os.path.abspath(args.tumor_outliers)) as f:
		with open(os.path.abspath(args.output_file), 'w') as o:
			header = f.readline().strip().split('\t')
			print('\t'.join(header), file=o)
			for line in f:
				tokens = line.strip().split('\t')
				if tokens[0] not in melanocyte_introns:
					if '1' in tokens[1:]:
						print(line.strip(), file=o)
						for i in range(1, len(tokens)):
							if tokens[i] == '1':
								intron_dict[header[i]].add(tokens[0])

