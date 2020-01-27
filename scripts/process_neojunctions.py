#!/usr/bin/env python

from __future__ import print_function
from collections import defaultdict
from numpy import median
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--jx-dir', '-j', type=str,required=True, help='path to directory containing junction files')
    parser.add_argument('--output-dir', '-o', type=str,required=True, help='path to output directory')
    parser.add_argument("--manifest", "-m", type=str, required=True, help="path to tumor-normal pair manifest file")
    args = parser.parse_args()

    # Get patient info
	patients = defaultdict(set)
	with open(os.path.abspath(args.manifest)) as f:
		for line in f:
			tokens = line.strip().split('\t')
			if tokens[2] != 'NA' and tokens[3] != 'NA':
				patients[tokens[0]].add((tokens[2], tokens[3]))

	# Tally junctions and write output
	out_file = os.path.join(os.path.abspath(args.output_dir), 'patient_jx_burdens.tsv')
	with open(out_file, 'w') as o:
		header = ['Patient', 'Tumor_ID', 'Jx_burden']
		print('\t'.join(header), file=o)
		for patient in patients:
			jx_counts = []
			for tumor in patients[patient]:
				jx_file = glob.glob(os.path.join(os.path.abspath(args.jx_dir), ''.join([tumor[1], '_queryresults_*.csv'])))[0]
				jx_count = 0
				with open(jx_file) as j:
					j.readline()
					for line in j:
						tokens = line.strip().split(',')
						if tokens[1] in ['1', '2']:
							# Not GENCODE annotated - store junction
							jx_count += 1
				jx_counts.append(jx_count)
			# Write output
			out_line = [patient, ';'.join(sorted([x[0] for x in patients[patient]])), str(median(jx_counts))]
			print('\t'.join(out_line), file=o)
