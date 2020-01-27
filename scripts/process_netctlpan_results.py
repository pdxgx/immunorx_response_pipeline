#!/usr/bin/env python

from __future__ import print_function
from collections import defaultdict
from numpy import median
import argparse
import glob
import os
import pickle

if __name__ == "__main__":

	# Parse command line options
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--manifest', type=str, required=True,
						help='path to tumor-normal pair manifest'
	)
	parser.add_argument('-i', '--input-dir', type=str, required=True,
						help='path to input directory with neCTLpan result dictionaries'
	)
	parser.add_argument('-o', '--output-file', type=str, required=True,
						help='path to output file'
	)
	args = parser.parse_args()

	# Iterate through manifest to process scores for each sample
	tumor_dict = defaultdict(int)
	with open(os.path.abspath(args.manifest)) as f:
		for line in f:
			tokens = line.strip().split('\t')
			# Grab all dictionaries for the sample
			dict_wildcard = os.path.join(
											os.path.abspath(args.input_dir), 
											'.'.join([tokens[0], tokens[2], '*', 'pickle'])
			)
			dictionaries = glob.glob(dict_wildcard)
			# Process scores to find best rank for each epitope
			score_dict = defaultdict(lambda:100.0)
			for dic in dictionaries:
				with open(dic, 'rb') as p:
					d = pickle.load(p)
				for entry in d:
					for score_set in d[entry]:
						rank = float(score_set[4])
						# Score retained is minumum rank for each peptide
						score_dict[entry[1]] = min(rank, score_dict[entry[1]])
			for peptide in score_dict:
				if score_dict[peptide] < 1:
					tumor_dict[(tokens[0], tokens[2])] += 1

	# Get median values for multi-sample patients
	multisample = set([x[0] for x in tumor_dict if len([y for y in tumor_dict if x[0] in y]) > 1])
	for patient in multisample:
		# Extract relevant keys/entries
		try:
			# Separate patients from Roh/Amaria cohorts (#s as pat. IDs)
			p = int(patient)
			relevant_keys = [x for x in tumor_dict if x[0] == patient and (
				  x[1][-1] in ['A', 'B', 'C', 'D', 'E'] or x[1] in [''.join([patient, 'D1']), ''.join([patient, 'D2'])]
				)
			]
			if len(relevant_keys) == 1:
				continue
		except ValueError:
			relevant_keys = [x for x in tumor_dict if x[0] == patient]
		relevant_entries = [tumor_dict[x] for x in relevant_keys]
		# Set up new combined key/entry
		new_key = (patient, ';'.join(sorted([x[1] for x in relevant_keys])))
		new_entry = median(relevant_entries)
		tumor_dict[new_key] = new_entry
		for key in relevant_keys:
			del tumor_dict[key]

	# Write output file
	with open(os.path.abspath(args.output_file), 'w') as f:
		header = ['Patient', 'Tumor_ID', 'NetCTLpan_epitopes']
		print('\t'.join(header), file=f)
		for patient in tumor_dict:
			out_line = [patient[0], patient[1], str(tumor_dict[patient])]
			print('\t'.join(out_line), file=f)

