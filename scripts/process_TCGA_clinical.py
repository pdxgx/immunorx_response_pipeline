#!/usr/bin/env python

from __future__ import print_function
import argparse
import glob
import os
import pandas as pd
from collections import defaultdict

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--maf-dir', '-m', type=str,required=True, help='path to directory containing MAFs')
    parser.add_argument('--clinical', '-c', type=str,required=True, help='path to file containing clinical_data')
    parser.add_argument("--output-file", "-o", type=str, required=True, help="path to output file")
    args = parser.parse_args()

	# Get mutational burden
	tmb_dict = defaultdict(int)
	maf_files = glob.glob(os.path.join(os.path.abspath(args.maf_dir), '*.maf.txt'))
	for maf in maf_files:
		# Get barcode
		barcode = maf.split('/')[-1].replace('.maf.txt', '').lower()
		df = pd.read_csv(maf, sep='\t', encoding='latin-1')
		# tally mutations
		for index, row in df.iterrows():
			tmb_dict[barcode] += 1

	# Get clinical data
	with open(os.path.abspath(args.clinical)) as f:
		for line in f:
			tokens = line.strip().split('\t')
			if tokens[0] == 'patient.days_to_birth':
				days_to_birth = tokens[1:]
			elif tokens[0] == 'patient.days_to_death':
				days_to_death = tokens[1:]
			elif tokens[0] == 'patient.days_to_last_followup':
				days_to_censor = tokens[1:]
			elif tokens[0] == 'patient.stage_event.pathologic_stage':
				cancer_stage = [x.replace('stage ', '').replace(' nos', '').upper() for x in tokens[1:]]
			elif tokens[0] == 'patient.samples.sample.bcr_sample_barcode':
				barcodes = [x[:-1] for x in tokens[1:]]

	# Write output file
	with open(os.path.abspath(args.output_file), 'w') as o:
		header = ['Patient', 'Stage', 'TMB', 'Overall_survival', 'Days_to_last_followup', 'OS_event', 'Censored']
		print('\t'.join(header), file=o)
		for i in range(len(barcodes)):
			if tmb_dict[barcodes[i]] > 0:
				out_line = [barcodes[i].upper(), cancer_stage[i], str(tmb_dict[barcodes[i]])]
				if days_to_death[i] != 'NA':
					out_line.extend([days_to_death[i], 'NA', '1', '0'])
					print('\t'.join(out_line), file=o)
				else:
					if days_to_censor[i] != 'NA':
						if int(days_to_censor[i]) > 0:
							out_line.extend([days_to_censor[i], days_to_censor[i], '0', '1'])
							print('\t'.join(out_line), file=o)
