#!/usr/bin/env python

import argparse
import gzip
import os
import pickle
from numpy import percentile

def create_map_dict(id_map, disease_types):
	''' Creates dictionary mapping TCGA CGHubAnalysis UUIDs to TCGA aliquot 
		barcodes for specified disease types

		id_map: path to TCGA_ID_MAP.csv file
		disease_types: list of TCGA disease types to include

		Return value: dictionary, where key is TCGA CGHubAnalysis UUID and 
					  value is TCGA aliquot barcode
	'''
	# Initialize dictionary
	map_dict = {}
	with open(id_map) as f:
		f.readline()
		# Iterate over non-header lines
		for line in f:
			tokens = line.strip().split(',')
			if tokens[3] in disease_types:
				# Relevant disease type, add UUID/barcode to dict
				map_dict[tokens[0]] = tokens[1]
	# Return dictionary
	return map_dict

def create_tpm_dict(input_dir, disease_types, map_dict):
	''' Creates a dictionary storing tumor and normal TPM values for
		different transcripts in different TCGA disease types

		input_dir: path to directory containing TCGA TPM files
		disease_types: list of TCGA disease types to include
		map_dict: dictinoary linking TCGA CGHubAnalysis UUIDs to 
				  TCGA aliquot barcodes (from create_map_dict())

		Return value: nested dictionary where keys are TCGA disease types
					  and values are dictionaries, where keys are transcript
					  IDs and values are nested lists of 
					  [tumor TPM values, normal TMP values]
	'''
	# Initialize dictionary
	tpm_dict = {}
	for disease in disease_types:
		# Create sub-dictionary for each disease
		tpm_dict[disease] = {}
		tpm_file = os.path.join(input_dir, ''.join(['TCGA_', disease, '_tpm.tsv']))
		with open(tpm_file) as f:
			# Process header to determine whether samples are tumor, normal, or control
			header = f.readline().strip().split('\t')
			sample_type_list = []
			for uuid in header:
				# Get barcode from uuid
				barcode = map_dict[uuid].split('-')
				if barcode[3].startswith('0'):
					# sample is tumor
					sample_type_list.append('tumor')
				elif int(barcode[3][0:2]) < 20:
					# sample is normal
					sample_type_list.append('normal')
				else:
					# sample is control
					sample_type_list.append('control')
			# Process TPM data for each transcript
			for line in f:
				tokens = line.strip().split('\t')
				# Get transcript ID
				transcript = tokens[0].split('|')[0]
				tpm_dict[disease][transcript] = []
				# Iterate over TPM values
				for i in range(1, len(tokens)):
					if sample_type_list[i-1] == 'tumor':
						# Value belongs to tumor sample
						tpm_dict[disease][transcript].append(float(tokens[i]))
	# Return dictionary
	return tpm_dict

def create_expression_dict(tpm_dict):
	''' Creates dictionary indicating whether a transcript is expressed if different 
		TCGA disease types

		tpm_dict: nested dictionary where keys are TCGA disease types
				  and values are dictionaries, where keys are transcript
				  IDs and values are nested lists of
				  [tumor TPM values, normal TMP values] (from create_tpm_dict())

		Return value: dictionary where keys are TCGA disease types and values 
					  are sets of transcript IDs for transcripts that are expressed 
					  at a TPM of at least 1 for the 75th percentile expression 
					  value for that transcript among tumor samples
	'''
	# Initialize dictionary
	expression_dict = {}
	for disease in tpm_dict:
		# Create subdictionary for each disease
		expression_dict[disease] = set()
		# Iterate over transcripts
		for transcript in tpm_dict[disease]:
			quantile75 = percentile(tpm_dict[disease][transcript], 0.75)
			if quantile75 >= 1:
				expression_dict[disease].add(transcript)
	# Return dictionary
	return expression_dict

# Set up command line parameters
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-dir', type=str, required=True,
					help='input directory containing TCGA ID map and TPM files'
				)
parser.add_argument('-s', '--store-intermediate', required=False, action='store_true',
					help='whether to store intermediate dictionaries'
				)
args = parser.parse_args()

# Establish diseases of interest and map file path
# Need comparators for melanoma, NSCLC, colon, endometrial, thyroid, prostate, renal cell carcinoma
# Using SKCM, LUAD/LUSC, COAD, UCEC, THCA, PRAD, KIRC
tcga_diseases = ['SKCM', 'LUAD', 'LUSC', 'COAD', 'UCEC', 'THCA', 'PRAD', 'KIRC']
map_file = os.path.join(args.input_dir, 'TCGA_ID_MAP.csv')

# Create CGHubAnalysis UUID to TCGA aliquot barcode dict
map_dict = create_map_dict(map_file, tcga_diseases)
if args.store_intermediate:
	map_path = os.path.join(args.input_dir, 'TCGA_ID_map.pickle')
	with open(map_path, 'wb') as p:
		pickle.dump(map_dict, p, protocol=2)

# Create disease type to transcript ID to expression lists dict
tpm_dict = create_tpm_dict(args.input_dir, tcga_diseases, map_dict)
if args.store_intermediate:
	tpm_path = os.path.join(args.input_dir, 'TCGA_TPMs.pickle')
	bytes_out = pickle.dumps(tpm_dict, protocol=2)
	max_bytes = 2**31 - 1
	with open(tpm_path, 'wb') as p:
	    for idx in range(0, len(bytes_out), max_bytes):
	        p.write(bytes_out[idx:idx+max_bytes])

# Create disease type to expressed transcripts dict
expression_dict = create_expression_dict(tpm_dict)
expression_path = os.path.join(args.input_dir, 'TCGA_expression.pickle')
with open(expression_path, 'wb') as p:
	pickle.dump(expression_dict, p, protocol=2)
