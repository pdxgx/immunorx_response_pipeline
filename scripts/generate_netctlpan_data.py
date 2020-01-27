#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import pickle
import subprocess
import sys
from Bio import SeqIO
from collections import defaultdict

def generate_fasta_entries(input_fasta, output_fasta):
	""" Converts neoepiscope output FASTAs new FASTA file with IDs as headers,
		storing IDs as a dictionary for mapping

		input_fasta: path to input neoepiscope FASTAs
		output_fasta: path to output fasta with alternate headers

		Return values: 
			mapper_dict: dictionary linking FASTA header IDs as keys to
						 original fasta headers as values
	"""
	# Set up dictionaries and sequence counter
	fasta_dict = {}
	mapper_dict = {}
	i = 0

	print('Parsing FASTA records', file=sys.stderr)
		
	# Iterate through all records in the FASTA file to extract data
	for record in SeqIO.parse(input_fasta, 'fasta'):
		
		# Extract sequence header and sequence
		header = str(record.id)
		sequence = str(record.seq)

		# Skip sequences that are too short for netCTLpan
		if len(sequence) < 8:
			continue

		# Create new header and save entry to dictionary
		new_header = ' '.join(['>', str(i)])
		fasta_dict[new_header] = sequence
		mapper_dict[str(i)] = header

		# Increment counter
		i += 1

	print(' '.join(['Stored', str(i), 'FASTA sequences.']), file=sys.stderr)

	# Write output fasta
	with open(output_fasta, 'w') as f:
		for header in fasta_dict:
			print(header, file=f)
			print(fasta_dict[header], file=f)

	# Return dictionaries
	return fasta_dict, mapper_dict

def run_netCTLpan(fasta, output_path, netCTLpan, n, allele):
	""" Runs netCTLpan on input FASTA file
		
		fasta: path to FASTA file with alternate headers (written by 
			   generate_fasta_entries() function)
		output_path: path to write netCTLpan output
		netCTLpan: path to netCTLpan executable
		n: peptide_length
		allele: MHC allele

		No return value.
	"""

	print('Running netCTLpan', file=sys.stderr)

	# Set netCTLpan command and run
	command = [
				netCTLpan, '-v', '-f', fasta, '-xls', '-xlsfile', output_path,
				'-l', n, '-a', allele
			  ]
	subprocess.check_call(command)

def extract_neoepitopes(neoepiscope_results):
	""" Extract neoepitopes to store for each sample/transcript

		neoepiscope_results: path to neoepiscope results
		
		Return value: dictionary with transcript IDs as keys with sets of 
					  relevant neoepitope sequences as values
	"""
	# Set up dictionary to store neoepitope information
	neoepitope_dict = defaultdict(set)

	print('Extracting neoepitope data', file=sys.stderr)

	# Process data
	with open(neoepiscope_results) as f:
		f.readline()
		f.readline()
		for line in f:
			tokens = line.strip().split('\t')
			# Extract transcripts to iterate through
			transcripts = tokens[9].split(';')
			for tx in transcripts:
				# Store neoepitope for each relevant transcript
				neoepitope_dict[tx].add(tokens[0])

	# Return neoepitope dictionary
	return neoepitope_dict

def process_netCTLpan_output(
								netCTLpan_output, mapping_dictionary,
								neoepitope_dictionary
	):
	""" Processes netCTL output to retain data only for neoepitope sequences,
		storing them in a dictionary
	
		netCTLpan_output: path to netCTLpan output
		mapping_dictionary: dictionary linking FASTA header IDs as keys to 
							original fasta headers as values
		neoepitope_dictionary: dictionary with transcript IDs keys with sets of 
							   relevant neoepitope sequences as values

		Return value: dictionary linking tuple of (transcript ID, peptide) as 
					  keys to set of tuples of (TAP score, cleavage score) as
					  values
	"""
	print('Processing netCTLpan output', file=sys.stderr)
	netCTL_score_dict = defaultdict(set)
	# Process file to extract data
	with open(netCTLpan_output) as f:
		for line in f:
			if line[0] != 'N':
				# Process peptide result
				tokens = line.strip().split('\t')
				# [N, Sequence Name, Peptide, Allele, MHC, TAP, Cle, Comb, %Rank]
				
				# Grab sequence identifier, removing initial underscore
				# Split by middle undescore to isolate just the ID
				# (large proteins will have been split into more entries)
				identifier = tokens[1].split('_')[1]
				
				# Use identifier to grab sample/transcript information
				original_header = mapping_dictionary[identifier]
				transcript = original_header.lstrip('>').split('_')[0]

				# Check if peptide is a neoepitope for that transcript
				peptide = tokens[2]
				if peptide in neoepitope_dictionary[transcript]:
					# Store MHC rank, TAP score, cleavage score, combined score, % rank
					scores = (tokens[4], tokens[5], tokens[6], tokens[7], tokens[8])
					netCTL_score_dict[(transcript, peptide)].add(scores)

	print("Done processing output", file=sys.stderr)

	# Return dictionary
	return netCTL_score_dict

if __name__ == "__main__":

	# Parse command line options
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', type=str, required=True,
						help='path to neoepiscope output file'
	)
	parser.add_argument('-o', '--output', type=str, required=True,
						help='path to write FASTAs and netCTL output files'
	)
	parser.add_argument('-e', '--executable', type=str, required=True,
						help='path to netCTLpan executable'
	)
	parser.add_argument('-a', '--allele', type=str, required=True,
						help='HLA allele'
	)
	parser.add_argument('-s', '--size', type=str, required=True,
						help='peptide size'
	)
	args = parser.parse_args()

	# Get absolute paths to input/output directories/files
	input_file = os.path.abspath(args.input)
	input_fasta = '.'.join([input_file, 'fasta'])
	sample_id = os.path.basename(input_file).replace('.neoepiscope.comprehensive.out', '')
	output_dir = os.path.abspath(args.output)
	output_fasta = os.path.join(
									output_dir, 
									'.'.join([sample_id, args.size, args.allele, 'netCTL.fasta'])
	)
	output_dictionary = os.path.join(
										output_dir, 
										'.'.join([sample_id, args.size, args.allele, 'pickle'])
	)
	netCTLpan = os.path.abspath(args.executable)

	# Generate FASTA dictionaries
	fasta_entries, mapper = generate_fasta_entries(input_fasta, output_fasta)
	
	# Extract neoepitopes
	neoepitopes = extract_neoepitopes(input_file)

	# Run netCTLpan
	complete_score_dict = {}
	netCTLpan_output = os.path.join(
									output_dir, 
									'.'.join(
												[
													sample_id, 
													args.size,
													args.allele, 
													'netCTL.out'
												]
									)
	)
	run_netCTLpan(output_fasta, netCTLpan_output, netCTLpan, args.size, args.allele)

	# Process netCTLpan data
	netCTLpan_scores = process_netCTLpan_output(
												netCTLpan_output, 
												mapper, 
												neoepitopes
	)
	complete_score_dict.update(netCTLpan_scores)

	print('Saving results', file=sys.stderr)

	# Store results in pickled dictionary
	with open(output_dictionary, 'wb') as p:
		pickle.dump(complete_score_dict, p)

	print('Done!', file=sys.stderr)
