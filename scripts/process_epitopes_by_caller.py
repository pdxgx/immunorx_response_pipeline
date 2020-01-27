#!/usr/bin/env python

from __future__ import print_function
import glob
import os
from collections import defaultdict
from numpy import median

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-file", "-p",
                        type=str,
                        required=True,
                        help="path to manifest file with tumor-normal pairs")
    parser.add_argument("--output-dir", "-o",
                        type=str,
                        required=True,
                        help="path to output directory")
    parser.add_argument("--caller-epitope-dir", "-e",
                        type=str,
                        required=True,
                        help="path to directory containing per-caller neoepitope predictions")
    parser.add_argument("--consensus-epitope-dir", "-c",
                        type=str,
                        help="path to directory containing consensus-variant neoepitope calls")
    args = parser.parse_args()

	summary_dict = {}
	patients = defaultdict(set)
	with open(pair_file) as p:
		for line in p:
			tokens = line.strip().split('\t')
			if tokens[2] != 'NA':
				# Get patient IDs
				patient_id = tokens[0]
				tumor = tokens[2]
				patients[patient_id].add(tumor)
	for patient in patients:
		epitope_dict = {}
		for tumor in patients[patient]:
			epitope_dict[tumor] = [set(), set(), set(), set(), set(), set(), set()]
			# Process MuSe data
			muse_file = ''.join([args.caller_epitope_dir, patient, '.', tumor, '.neoepiscope.muse.out'])
			with open(muse_file, 'r') as f:
				f.readline()
				for line in f:
					tokens = line.strip().split('\t')
					epitope_dict[tumor][0].add(tokens[0])
					if tokens[0] not in summary_dict:
						summary_dict[tokens[0]] = [0, 0, 0, 0, 0, 0, 0]
					summary_dict[tokens[0]][0] += 1
			# Process MuTect data
			mutect_file = ''.join([args.caller_epitope_dir, patient, '.', tumor, '.neoepiscope.mutect.out'])
			with open(mutect_file, 'r') as f:
				f.readline()
				for line in f:
					tokens = line.strip().split('\t')
					epitope_dict[tumor][1].add(tokens[0])
					if tokens[0] not in summary_dict:
						summary_dict[tokens[0]] = [0, 0, 0, 0, 0, 0, 0]
					summary_dict[tokens[0]][1] += 1
			# Process Pindel data
			if tumor not in ['SRR6431469', 'SRR6431473']: # Skip tumors with no Pindel variants
				pindel_file = ''.join([args.caller_epitope_dir, patient, '.', tumor, '.neoepiscope.pindel.out'])
				with open(pindel_file, 'r') as f:
					f.readline()
					for line in f:
						tokens = line.strip().split('\t')
						epitope_dict[tumor][2].add(tokens[0])
						if tokens[0] not in summary_dict:
							summary_dict[tokens[0]] = [0, 0, 0, 0, 0, 0, 0]
						summary_dict[tokens[0]][2] += 1
			# Process RADIA data
			if tumor not in ['SRR3341247', 'SRR3341245']: # Skip tumors w/ no RADIA variants
				radia_file = ''.join([args.caller_epitope_dir, patient, '.', tumor, '.neoepiscope.radia.out'])
				with open(radia_file, 'r') as f:
					f.readline()
					for line in f:
						tokens = line.strip().split('\t')
						epitope_dict[tumor][3].add(tokens[0])
						if tokens[0] not in summary_dict:
							summary_dict[tokens[0]] = [0, 0, 0, 0, 0, 0, 0]
						summary_dict[tokens[0]][3] += 1
			# process SomaticSniper data
			if tumor not in ['SRR6431481']: # Skip tumors with no SomaticSniper variants
				somatic_sniper_file = ''.join([args.caller_epitope_dir, patient, '.', tumor, '.neoepiscope.somatic_sniper.out'])
				with open(somatic_sniper_file, 'r') as f:
					f.readline()
					for line in f:
						tokens = line.strip().split('\t')
						epitope_dict[tumor][4].add(tokens[0])
						if tokens[0] not in summary_dict:
							summary_dict[tokens[0]] = [0, 0, 0, 0, 0, 0, 0]
						summary_dict[tokens[0]][4] += 1
			# process VarScan data
			varscan_file = ''.join([args.caller_epitope_dir, patient, '.', tumor, '.neoepiscope.varscan.out'])
			with open(varscan_file, 'r') as f:
				f.readline()
				for line in f:
					tokens = line.strip().split('\t')
					epitope_dict[tumor][5].add(tokens[0])
					if tokens[0] not in summary_dict:
						summary_dict[tokens[0]] = [0, 0, 0, 0, 0, 0, 0]
					summary_dict[tokens[0]][5] += 1
			# process consensus data
			consensus_file = ''.join([args.consensus_epitope_dir, patient, '.', tumor, '.neoepiscope.somatic.out'])
			with open(consensus_file, 'r') as f:
				f.readline()
				for line in f:
					tokens = line.strip().split('\t')
					if len(tokens[0]) >= 8 and len(tokens[0]) <= 11:
						epitope_dict[tumor][6].add(tokens[0])
						if tokens[0] not in summary_dict:
							summary_dict[tokens[0]] = [0, 0, 0, 0, 0, 0, 0]
						summary_dict[tokens[0]][6] += 1

	# Find number of peptides for each caller combination
	combo_dict = {}
	for epitope in summary_dict:
		combo_dict[','.join(out_line[1:])].add(epitope)
	
	# Write combinations to output file
	upset_file = ''.join([args.output_dir, 'epitope_upset_table.tsv'])
	with open(upset_file, 'w') as f:
		for combo in combo_dict:
			print('\t'.join([combo, len(combo_dict[combo])]), file=f)
