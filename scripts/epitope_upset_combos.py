#!/usr/bin/env python

from __future__ import print_function
from collections import defaultdict
import argparse
import os

if __name__ == "__main__":

	# Parse command line options
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input-file', type=str, required=True,
						help='path to input file with caller combinations'
	)
	parser.add_argument('-o', '--output-dir', type=str, required=True,
						help='path to output directory'
	)
	args = parser.parse_args()

	combo_dict = {}

	callers = ['MuSE', 'MuTect', 'Pindel', 'RADIA', 'SomaticSniper', 'VarScan', 'Consensus']

	# Get caller combinations
	with open(os.path.abspath(args.input_file)) as f:
		for line in f:
			tokens = line.strip().split('\t')
			caller_combo = tokens[0].split(',')
			caller_list = []
			for i in range(7):
				if caller_combo[i] == '1':
					caller_list.append(callers[i])
			combo_dict[tuple(caller_list)] = int(tokens[1])

	caller_count = defaultdict(int)

	# Write output files
	with open(os.path.join(os.path.abspath(args.output_dir), 'epitope_upset.tsv'), 'w') as o1:
		other_count = 0
		header = ['overlap', 'tool', 'classification', 'proportion']
		print('\t'.join(header), file=o1)
		with open(os.path.join(os.path.abspath(args.output_dir), 'epitope_upset_counts.tsv'), 'w') as o2:
			header = ['overlap', 'count']
			print('\t'.join(header), file=o2)
			for i in range(1, 8):
				# Sort combos of same caller count by decreasing order
				combo_ids = sorted([x for x in combo_dict if len(x) == i], key=lambda k: combo_dict[k], reverse=True)
				if i == 1:
					# Retain all single-caller combos
					reduced_combo_ids = [x for x in combo_ids]
					skipped_combos = []
				else:
					# Retain only top two combos for other cases
					reduced_combo_ids = combo_ids[0:2]
					skipped_combos = combo_ids[2:]
				# Write output lines for combos
				for combo in reduced_combo_ids:
					for caller in callers:
						if caller in combo:
							caller_count[caller] += combo_dict[combo]
							out_line1 = ['_'.join(combo), caller, 'present', '1']
							out_line2 = ['_'.join(combo), caller, 'absent', '0']
							out_line3 = ['_'.join(combo), caller, 'other', '0']
						else:
							out_line1 = ['_'.join(combo), caller, 'present', '0']
							out_line2 = ['_'.join(combo), caller, 'absent', '1']
							out_line3 = ['_'.join(combo), caller, 'other', '0']
						print('\t'.join(out_line1), file=o1)
						print('\t'.join(out_line2), file=o1)
						print('\t'.join(out_line3), file=o1)
					out_line = ['_'.join(combo), str(combo_dict[combo])]
					print('\t'.join(out_line), file=o2)
				for combo in skipped_combos:
					other_count += combo_dict[combo]
			# Write summary output lines for skipped combos
			for caller in callers:
				out_line1 = ['Other', caller, 'present', '0']
				out_line2 = ['Other', caller, 'absent', '0']
				out_line3 = ['Other', caller, 'other', '1']
				print('\t'.join(out_line1), file=o1)
				print('\t'.join(out_line2), file=o1)
				print('\t'.join(out_line3), file=o1)
			out_line = ['Other', str(other_count)]
			print('\t'.join(out_line), file=o2)
