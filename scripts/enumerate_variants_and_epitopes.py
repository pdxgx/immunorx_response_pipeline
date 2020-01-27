#!/usr/bin/env python

from __future__ import print_function
import argparse
import pandas as pd 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--immunorx-table', '-i', type=str,required=True, help='path to immunotherapy data table')
    parser.add_argument('--jx-file', '-j', type=str,required=True, help='path to file containing junction data')
    parser.add_argument('--ri-file', '-r', type=str,required=True, help='path to file containing retained intron data')
    parser.add_argument("--output-dir", "-o", type=str, required=True, help="path to output directory")
    args = parser.parse_args()

	mut_dict = {}
	ep_dict = {}

	# Get DNA variants/epitopes
	immunorx = pd.read_csv(os.path.abspath(args.immunorx_table), sep='\t')
	for idx, row in immunorx.iterrows():
		patient = row['Patient']
		tumor = row['Tumor_ID']
		total_mutations = int(float(row['Total_mutations']))
		snvs = int(float(row['SNVs']))
		indels = int(float(row['Frameshift_deletions']) + float(row['Frameshift_insertions']) + float(row['Inframe_deletions']) + float(row['Inframe_insertions']))
		epitopes = int(float(row['MHCnuggets_eps']))
		mut_dict[(patient, tumor)] = [total_mutations, snvs, indels]
		ep_dict[(patient, tumor)] = [epitopes]

	# Get junctions
	jx = pd.read_csv(os.path.abspath(args.jx_file), sep='\t')
	for idx, row in jx.iterrows():
		patient = row['Patient']
		tumor = row['Tumor']
		jx = int(float(row['Jx_burden']))
		mut_dict[(patient, tumor)].append(jx)

	# Get retained introns/RI epitopes
	ri = pd.read_csv(os.path.abspath(args.ri_file), sep='\t')
	for idx, row in ri.iterrows():
		patient = row['Patient']
		tumor = row['Tumor_ID']
		ri = int(float(row['Intron_burden']))
		ri_eps = int(float(row['Binding_intron_epitope_burden']))
		mut_dict.append(ri)
		ep_dict.append(ri_eps)

	# Write mutation output file
	with open(os.path.join(os.path.abspath(args.output_dir), 'rna_enumerated_mutations.tsv'), 'w') as o:
		print('\t'.join(['Patient', 'Class', 'Mut_count', 'TVB']), file=o)
		for patient in pat_dict:
			data = pat_dict[patient]
			if len(data) > 3:
				tvb = str(data[0] + data[3] + data[4])
				# Enumerate junctions
				for i in range(data[3]):
					out_line = ['_'.join([patient[0], patient[1]]), 'Jx', str(data[3]), tvb]
					print('\t'.join(out_line), file=o)
				# Enumerate RIs
				for i in range(data[4]):
					out_line = ['_'.join([patient[0], patient[1]]), 'RI', str(data[4]), tvb]
					print('\t'.join(out_line), file=o)
				# Enumerate SNVs
				for i in range(data[1]):
					out_line = ['_'.join([patient[0], patient[1]]), 'SNV', str(data[1]), tvb]
					print('\t'.join(out_line), file=o)
				# Enumerate indels
				for i in range(data[2]):
					out_line = ['_'.join([patient[0], patient[1]]), 'Indel', str(data[2]), tvb]
					print('\t'.join(out_line), file=o)

	# Write epitope output file
	with open(os.path.join(os.path.abspath(args.output_dir), 'rna_enumerated_epitopes.tsv'), 'w') as o:
		print('\t'.join(['Patient', 'Class', 'Mut_count', 'TVB']), file=o)
		for patient in pat_dict:
			data = pat_dict[patient]
			if len(data) > 1:
				tvb = str(data[0] + data[1])
				# Enumerate RI epitopes
				for i in range(data[1]):
					out_line = ['_'.join([patient[0], patient[1]]), 'RI', str(data[1]), tvb]
					print('\t'.join(out_line), file=o)
				# Enumerate DNA epitopes
				for i in range(data[0]):
					out_line = ['_'.join([patient[0], patient[1]]), 'SNV', str(data[0]), tvb]
					print('\t'.join(out_line), file=o)



