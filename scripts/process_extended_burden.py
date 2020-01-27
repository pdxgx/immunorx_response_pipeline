#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
from collections import defaultdict
from numpy import median
from operator import itemgetter

def process_extended_burden_file(burden_file):
	""" Processes a per-patient extended burden file to summarize it

		burden_file: path to per-patient burden file

		Return value: list of summary data for burden file
	"""
	epitope_dict = defaultdict(dict)
	# Establish counts for burden
	ep_by_mismatches = 0.0
	ep_by_anchor = 0.0
	ep_by_alleles = 0.0
	ep_by_expression = 0.0
	ep_by_tcga = 0.0
	ep_by_mismatches_and_alleles = 0.0
	ep_by_mismatches_and_expression = 0.0
	ep_by_mismatches_and_tcga = 0.0
	ep_by_anchor_and_alleles = 0.0
	ep_by_anchor_and_expression = 0.0
	ep_by_anchor_and_tcga = 0.0
	ep_by_alleles_and_expression = 0.0
	ep_by_alleles_and_tcga = 0.0
	ep_by_mismatches_alleles_and_expression = 0.0
	ep_by_mismatches_alleles_and_tcga = 0.0
	ep_by_anchor_alleles_and_expression = 0.0
	ep_by_anchor_alleles_and_tcga = 0.0
	alleles = set()
	alleles_per_mutation = defaultdict(set)
	# Run initial process through file
	with open(burden_file) as f:
		f.readline()
		for line in f:
			tokens = line.strip().split('\t')
			mut_id = tuple(tokens[3].split(','))
			if mut_id in epitope_dict[tokens[2]]:
				if tokens[5:] not in epitope_dict[tokens[2]][mut_id]:
					 epitope_dict[tokens[2]][mut_id].append(tokens[5:])
			else:
				epitope_dict[tokens[2]][mut_id] = [tokens[5:]]
	# Process through epitopes
	for epitope in epitope_dict:
		if len(epitope_dict[epitope]) != 1:
			# Reduce redudant data if applicable
			mutations = list(epitope_dict[epitope].keys())
			mutations.sort(key=itemgetter(0,1))
			mutation_track = []
			mutation_info = []
			for i in range(len(mutations) - 1):
				if epitope_dict[epitope][mutations[i]] == epitope_dict[epitope][mutations[i+1]]:
					if mutations[i] in mutation_track:
						mutation_track.append(mutations[i+1])
					elif mutation_track == []:
						mutation_track = [mutations[i], mutations[i+1]]
						mutation_info = [tokens[5:]]
					else:
						epitope_dict[epitope][tuple([mut for mut in mutation_track])] = mutation_info
						for mut in mutation_track:
							del epitope_dict[epitope][mut]
						mutation_track = [mutations[i], mutations[i+1]]
						mutation_info = [tokens[5:]]
				else:
					if mutation_track != []:
						epitope_dict[epitope][tuple([mut for mut in mutation_track])] = mutation_info
						for mut in mutation_track:
							del epitope_dict[epitope][mut]
						mutation_track = []
						mutation_info = []
		# Process data for each mutation set
		for mutation_set in epitope_dict[epitope]:
			for data in epitope_dict[epitope][mutation_set]:
				# Process allele data
				alleles_presenting = data[8].split(',')
				for hla in alleles_presenting:
					alleles.add(hla)
					alleles_per_mutation[mutation_set].add(hla)
				ep_by_alleles += len(alleles_presenting)
				# Add mismatch info if available
				if data[3] != 'NA':
					ep_by_mismatches += float(data[3])
					ep_by_mismatches_and_alleles += (len(alleles_presenting)*float(data[3]))
					ep_by_anchor += float(data[4])
					ep_by_anchor_and_alleles += (len(alleles_presenting)*float(data[4]))
				# Add TCGA expression info
				ep_by_tcga += float(data[13])
				ep_by_alleles_and_tcga += (len(alleles_presenting)*float(data[13]))
				if data[3] != 'NA' :
					ep_by_mismatches_and_tcga += (float(data[3])*float(data[13]))
					ep_by_mismatches_alleles_and_tcga += (len(alleles_presenting)*float(data[3])*float(data[13]))
					ep_by_anchor_and_tcga += (float(data[4])*float(data[13]))
					ep_by_anchor_alleles_and_tcga += (len(alleles_presenting)*float(data[4])*float(data[13]))
				# Add patient-specific expression info if available
				if data[11] != 'NA':
					ep_by_expression += float(data[11])
					ep_by_alleles_and_expression += (len(alleles_presenting)*float(data[11]))
					if data[3] != 'NA':
						ep_by_mismatches_and_expression += (float(data[3])*float(data[11]))
						ep_by_mismatches_alleles_and_expression += (len(alleles_presenting)*float(data[3])*float(data[11]))
						ep_by_anchor_and_expression += (float(data[4])*float(data[11]))
						ep_by_anchor_alleles_and_expression += (len(alleles_presenting)*float(data[4])*float(data[11]))
	# Summarize and return data
	data_list = [ep_by_mismatches, ep_by_anchor, ep_by_alleles, ep_by_expression, ep_by_tcga, ep_by_mismatches_and_alleles,
				 ep_by_mismatches_and_expression, ep_by_mismatches_and_tcga, ep_by_anchor_and_alleles,
				 ep_by_anchor_and_expression, ep_by_anchor_and_tcga, ep_by_alleles_and_expression,
				 ep_by_alleles_and_tcga, ep_by_mismatches_alleles_and_expression, ep_by_mismatches_alleles_and_tcga,
				 ep_by_anchor_alleles_and_expression, ep_by_anchor_alleles_and_tcga,
				 len(alleles), sum([len(alleles_per_mutation[x]) for x in alleles_per_mutation])/len(alleles_per_mutation)]
	return data_list

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--burden-dir', '-b', type=str,required=True, help='path to directory containing burden files')
    parser.add_argument('--rna-dir', '-r', type=str,required=True, help='path to directory containing RNA bams')
    parser.add_argument("--manifest", "-m", type=str, required=True, help="path to tumor-normal pair manifest file")
    args = parser.parse_args()

	# Studies used
	pairs = ['carreno', 'gao', 'le', 'le2', 'hugo', 'miao', 'miao2', 'rizvi', 'roh', 'snyder', 
			 'graff', 'vanallen', 'zaretsky', 'amaria', 'eroglu']

	mutation_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
	neoepitope_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

	combined_file = os.path.join(os.path.abspath(args.burden_dir), 'summarized_epitope_burden_updated.tsv')
	raw_summary = os.path.join(os.path.abspath(args.burden_dir), 'raw_epitope_table_updated.tsv')
	with open(combined_file, 'w') as o:
		header = ['Patient', 'Tumor', 'Epitope_by_mismatch_burden', 'Epitope_by_anchor_burden', 'Epitope_by_allele_burden', 
				  'Epitope_by_expression_burden', 'Epitope_by_TCGA_burden', 'Epitope_by_mismatch_and_allele_burden', 
				  'Epitope_by_mismatch_and_expression_burden', 'Epitope_by_mismatch_and_TCGA_burden', 
				  'Epitope_by_anchor_and_allele_burden', 'Epitope_by_anchor_and_expression_burden', 
				  'Epitope_by_anchor_and_TCGA_burden','Epitope_by_allele_and_expression_burden', 'Epitope_by_allele_and_TCGA_burden', 
				  'Epitope_by_mismatches_alleles_and_expression_burden', 'Epitope_by_mismatches_alleles_and_TCGA_burden', 
				  'Epitope_by_anchor_alleles_and_expression_burden', 'Epitope_by_anchor_alleles_and_TCGA_burden', 
				  'Alleles_presenting_epitopes', 'Alleles_presenting_per_mutation']
		print('\t'.join(header), file=o)
		with open(raw_summary, 'w') as o2:
			header = ['Patient', 'Tumor_ID', 'Epitope', 'Mutation', 'Mutation_type', 'Allele_count', 'Alleles',  'Disease', 'Adj_disease']
			print('\t'.join(header), file=o2)
			# Iterate through all patients to store samples
			patients = defaultdict(set)
			with open(os.path.abspath(args.manifest)) as p:
				for line in p:
					tokens = line.strip().split('\t')
					if tokens[2] != 'NA':
						# Get patient IDs
						patient_id = tokens[0]
						tumor = tokens[2]
						rna = tokens[3]
						if tumor == 'PGDX6340T_Ex_GCCAAT':
							patients[patient_id].add((tumor, rna, 'thyroid', 'MMR-deficient'))
						elif tumor in ['PGDX4311T_Ex', 'PGDX3787T_Ex']:
							patients[patient_id].add((tumor, rna, 'endometrial', 'MMR-deficient'))
						elif tumor.startswith('PGDX'):
							patients[patient_id].add((tumor, rna, 'colon', 'MMR-deficient'))
						elif tumor.startswith('Pem'):
							patients[patient_id].add((tumor, rna, 'prostate', 'Prostate'))
						elif tumor.startswith('SRR264'):
							patients[patient_id].add((tumor, rna, 'NSCLC', 'NSCLC'))
						elif tumor.startswith('SRR6') and not patient_id.startswith('DM0'):
							patients[patient_id].add((tumor, rna, 'RCC', 'RCC'))
						else:
							patients[patient_id].add((tumor, rna, 'melanoma', 'Melanoma'))
			# Iterate through all patients for study
			for patient_id in patients:
				tumor_data = []
				tumor_bams = []
				for tumor in patients[patient_id]:
					burden_file = os.path.join(os.path.abspath(args.burden_dir), ''.join([patient_id, '.', tumor[0], '.extended_epitope_burden.tsv']))
					tumor_data.append(process_extended_burden_file(burden_file))
					# Check if patient has expression data
					rna_bam = os.path.join(os.path.abspath(args.rna_dir), tumor[1], ''.join([tumor[1], '.sorted.bam']))
					if os.path.isfile(rna_bam):
						tumor_bams.append(rna_bam)
					with open(burden_file) as b:
						b.readline()
						for line in b:
							tokens = line.strip().split('\t')
							mutation_dict[tumor[3]][patient_id][tokens[3]].add('testtesttest')
							for hla in tokens[13].split(','):
								if hla != 'NA':
									mutation_dict[tumor[3]][patient_id][tokens[3]].add(hla)
									neoepitope_dict[tumor[3]][patient_id][tokens[2]].add(hla)
							out_line2 = tokens[0:5] + tokens[12:14] + [tumor[2], tumor[3]]
							print('\t'.join(out_line2), file=o2)
				out_line = [patient_id, ';'.join(sorted([tum[0] for tum in patients[patient_id]]))]
				for i in range(len(tumor_data[0])):
					# Check whether this is an patient-specific expression entry
					if i in [3, 6, 9, 11, 13, 15] and tumor_bams == []:
						# Patient has no expression data
						med_val = 'NA'
					else:
						try:
							med_val = str(median([x[i] for x in tumor_data if type(x[i]) != str]))
						except TypeError:
							med_val = 'NA'
					out_line.append(med_val)
				print('\t'.join(out_line), file=o)
