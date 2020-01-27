#!/usr/bin/env python

from __future__ import print_function
from collections import defaultdict
from numpy import median
import argparse
import copy
import glob
import os
import pickle

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", "-m",
                        type=str,
                        required=True,
                        help="path to manifest file with tumor-normal pair info")
    parser.add_argument("--output-dir", "-o",
                        type=str,
                        required=True,
                        help="path to output directory")
    parser.add_argument("--vcf-dir", "-v",
                        type=str,
                        required=True,
                        help="path to directory with consensus vcfs")
    parser.add_argument("--epitope-dir", "-e",
                        type=str,
                        required=True,
                        help="path to directory with neoepitope predictions")
    parser.add_argument("--binding-dir", "-b",
                        type=str,
                        required=True,
                        help="path to directory with binding prediction dictionaries")
    args = parser.parse_args()

	out_dir = os.path.abspath(args.output_dir)

	tumor_dict = {}
	mut_binders = defaultdict(lambda: defaultdict(set))
	driver_mut_binders = defaultdict(lambda: defaultdict(set))
	ep_binders = defaultdict(lambda: defaultdict(set))
	driver_ep_binders = defaultdict(lambda: defaultdict(set))

	mutation_lists = {}
	driver_mutation_lists = {}
	epitope_lists = defaultdict(set)
	driver_epitope_lists = defaultdict(set)

	with open(os.path.abspath(args.manifest)) as p:
		for line in p:
			# Extract sample identifiers
			tokens = line.strip().split('\t')
			patient = tokens[0]
			tumor = tokens[2]
			if tumor.startswith('PGDX'):
				disease = 'MMR-deficient'
			elif tumor.startswith('Pem'):
				disease = 'prostate'
			elif tumor.startswith('SRR264'):
				disease = 'NSCLC'
			elif tumor.startswith('SRR6') and not patient_id.startswith('DM0'):
				disease = 'RCC'
			else:
				'melanoma'
			# Set file paths
			vcf = os.path.join(os.path.abspath(args.vcf_dir), '.'.join([patient, tumor, 'consensus.vcf']))
			hg19_vars = ''.join([vcf, '.crv.hg19.var'])
			clinvar = ''.join([vcf, '.clinvar.var'])
			neoepitopes = os.path.join(os.path.abspath(args.epitope_dir), '.'.join([patient, tokens[2], 'neoepiscope.comprehensive.out']))
			# Find driver variant IDs
			driver_uids = set()
			with open(clinvar) as c:
				for line in c:
					if not line.startswith('#'):
						ctokens = line.strip().split('\t')
						if len(ctokens) > 1:
							if "Pathogenic" in ctokens[1] or "Likely pathogenic" in ctokens[1]:
								driver_uids.add(ctokens[0])
			# Get variant locations
			driver_variant_positions = set()
			with open(hg19_vars) as v:
				for line in v:
					if not line.startswith('#'):
						vtokens = line.strip().split('\t')
						if vtokens[0] in driver_uids:
							driver_variant_positions.add((vtokens[1], vtokens[2]))
			# Get full variant identities
			driver_variants = set()
			all_variants = set()
			with open(vcf) as v:
				for line in v:
					if not line.startswith('#'):
						vtokens = line.strip().split('\t')
						contig = ''.join(['chr', vtokens[0]])
						if vtokens[4] == "<DEL>" or vtokens[4] == '*':
							mutation_type = "D"
							pos = vtokens[1]
							ref = vtokens[3]
							alt = "*"
							all_variants.add((contig, pos, ref, alt, mutation_type))
							if (contig, vtokens[1]) in driver_variant_positions:
								driver_variants.add((contig, pos, ref, alt, mutation_type))
						elif len(vtokens[3]) == len(vtokens[4]):
							mutation_type = "V"
							pos = vtokens[1]
							ref = vtokens[3]
							alt = vtokens[4]
							all_variants.add((contig, pos, ref, alt, mutation_type))
							if (contig, vtokens[1]) in driver_variant_positions:
								driver_variants.add((contig, pos, ref, alt, mutation_type))
						elif len(vtokens[3]) > len(vtokens[4]):
							if vtokens[3].startswith(vtokens[4]):
								# Simple deletion
								mutation_type = "D"
								deletion_size = len(vtokens[3]) - len(vtokens[4])
								pos = str(int(vtokens[1]) + (len(vtokens[3]) - deletion_size))
								ref = vtokens[3][len(vtokens[4]):]
								alt = "*"
								all_variants.add((contig, pos, ref, alt, mutation_type))
								if (contig, vtokens[1]) in driver_variant_positions:
									driver_variants.add((contig, pos, ref, alt, mutation_type))
							else:
								# Complex indel
								# Add deletion first
								mutation_type = "D"
								pos = vtokens[1]
								ref = vtokens[3]
								alt = "*"
								all_variants.add((contig, pos, ref, alt, mutation_type))
								if (contig, vtokens[1]) in driver_variant_positions:
									driver_variants.add((contig, pos, ref, alt, mutation_type))
								# Then add insertion
								mutation_type = "I"
								pos = str(int(vtokens[1]) + len(vtokens[3]) - 1)
								ref = "*"
								alt = vtokens[4]
								all_variants.add((contig, pos, ref, alt, mutation_type))
								if (contig, vtokens[1]) in driver_variant_positions:
									driver_variants.add((contig, pos, ref, alt, mutation_type))
						elif len(vtokens[3]) < len(vtokens[4]):
							if vtokens[4].startswith(vtokens[3]):
								# Simple insertion
								mutation_type = "I"
								pos = str(int(vtokens[1]) + len(vtokens[3]) - 1)
								ref = "*"
								alt = vtokens[4][len(vtokens[3]) :]
								all_variants.add((contig, pos, ref, alt, mutation_type))
								if (contig, vtokens[1]) in driver_variant_positions:
									driver_variants.add((contig, pos, ref, alt, mutation_type))
							else:
								# Complex indel - add deletion first
								mutation_type = "D"
								pos = vtokens[1]
								ref = vtokens[3]
								alt = "*"
								all_variants.add((contig, pos, ref, alt, mutation_type))
								if (contig, vtokens[1]) in driver_variant_positions:
									driver_variants.add((contig, pos, ref, alt, mutation_type))
								# Then add insertion
								mutation_type = "I"
								pos = str(int(vtokens[1]) + len(vtokens[3]) - 1)
								ref = "*"
								alt = vtokens[4]
								all_variants.add((contig, pos, ref, alt, mutation_type))
								if (contig, vtokens[1]) in driver_variant_positions:
									driver_variants.add((contig, pos, ref, alt, mutation_type))
			# Get driver epitope sequences
			driver_epitopes = defaultdict(set)
			all_epitopes = defaultdict(set)
			driver_neopeptide_count = 0
			neopeptide_count = 0
			with open(neoepitopes) as n:
				n.readline()
				n.readline()
				for line in n:
					ntokens = line.strip().split('\t')
					if tuple(ntokens[1:6]) in driver_variants:
						driver_epitopes[tuple(ntokens[1:6])].add(ntokens[0])
						driver_neopeptide_count += 1
					neopeptide_count += 1
					all_epitopes[tuple(ntokens[1:6])].add(ntokens[0])
			# Store all mutations/epitopes
			if patient not in mutation_lists:
				mutation_lists[(patient, disease)] = copy.copy(all_variants)
				driver_mutation_lists[(patient, disease)] = copy.copy(driver_variants)
			else:
				mutation_lists[(patient, disease)].update(all_variants)
				mutation_lists[(patient, disease)].update(driver_variants)
			for variant in all_epitopes:
				epitope_lists[(patient, disease)].update(all_epitopes[variant])
			for variant in driver_epitopes:
				driver_epitope_lists[(patient, disease)].update(driver_epitopes[variant])
			# Get binding scores to filter epitopes
			filtered_epitopes = defaultdict(set)
			filtered_epitopes_count = 0
			driver_filtered_epitopes = defaultdict(set)
			driver_filtered_epitopes_count = 0
			binding_dicts = glob.glob(os.path.join(os.path.abspath(args.binding_dir), '.'.join([patient, tokens[2], 'mhcnuggets.*.pickle'])))
			scores = defaultdict(dict)
			for dic in binding_dicts:
				allele = dic.split('.')[-2]
				with open(dic, 'rb') as d:
					b_scores = pickle.load(d)
				for peptide in b_scores:
					if b_scores[peptide] != 'NA':
						scores[peptide][allele] = float(b_scores[peptide])
			for variant in all_epitopes:
				for ep in list(all_epitopes[variant]):
					binding_alleles = [x for x in scores[ep] if scores[ep][x] < 500.0]
					if binding_alleles:
						filtered_epitopes[variant].add(ep)
						filtered_epitopes_count += 1
						if ep in driver_epitopes[variant]:
							driver_filtered_epitopes[variant].add(ep)
							driver_filtered_epitopes_count += 1
					for allele in binding_alleles:
						mut_binders[(patient, disease)][variant].add(allele)
						ep_binders[(patient, disease)][ep].add(allele)
						if ep in driver_epitopes[variant]:
							driver_mut_binders[(patient, disease)][variant].add(allele)
							driver_ep_binders[(patient, disease)][ep].add(allele)
			tumor_dict[(patient, tokens[2], disease)] = [len(all_variants), neopeptide_count, len(filtered_epitopes), filtered_epitopes_count,
														 len(driver_variants), driver_neopeptide_count, len(driver_filtered_epitopes), driver_filtered_epitopes_count]


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
		# Set up new combined key/entrie
		new_key = (patient, ';'.join(sorted([x[1] for x in relevant_keys])), relevant_keys[0][2])
		new_entry = []
		for i in range(len(relevant_entries[0])):
			new_entry.append(median([float(x[i]) for x in relevant_entries]))
		assert len(new_entry) == len(relevant_entries[0])
		tumor_dict[new_key] = new_entry
		for key in relevant_keys:
			del tumor_dict[key]

	with open(os.path.join(out_dir, 'drivers.tsv'), 'w') as o:
		header = ['Patient', 'Tumor', 'Disease', 'Total_variants', 'Total_neopeptides', 'Presented_variants', 'Presented_epitopes',
					 'Total_clinvar_variants', 'Total_clinvar_neopeptides', 'Presented_clinvar_variants', 'Presented_clinvar_epitopes']
		print('\t'.join(header), file=o)
		for tumor in tumor_dict:
			out_line = list(tumor) + tumor_dict[tumor]
			print('\t'.join([str(x) for x in out_line]), file=o)

	with open(os.path.join(out_dir, 'mut_binders.tsv'), 'w') as o:
		header = ['Patient', 'Disease', 'Mutation', 'Allele_count']
		print('\t'.join(header), file=o)
		for patient in mutation_lists:
			for mutation in mutation_lists[patient]:
				out_line = [patient[0], patient[1], ','.join(list(mutation)), str(len(mut_binders[patient][mutation]))]
				print('\t'.join(out_line), file=o)

	with open(os.path.join(out_dir, 'driver_mut_binders.tsv', 'w')) as o:
		header = ['Patient', 'Disease', 'Mutation', 'Allele_count']
		print('\t'.join(header), file=o)
		for patient in driver_mutation_lists:
			for mutation in driver_mutation_lists[patient]:
				out_line = [patient[0], patient[1], ','.join(list(mutation)), str(len(driver_mut_binders[patient][mutation]))]
				print('\t'.join(out_line), file=o)

	with open(os.path.join(out_dir, 'ep_binders.tsv'), 'w') as o:
		header = ['Patient', 'Disease', 'Epitope', 'Allele_count']
		print('\t'.join(header), file=o)
		for patient in epitope_lists:
			for ep in epitope_lists[patient]:
				out_line = [patient[0], patient[1], ep, str(len(ep_binders[patient][ep]))]
				print('\t'.join(out_line), file=o)

	with open(os.path.join(out_dir, 'driver_ep_binders.tsv'), 'w') as o:
		header = ['Patient', 'Disease', 'Epitope', 'Allele_count']
		print('\t'.join(header), file=o)
		for patient in driver_epitope_lists:
			for ep in driver_epitope_lists[patient]:
				out_line = [patient[0], patient[1], ep, str(len(driver_ep_binders[patient][ep]))]
				print('\t'.join(out_line), file=o)

