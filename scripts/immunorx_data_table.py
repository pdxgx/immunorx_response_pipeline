#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
from collections import defaultdict
from neoepiscope import bowtie_index, paths
from neoepiscope.transcript import Transcript
from numpy import median
import argparse
import glob
import os
import pickle
import copy
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', '-o', type=str,required=True, help='path to directory to write output file')
    parser.add_argument('--coverage-file', '-c', type=str, required=True, help='path to file with coverage data')
    parser.add_argument('--sample-data', '-s', type=str, required=True, help='path to file with sample data')
    parser.add_argument("--manifest", "-m", type=str, required=True, help="path to tumor-normal pair manifest file")
   	parser.add_argument(
    						"--consensus-vcf-dir", "-v", type=str, required=True, 
    						help="path to directory containing consensus VCFs for each patient"
    )
    parser.add_argument(
    						"--raw-vcf-dir", "-r", type=str, required=True, 
    						help="path to directory containing raw vcf subdirectories for each patient"
    )
    parser.add_argument(
    						"--hla-type-dir", "-t", type=str, required=True, 
    						help="path to directory containing HLA typing information for each patient"
    )
    parser.add_argument(
    						"--msi-dir", "-i", type=str, required=True, 
    						help="path to directory containing MSI information for each patient"
    )
    parser.add_argument(
    						"--neoepitope-dir", "-n", type=str, required=True, 
    						help="path to directory containing neoepitope predictions"
    )
    parser.add_argument(
    						"--binding-dir", "-n", type=str, required=True, 
    						help="path to directory containing binding affinity prediction dictionaries"
    )
    parser.add_argument(
    						"--unphased-binding-dir", "-u", type=str, required=True, 
    						help="path to directory containing binding affinity prediction dictionaries"
    )
    parser.add_argument(
    						"--hla-reference-fasta", "-f", type=str, required=True, 
    						help="path to HLA reference FASTA file from Optitype"
    )
    args = parser.parse_args()

	# Load transcript dictionaries and bowtie index from neoepiscope
	with open(os.path.join(paths.gencode_v19,'intervals_to_transcript.pickle'), 'rb') as interval_stream:
		interval_dict = pickle.load(interval_stream)
	with open(os.path.join(paths.gencode_v19, 'transcript_to_CDS.pickle'), 'rb') as cds_stream:
		cds_dict = pickle.load(cds_stream)
	reference_index = bowtie_index.BowtieIndexReference(paths.bowtie_hg19)

	# Get HLA decoder
	hla_id_linker = {}
	hla_fasta = os.path.abspath(args.hla_reference_fasta)
	for record in SeqIO.parse(hla_fasta, 'fasta'):
		header = str(record.description).split()
		hla_id_linker[header[0]] = ':'.join(header[1].split(':')[0:2])

	# Dictionary to account for changes in naming conventions
	hopkins2_dict = {
						'Hopkins_6_101' : 'PGDX6340',
						'Hopkins_4_002_recurrent' : 'PGDX6334',
						'Hopkins_4_002' : 'PGDX4947',
						'Hopkins_44' : 'PGDX4398',
						'Hopkins_52' : 'PGDX4397',
						'Hopkins_29' : 'PGDX4311',
						'Hopkins_25' : {
											'Bowel (small)' : 'PGDX4137T_Ex', 
											'Brain' : 'PGDX4137T1_Ex', 
											'Blood': 'PGDX4137N_Ex'
						},
						'Hopkins_42' : 'PGDX4127',
						'Hopkins_36' : 'PGDX3799',
						'Hopkins_30' : 'PGDX3787',
						'Hopkins_20' : 'PGDX3785',
						'Hopkins_33' : 'PGDX3784',
						'Hopkins_31' : 'PGDX3783'
	}

	# Produce coverage dictionary - Mbp w/ min 6 reads coverage
	coverage_dict = {}
	with open(os.path.abspath(args.coverage_file)) as c:
		c.readline()
		for line in c:
			tokens = line.strip('\n').split('\t')
			coverage_dict[(tokens[0], tokens[1])] = tokens[2]

	# Produce dictionary with data from immunotherapy sample key
	key_dict = defaultdict(dict)
	with open(os.path.abspath(args.sample_data)) as k:
		k.readline()
		for line in k:
			tokens = line.strip('\n').split('\t')
			# Determine patient/tumor IDs
			if 'Hopkins_' in tokens[0] and 'PGD' in tokens[13]:
				sample = tokens[13].replace('_R', '')
				patient = sample.replace('_Ex', '')[:-1]
			elif 'Hopkins_' in tokens[0]:
				if tokens[0] == 'Hopkins_25':
					patient = 'PGDX4137'
					sample = hopkins2_dict['Hopkins_25'][tokens[3]]
				else:
					patient = hopkins2_dict[tokens[0]]
					if patient == 'PGDX4947':
						if tokens[2] == 'T':
							sample = ''.join([patient, tokens[2], '_Ex_ATGTCA'])
						else:
							sample = ''.join([patient, tokens[2], '_Ex_AGTTCC'])
					elif patient == 'PGDX6334':
						if tokens[2] == 'T':
							sample = ''.join([patient, tokens[16], '_Ex_TGACCA'])
						else:
							sample = ''.join([patient, tokens[2], '_Ex'])
					elif patient == 'PGDX6340':
						if tokens[2] == 'T':
							sample = ''.join([patient, tokens[2], '_Ex_GCCAAT'])
						else:
							sample = ''.join([patient, tokens[2], '_Ex_ACAGTG'])
					else:
						sample = ''.join([patient, tokens[2], '_Ex'])
			elif '11025_' in tokens[0]:
				patient = tokens[0].replace('_', '-')
				if patient[-2:] == '01':
					patient = '11025-1'
				patient_number = tokens[0][-2:]
				if tokens[2] == 'T':
					tn = 'T'
				else:
					tn = 'N'
				if patient_number == '01':
					patient_number = '1'
				sample = ''.join(['Pem', patient_number, tn])
			else:
				patient = tokens[0]
				sample = tokens[13]
			# Cancer stage
			if 'I' in tokens[1] or 'M' in tokens[1]:
				stage = tokens[1]
			elif tokens[12] == 'PRJNA420786':
				stage = 'IV'
			else:
				stage = 'NA'
			# Treatment
			if tokens[4] == 'NA' or tokens[4] == '':
				treatment = 'NA'
				pd1 = '0'
				actla4 = '0'
				pdl1 = '0'
				other = '0'
			else:
				if 'Pembrolizumab' in tokens[4]:
					pd1 = '1'
					if 'Ipilimumab' in tokens[4]:
						actla4 = '1'
					else:
						actla4 = '0'
					pdl1 = '0'
					other = '0'
				elif 'Nivolumab' in tokens[4]:
					pd1 = '1'
					pdl1 = '0'
					other = '0'
					if 'Ipilimumab' in tokens[4]:
						actla4 = '1'
					else:
						actla4 = '0'
				elif 'nivolumab' in tokens[4]:
					pd1 = '1'
					actla4 = '0'
					pdl1 = '0'
					other = '0'
				elif 'Ipilimumab' in tokens[4]:
					tokens[4] = tokens[4].strip('"')
					pd1 = '0'
					actla4 = '1'
					pdl1 = '0'
					if '+dacarbazine' in tokens[4]:
						other = '1'
					elif 'dacarbazine' in tokens[4]:
						other = '1'
					elif '+vemurafenib' in tokens[24]:
						other = '1'
					elif 'vemurafenib' in tokens[4]:
						other = '1'
					else:
						other = '0'
				elif 'Tremelimumab' in tokens[4]:
					pd1 = '0'
					actla4 = '1'
					pdl1 = '0'
					other = '0'
				elif 'Atezolizumab' in tokens[4]:
					pd1 = '0'
					actla4 = '0'
					pdl1 = '1'
					other = '0'
				elif tokens[4] == 'aCTLA-4 + aPD-1':
					pd1 = '1'
					actla4 = '1'
					pdl1 = '0'
					other = '0'
				elif 'aCTLA-4' in tokens[4]:
					pd1 = '0'
					actla4 = '1'
					pdl1 = '0'
					other = '0'
				elif 'aPD-1' in tokens[4]:
					pd1 = '1'
					actla4 = '0'
					pdl1 = '0'
					other = '0'
				elif 'Control Ab' in tokens[4]:
					pd1 = 'NA'
					actla4 = 'NA'
					pdl1 = 'NA'
					other = 'NA'
				elif 'Interferon-gamma' in tokens[4]:
					pd1 = '0'
					actla4 = '0'
					pdl1 = '0'
					other = '1'
			# CTLA4 response
			if tokens[5] == 'NA' or tokens[5] == '':
				ctla4_response = 'NA'
				ctla4 = '0'
			else:
				if 'Y' in tokens[5]:
					ctla4_response = '1'
					actla4 = '1'
				else:
					ctla4_response = '0'
					actla4 = '1'
			# PD1 response
			if tokens[6] == 'NA' or tokens[6] == '':
				pd1_response = 'NA'
				pd1 = '0'
			else:
				if tokens[6] == 'Y' or 'intermediate' in tokens[6]:
					pd1_response = '1'
					pd1 = '1'
				else:
					pd1_response = '0'
					pd1 = '1'
			# RECIST
			if tokens[17] == 'NA' or tokens[17] == '':
				recist = 'NA'
			else:
				if 'Stable Disease' in tokens[17] or 'SD' in tokens[17]:
					recist = 'SD'
				elif 'Partial Response' in tokens[17] or 'PR' in tokens[17]:
					recist = 'PR'
				elif 'Complete Response' in tokens[17] or 'CR' in tokens[17]:
					recist = 'CR'
				elif 'Progressive Disease' in tokens[17] or 'PD' in tokens[17]:
					recist = 'PD'
				else:
					recist = 'X'
			# Combined response
			if ctla4_response == '1' or pd1_response == '1':
				combined_response = '1'
			elif recist == 'PR' or recist == 'CR':
				combined_response = '1'
			elif ctla4_response == '0' or pd1_response == '0':
				combined_response = '0'
			elif recist == 'SD' or recist == 'PD':
				combined_response = '0'
			else:
				combined_response = 'NA'
			# Response duration
			if tokens[7] == 'NA' or tokens[7] == '':
				response_duration = 'NA'
			else:
				response_duration = tokens[7]
			# PFS
			if tokens[8] == 'NA' or tokens[8] == '':
				if response_duration != 'NA':
					pfs = copy.copy(response_duration)
				else:
					pfs = 'NA'
			else:
				pfs = tokens[8]
			# OS
			if tokens[9] == 'NA' or tokens[9] == '':
				OS = 'NA'
			else:
				OS = tokens[9]
			# Cesnsoring
			if tokens[10] == 'NA' or tokens[10] == '':
				censoring = 'NA'
			else:
				censoring = tokens[10]
			# Check PFS event status		
			if pfs == 'NA':
				pfs_event = 'NA'
			elif pfs == OS:
				pfs_event = '0'
			else:
				if OS != 'NA':
					pfs_event = '1'
				elif censoring != 'NA':
					pfs_event = '0'
				else:
					pfs_event = '1'
			# Check OS event status, vital status, censoring
			if OS == 'NA':
				os_event = 'NA'
				if 'Alive' in tokens[11]:
					vital = '0'
					if pfs != 'NA':
						censor_stat = '1'
					else:
						censor_stat = 'NA'
				elif 'Dead' in tokens[11]:
					vital = '1'
					if pfs != 'NA':
						censor_stat = '0'
					else:
						censor_stat = 'NA'
				else:
					if censoring != 'NA':
						vital = '0'
						censor_stat = '1'
					elif pfs != 'NA':
						vital = 'NA'
						censor_stat = '0'
					else:
						vital = 'NA'
						censor_stat = 'NA'
			else:
				if censoring != 'NA' or 'Alive' in tokens[11]:
					os_event = '0'
					vital = '0'
					censor_stat = '1'
				elif 'Dead' in tokens[11]:
					vital = '1'
					censor_stat = '0'
					if 'unrelated' in tokens[11]:
						os_event = '0'
					else:
						os_event = '1'
				else:
					os_event = '1'
					vital = '1'
					censor_stat = '0'
			# Get author-reported burdens
			non_synon = tokens[14]
			all_muts = tokens[15]
			neoantigens = tokens[16]
			# Get disease/study of origin info
			if tokens[12] not in ['PRJNA420786', 'PRJNA293912', 'Spellman11025', 'Hopkins', 'Hopkins2']:
				disease = 'melanoma'
				if tokens[12] in ['PRJNA312948', 'PRJNA307199', 'PRJNA343789']:
					study = 'hugo'
				elif tokens[12] in ['PRJNA306070', 'PRJNA305077']:
					study = 'snyder'
				elif tokens[12] == ['PRJNA278450']:
					study = 'carreno'
				elif tokens[12] == ['PRJNA82745']:
					study = 'vanallen'
				elif tokens[12] == ['PRJNA324705']:
					study = 'zaretsky'
				elif tokens[12] == ['PRJNA357321']:
					study = 'gao'
				elif tokens[12] == ['PRJNA369259']:
					study = 'roh'
				elif tokens[12] == 'EGAD00001004352':
					study = 'amaria'
				elif tokens[12] == 'PRJNA414014':
					study = 'eroglu'
			elif tokens[12] == 'PRJNA420786':
				study = 'miao'
				disease = 'RCC'
			elif tokens[12] == 'PRJNA293912':
				study = 'rizvi'
				disease = 'NSCLC'
			elif tokens[12] == 'Spellman11025':
				study = 'graff'
				disease = 'prostate'
			elif tokens[12] == 'Hopkins':
				study = 'le'
				disease = 'colon'
			elif tokens[12] == 'Hopkins2':
				study = 'le'
				if tokens[13] not in ['PGDX6340', 'PGDX4311', 'PGDX3787']:
					disease = 'colon'
				elif tokens[13] in ['PGDX4311', 'PGDX3787']:
					disease = 'endometrial'
				else:
					disease = 'thyroid'
			key_dict[patient][sample] = [
											stage, pd1, pdl1, actla4, other, ctla4_response, pd1_response, 
											combined_response, pfs, OS, censoring, vital, os_event, pfs_event, 
											censor_stat, non_synon, all_muts, neoantigens, study, disease
										]

	# Gather genomic data and link to metadata
	tumor_dict = {}
	with open(os.path.abspath(args.manifest)) as f:
		for line in f:
			tokens = line.strip('\n').split('\t')
			if tokens[2] != 'NA':
				# Get study ID
				study = key_dict[tokens[0]][tokens[2]][-2]
				# Gather coverage data
				coverage = coverage_dict[(tokens[0], tokens[2])]
				# Determine VCF path
				tumor = copy.copy(tokens[2])
				vcf = os.path.join(
								os.path.abspath(args.consensus_vcf_dir), '.'.join([tokens[0], tumor, 'consensus.vcf'])
				)
				# Set up counters for different variant types
				total_muts, snvs, inframe_insertions, inframe_deletions, fs_insertions, fs_deletions = 0, 0, 0, 0, 0, 0
				nonsynonymous = 0
				muse, muse_snv, muse_del, muse_ins = 0, 0, 0, 0,
				mutect, mutect_snv, mutect_del, mutect_ins =  0, 0, 0, 0
				pindel, pindel_snv, pindel_del, pindel_ins = 0, 0, 0, 0
				radia, radia_snv, radia_del, radia_ins = 0, 0, 0, 0, 
				somatic_sniper, somatic_sniper_snv = 0, 0
				somatic_sniper_del, somatic_sniper_ins = 0, 0
				varscan, varscan_snv, varscan_del, varscan_ins = 0, 0, 0, 0
				# Iterate through VCF for variant data
				with open(vcf) as v:
					for line in v:
						if line[0] != '#':
							total_muts += 1
							columns = line.strip('\n').split('\t')
							# SNVs
							if len(columns[3]) == len(columns[4]):
								snvs += 1
								# Determine mutation consequences
								contig = ''.join(['chr', columns[0]])
								if contig in interval_dict:
									if len(interval_dict[contig][int(columns[1])]) > 0:
										peptide_dict = defaultdict(list)
										for interval in interval_dict[contig][int(columns[1])]:
											transcript = interval.data
											tx = Transcript(
												reference_index, 
												[[str(chrom), '.', seq_type, str(start), str(end), '.', strand] for (
																		chrom, seq_type, start, end, strand, tx_type
															) in cds_dict[transcript]], 
	                      						transcript
		                    				)
											tx.edit(
														columns[4], int(columns[1]), mutation_type='V',
														 mutation_class='S', vaf=None
											)
											peps = tx.neopeptides()
											for pep in peps:
												peptide_dict[pep].extend(peps[pep])
										if peptide_dict:
											nonsynonymous += 1
							# Deletion
							elif len(columns[3]) > len(columns[4]):
								if ((len(columns[3]) - len(columns[4])) % 3):
									fs_deletions += 1
								else:
									inframe_deletions += 1
							# Insertion
							elif len(columns[3]) < len(columns[4]):
								if ((len(columns[4]) - len(columns[3])) % 3):
									fs_insertions += 1
								else:
									inframe_insertions += 1
				# Get VCF paths from individual callers
				muse_vcf = os.path.join(os.path.abspath(args.raw_vcf_dir), tokens[2], 'muse.reheadered.vcf')
				mutect_vcf = os.path.join(os.path.abspath(args.raw_vcf_dir), tokens[2], 'mutect.reheadered.vcf')
				pindel_vcf = os.path.join(os.path.abspath(args.raw_vcf_dir), tokens[2], 'pindel.reheadered.vcf')
				radia_vcf = os.path.join(os.path.abspath(args.raw_vcf_dir), tokens[2], 'radia_filtered.reheadered.vcf')
				somatic_sniper_vcf = os.path.join(
													os.path.abspath(args.raw_vcf_dir), tokens[2], 
													'somatic_sniper_fpfilter.reheadered.vcf'
				)
				varscan_vcf = os.path.join(
											os.path.abspath(args.raw_vcf_dir), tokens[2], 
											'varscan_fpfilter.reheadered.vcf'
				)
				# Process MuSE VCF
				with open(muse_vcf, 'r') as v:
					for line in v:
						if line[0] != '#':
							columns = line.strip('\n').split('\t')
							if columns[6] != 'Tier5':
								muse += 1
								if len(columns[3]) == len(columns[4]):
									muse_snv += 1
								elif len(columns[3]) > len(columns[4]):
									muse_del += 1
								elif len(columns[3]) < len(columns[4]):
									muse_ins += 1
				# Process MuTect VCF
				with open(mutect_vcf, 'r') as v:
					for line in v:
						if line[0] != '#':
							columns = line.strip('\n').split('\t')
							if columns[6] != 'REJECT':
								mutect += 1
								if len(columns[3]) == len(columns[4]):
									mutect_snv += 1
								elif len(columns[3]) > len(columns[4]):
									mutect_del += 1
								elif len(columns[3]) < len(columns[4]):
									mutect_ins += 1
				# Process Pindel VCF
				with open(pindel_vcf, 'r') as v:
					for line in v:
						if line[0] != '#':
							pindel += 1
							columns = line.strip('\n').split('\t')
							if len(columns[3]) == len(columns[4]):
								pindel_snv += 1
							elif len(columns[3]) > len(columns[4]):
								pindel_del += 1
							elif len(columns[3]) < len(columns[4]):
								pindel_ins += 1
				# Process RADIA VCF
				with open(radia_vcf, 'r') as v:
					for line in v:
						if line[0] != '#':
							radia += 1
							columns = line.strip('\n').split('\t')
							if len(columns[3]) == len(columns[4]):
								radia_snv += 1
							elif len(columns[3]) > len(columns[4]):
								radia_del += 1
							elif len(columns[3]) < len(columns[4]):
								radia_ins += 1
				# Process SomaticSniper VCF
				with open(somatic_sniper_vcf, 'r') as v:
					for line in v:
						if line[0] != '#':
							columns = line.strip('\n').split('\t')
							if columns[6] == 'PASS':
								somatic_sniper += 1
								if len(columns[3]) == len(columns[4]):
									somatic_sniper_snv += 1
								elif len(columns[3]) > len(columns[4]):
									somatic_sniper_del += 1
								elif len(columns[3]) < len(columns[4]):
									somatic_sniper_ins += 1
				# Process VarScan VCF
				with open(varscan_vcf, 'r') as v:
					for line in v:
						if line[0] != '#':
							columns = line.strip('\n').split('\t')
							if columns[6] == 'PASS':
								varscan += 1
								if len(columns[3]) == len(columns[4]):
									varscan_snv += 1
								elif len(columns[3]) > len(columns[4]):
									varscan_del += 1
								elif len(columns[3]) < len(columns[4]):
									varscan_ins += 1
				# MHC I allele data
				optitype_tumor = os.path.join(os.path.abspath(args.hla_type_dir), '_'.join([tokens[2], 'result.tsv']))
				with open(optitype_tumor, 'r') as o:
					o.readline()
					columns = o.readline().strip('\n').split('\t')
					tumor_hla1 = set()
					for i in range(1, 7):
						if columns[i] != '':
							if 'HLA' in columns[i]:
								hla = hla_id_linker[columns[i]]
							else:
								hla = '-'.join(['HLA', columns[i]])
							tumor_hla1.add(hla)
				tumor_hla1_count = len(tumor_hla1)
				# MHC II allele data
				tumor_hla2 = set()
				seq2hla_tumor = os.path.join(
												os.path.abspath(args.hla_type_dir), 
												'-'.join([tumor, 'ClassII.HLAgenotype4digits'])
				)
				with open(seq2hla_tumor, 'r') as o:
					o.readline()
					for line in o:
						columns = line.strip('\n').split('\t')
						for i in [1,3]:
							allele_opt = columns[i].split(',')
							for allele in allele_opt:
								if allele != 'no':
									tumor_hla2.add(''.join(['HLA-', allele.strip("'")]))
				tumor_hla2_count = len(tumor_hla2)
				patient_alleles = [x for x in tumor_hla_1] + [x for x in tumor_hla2]
 				# Get MSI data
				msi_file = glob.glob(
								os.path.join(
												os.path.abspath(args.msi_dir), ''.join([tumor, '*']), 
												'.'.join([tumor, 'reheadered.realigned.cleaned.MSI_Analysis.txt'])
								)
				)[0]
				with open(msi_file) as msi_stream:
					for i in range(4):
						msi_stream.readline()
					stat_info = msi_stream.readline().strip().split('\t')
				if stat_info == 'POS':
					msi_status = '1'
				else:
					msi_status = '0'
				# Prep dictionary with neoepitope data
				unphased_neoepitopes = os.path.join(
														os.path.abspath(args.neoepitope_dir),
														'.'.join([tokens[0], tokens[2], 'neoepiscope.somatic.out'])
				)
				comprehensive_neoepitopes = os.path.join(
													os.path.abspath(args.neoepitope_dir),
													'.'.join([tokens[0],tokens[2],'neoepiscope.comprehensive.out'])
				)
				ep_dict = defaultdict(lambda: defaultdict(set))
				# Extract unphased epitope data
				with open(unphased_neoepitopes, 'r') as n:
					n.readline()
					n.readline()
					for line in n:
						columns = line.strip('\n').split('\t')
						data = tuple(columns[1:6])
						ep_dict[columns[0]]['unphased'].add(data)
				# Extract phased epitope data
				with open(comprehensive_neoepitopes, 'r') as n:
					n.readline()
					n.readline()
					for line in n:
						columns = line.strip('\n').split('\t')
						data = tuple(columns[1:6])
						ep_dict[columns[0]]['comprehensive'].add(data)
				# Load binding score data
				mhcnuggets = {}
				for allele in patient_alleles:
					mhcnuggets[allele] = {}
					for ending in ['pickle', 'unphased.pickle']:
						allele_dict = os.path.join(
													os.path.abspath(args.binding_dir),
													'.'.join([tokens[0], tokens[2], 'mhcnuggets', allele, '.pickle'])

						)
						with open(allele_dict, 'rb') as pic:
							d = pickle.load(pic)
						for peptide in d:
							mhcnuggets[allele][peptide] = d[peptide]
				netMHCpan = {}
				for allele in patient_alleles:
					netMHCpan[allele] = {}
					for ending in ['pickle', 'unphased.pickle']:
						allele_dict = os.path.join(
													os.path.abspath(args.unphased_binding_dir),
													'.'.join([tokens[0], tokens[2], 'netMHCpan', allele, '.pickle'])

						)
						with open(allele_dict, 'rb') as pic:
							d = pickle.load(pic)
						for peptide in d:
							netMHCpan[allele][peptide] = d[peptide]
				# Gather counts of different neoepitope types
				total_u_eps, total_c_eps = 0, 0
				mhcnuggets_binding_eps, mhcnuggets_class1_binding_eps, mhcnuggets_class2_binding_eps = 0, 0, 0
				binding_man_eps = 0
				# Process epitope data
				for ep in ep_dict:
					mut_list = []
					for k in ep_dict[ep].keys():
						mut_list.extend(list(ep_dict[ep][k]))
					mut_list = list(set(mut_list))
					# Process unphased epitopes
					if 'unphased' in ep_dict[ep]:
						total_u_eps += 1
						if study in ['rizvi', 'hugo', 'roh', 'vanallen', 'carreno']:
							# Patient is from a study that has epitope burdens reported
							proceed = False
							if len(ep) == 9:
								proceed = True
							elif len(ep) == 10 and study in ['hugo', 'roh', 'vanallen']:
								proceed = True
							elif len(ep) == 11 and study in ['hugo', 'roh']:
								proceed = True
							if proceed:
								# Epitope is of acceptable size for this study
								if study == 'carreno':
									# Must check that epitope binds to HLA-A*02:01
									if 'HLA-A*02:01' in netMHCpan:
										if netMHCpan['HLA-A*02:01'][ep] != 'NA':
											if float(netMHCpan['HLA-A*02:01'][ep]) < 500.0:
												binding_man_eps += 1
								else:
									# Must check that epitope binds to any MHC class I allele
									binds = False
									for allele in patient_alleles:
										if 'HLA-A' in allele or 'HLA-B' in allele or 'HLA-C' in allele:
											if allele in netMHCpan:
												if netMHCpan[allele][ep] != 'NA':
													if float(netMHCpan[allele][ep]) < 500.0:
														binds = True
									if binds == True:
										binding_man_eps += 1
					# Process phased epitopes
					if 'comprehensive' in ep_dict[ep]:
						total_c_eps += 1
						# Process mhcnuggets binding score data
						mhcnuggets_mhc, mhcnuggets_mhc1, mhcnuggets_mhc2 = set(), set(), set()
						for allele in mhcnuggets:
							if ('HLA-A' in allele) or ('HLA-B' in allele) or ('HLA-C' in allele):
								if mhcnuggets[allele][ep] != 'NA':
									mhcnuggets_mhc1.add(mhcnuggets[allele][ep])
									mhcnuggets_mhc.add(mhcnuggets[allele][ep])
							else:
								if mhcnuggets[allele][ep] != 'NA':
									mhcnuggets_mhc2.add(mhcnuggets[allele][ep])
									mhcnuggets_mhc.add(mhcnuggets[allele][ep])
						if [x for x in mhcnuggets_mhc1 if float(x) <= 500.0] != []:
							mhcnuggets_class1_binding_eps += 1
						if [x for x in mhcnuggets_mhc2 if float(x) <= 500.0] != []:
							mhcnuggets_class2_binding_eps += 1
						if [x for x in mhcnuggets_mhc if float(x) <= 500.0] != []:
							mhcnuggets_binding_eps += 1
				# Gather data from the sample key
				metadata = []
				tumor_data = key_dict[tokens[0]][tokens[2]]
				normal_data = key_dict[tokens[0]][tokens[1]]
				for i in range(0, len(tumor_data)):
					if tumor_data[i] == normal_data[i]:
						metadata.append(tumor_data[i])
					elif tumor_data[i] != 'NA' and normal_data[i] == 'NA':
						metadata.append(tumor_data[i])
					elif tumor_data[i] == 'NA' and normal_data[i] != 'NA':
						metadata.append(normal_data[i])
					else:
						metadata.append(tumor_data[i])
				# Store all tumor data in dictionary
				tumor_dict[(tokens[0], tokens[1], tokens[2])] = [
														   			coverage, total_muts, snvs, inframe_insertions, 
														   			inframe_deletions, fs_insertions, fs_deletions, 
														   			nonsynonymous, muse, muse_snv, muse_del, muse_ins, 
														   			mutect, mutect_snv, mutect_del, mutect_ins, pindel, 
														   			pindel_snv, pindel_del, pindel_ins, radia, 
														   			radia_snv, radia_del, radia_ins, somatic_sniper, 
														   			somatic_sniper_snv, somatic_sniper_del, 
														   			somatic_sniper_ins, varscan, varscan_snv, 
														   			varscan_del, varscan_ins, tumor_hla1_count,
														   			tumor_hla2_count, msi_status, total_u_eps, 
														   			total_c_eps, mhcnuggets_binding_eps, 
														   			mhcnuggets_class1_binding_eps, 
														   			mhcnuggets_class2_binding_eps, binding_man_eps
														   		]
				tumor_dict[(tokens[0], tokens[1], tokens[2])].extend(metadata)

	# Get median values for multi-sample patients
	multisample = set([x[0] for x in tumor_dict if len([y for y in tumor_dict if x[0] in y]) > 1])
	for patient in multisample:
		# Extract relevant keys/entries
		try:
			# Separate patients from Roh/Amaria cohorts (#s as pat. IDs)
			p = int(patient)
			relevant_keys = [x for x in tumor_dict if x[0] == patient and (
				  x[2][-1] in ['A', 'B', 'C', 'D', 'E'] or x[2] in [''.join([patient, 'D1']), ''.join([patient, 'D2'])]
				)
			]
			if len(relevant_keys) == 1:
				continue
		except ValueError:
			relevant_keys = [x for x in tumor_dict if x[0] == patient]
		relevant_entries = [tumor_dict[x] for x in relevant_keys]
		# Set up new combined key/entry
		new_key = (
			patient, ';'.join(sorted([x[1] for x in relevant_keys])), ';'.join(sorted([x[2] for x in relevant_keys]))
		)
		new_entry = []
		# Get medians for genomic coverage and for variant burdens
		for i in range(0, 34):
			new_entry.append(median([float(x[i]) for x in relevant_entries]))
		# Get processed MSI status
		if len([x[34] for x in relevant_entries if x[34] == '1']) >= 1:
			new_entry.append('1')
		else:
			new_entry.append('0')
		# Get medians for epitope burdens
		for i in range(35, 41):
			new_entry.append(median([float(x[i]) for x in relevant_entries]))
		# Get processed cancer stage
		new_entry.append(';'.join(list(set([x[41] for x in relevant_entries if x]))))
		# Get processed treatment/response values
		for i in range(42, 49):
			if len([x[i] for x in relevant_entries if x[i] == '1']) >= 1:
				new_entry.append('1')
			else:
				new_entry.append('0')
		# Get medians for PFS/OS/censoring
		for i in range(49, 52):
			new_entry.append(median([float(x[i]) for x in relevant_entries]))
		# Get processed event statuses
		for i in range(52, 56):
			if len([x[i] for x in relevant_entries if x[i] == '1']) >= 1:
				new_entry.append('1')
			else:
				new_entry.append('0')
		# Get median author-reported burdens
		for i in range(56, 59):
			new_entry.append(median([float(x[i]) for x in relevant_entries]))
		# Get processed disease/study
		new_entry.append(';'.join(list(set([x[59] for x in relevant_entries if x]))))
		new_entry.append(';'.join(list(set([x[60] for x in relevant_entries if x]))))
		# Replace dictionary entries
		assert len(new_entry) == len(relevant_entries[0])
		tumor_dict[new_key] = new_entry
		for key in relevant_keys:
			del tumor_dict[key]

	# Write to file
	with open(os.path.join(os.path.abspath(args.output_dir), 'immunotherapy_data_table.tsv'), 'w') as f:
		print(
				'\t'.join(
							[
								'Patient', 'Normal_ID', 'Tumor_ID', 'Coverage', 'Total_mutations', 'SNVs', 
								'Inframe_insertions', 'Inframe_deletions', 'Frameshift_insertions', 
								'Frameshift_deletions',  'Nonsynonymous_SNVs', 'Muse_variants', 'Muse_SNVs',
								'Muse_deletions', 'Muse_insertions', 'Mutect_variants', 'Mutect_SNVs', 
								'Mutect_deletions',  'Mutect_insertions', 'Pindel_variants', 'Pindel_SNVs', 
								'Pindel_deletions', 'Pindel_insertions', 'Radia_variants', 'Radia_SNVs', 
								'Radia_deletions', 'Radia_insertions', 'Somaticsniper_variants', 'Somaticsniper_SNVs', 
								'Somaticsniper_deletions', 'Somaticsniper_insertions', 'Varscan_variants', 
								'Varscan_SNVs', 'Varscan_deletions', 'Varscan_insertions', 'Tumor_HLA1_count', 
								'Tumor_HLA2_count', 'MSI_status', 'Total_unphased_neoepitopes', 
								'Total_comprehensive_neoepitopes', 'MHCnuggets_eps', 'MHCnuggets_ClassI_eps', 
								'MHCnuggets_ClassII_eps', 'Manuscript_binding_eps' 'Cancer_stage', 'aPD1_treatment', 
								'aPDL1_treatment', 'aCTLA4_treatment',  'Other_treatment', 'aCTLA4_response', 
								'aPD1_response',  'Combined_response', 'PFS', 'OS', 'Censoring_days', 'Vital_status', 
								'OS_event_status', 'PFS_event_status', 'Censoring_status', 
								'Original_nonsynonymous_mutations', 'Original_total_mutations', 'Original_neoantigens', 
								'Study', 'Disease', 
							]
						), 
				file=f
		)
		for tumor in tumor_dict:
			print('\t'.join(['\t'.join(list(tumor)), '\t'.join([str(x) for x in tumor_dict[tumor]])]), file=f)

