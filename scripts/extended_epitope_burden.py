#!/usr/bin/env python

from __future__ import print_function
from Bio.SubsMat import MatrixInfo
from collections import defaultdict
from mhcnuggets.src.predict import predict
from neoepiscope import bowtie_index
from neoepiscope import paths
from neoepiscope.transcript import get_transcripts_from_tree
from pyliftover import LiftOver
import argparse
import ast
import glob
import os
import pickle
import pysam
import re
import subprocess

def liftover_coordinates(neoepiscope_file, hg19_bed_path, hg38_bed_path, liftover, chain,
						 remove_bed_files=True):
	''' Creates a dictionary linking hg19 variants to their hg38 coordinates, if available

		neoepiscope_file: path to file containing neoepiscope output (comprehensive mode)
		hg19_bed_path: path to write hg19 BED file
		hg38_bed_path: path to write hg38 BED file
		liftover: path to liftOver executable
		chain: path to hg19 to hg38 chain file

		Return value: dictionary linking hg19 variant IDs to their hg38 coordinates
	'''
	# Establish dictionaries
	hg19_position_dict = {}
	hg38_position_dict = {}
	# Process neoepiscope output to ID variants
	with open(neoepiscope_file) as f:
		f.readline()
		f.readline()
		# Iterate through neoepitope entry lines
		for line in f:
			tokens = line.strip().split('\t')
			# Generate variant position ID
			name = ','.join([tokens[1], tokens[2], tokens[3], tokens[4], tokens[5]])
			if name not in hg19_position_dict:
				# Find end position and add info to dictionary
				start_pos = int(tokens[2]) - 1
				# Find end position
				if tokens[5] == 'V' or tokens[5] == 'D':
					end_pos = start_pos + len(tokens[3])
				elif tokens[5] == 'I':
					end_pos = start_pos + 1
				hg19_position_dict[name] = [tokens[1], str(start_pos), str(end_pos)]
	# Write variants to BED file
	with open(hg19_bed_path, 'w') as f:
		for variant in hg19_position_dict:
			out_line = hg19_position_dict[variant] + [variant]
			print('\t'.join(out_line), file=f)
	# Run liftOver
	subprocess.call([liftover, hg19_bed_path, chain, hg38_bed_path, '.'.join([hg38_bed_path, 'unmapped'])])
	# Process new BED file and add variants to dictionary
	with open(hg38_bed_path) as f:
		for line in f:
			tokens = line.strip().split('\t')
			hg38_position_dict[tokens[3]] = [tokens[0], int(tokens[1]), int(tokens[2])]
	# Remove BED files
	if remove_bed_files:
		for file_to_remove in [hg19_bed_path, hg38_bed_path, '.'.join([hg38_bed_path, 'unmapped'])]:
			os.remove(file_to_remove)
	# Return dictionary
	return hg38_position_dict

def process_epitope_file(neoepiscope_file, interval_dict, hg38_coordinates, binding_dict_search):
	''' Process neoepiscope output to generate allele list and dictionary

		neoepiscope_file: path to file containing neoepiscope output (comprehensive mode)
		interval_dict: dictionary linking contigs to IntervalTree objects
					   connecting genomic intervals to overlapping transcripts
		hg38_coordinates: dictionary linking hg19 variant IDs to their hg38 coordinates
		binding_dict_search: wildcard path for binding score dictionaries

		Return values: epitope_dict is a nested dictionary linking epitopes to mutations
					   they derive from, and a list of [paired nomal allele, list of hg19 transcripts
					   the epitope could derive from, list of hg38 transcripts the epitope could 
					   derive from, and binary binding score for each relevant allele (0 = doesn't 
					   bind allele, 1 = does bind allele)]; allele list is a list of HLA alleles for 
					   the patient; tx_set is a set of transcripts that neoepitopes might derive from
	'''
	epitope_dict = defaultdict(lambda: defaultdict(list))
	transcript_set = set()
	# Get binding scores and alleles
	binding_scores = {}
	binding_dicts = glob.glob(binding_dict_search)
	print(binding_dicts)
	allele_list = []
	for binding_dict in binding_dicts:
		allele = binding_dict.split('.')[-2]
		allele_list.append(allele)
		binding_scores[allele] = {}
		with open(binding_dict, 'rb') as p:
			d = pickle.load(p)
		for peptide in d:
			binding_scores[allele][peptide] = d[peptide]
	allele_list.sort()
	# Process peptides
	with open(neoepiscope_file) as f:
		f.readline()
		# Process header line to get relevant alleles
		header = f.readline().strip().split('\t')
		# Iterate through neoepitope entry lines
		for line in f:
			tokens = line.strip().split('\t')
			tx_list = tokens[9].split(';')
			# Extract paired normal and transcript list
			info_list = [tokens[7], tx_list]
			# Search for hg38 transcripts
			variant_id = ','.join([tokens[1], tokens[2], tokens[3], tokens[4], tokens[5]])
			if variant_id in hg38_coordinates:
				coords = hg38_coordinates[variant_id]
				overlapping_transcripts = get_transcripts_from_tree(coords[0], coords[1]+1, coords[2]+1, interval_dict)
			else:
				overlapping_transcripts = []
			info_list.append(overlapping_transcripts)
			# Iterate through binding scores
			for allele in allele_list:
				try:
					if float(binding_scores[allele][tokens[0]]) <= 500:
						# Epitope binds to this allele - tally as binder
						info_list.append(1)
					else:
						# Epitope doesn't bind to this allele - tally as non-binder
						info_list.append(0)
				except ValueError:
					info_list.append('NA')
			# Add info to mutation for epitope in dictionary
			epitope_dict[tokens[0]][(tokens[1], tokens[2], tokens[3], tokens[4], tokens[5])].append(info_list)
			for tx in tx_list:
				transcript_set.add(tx)
	# Return dictionary and allele list
	return epitope_dict, allele_list, transcript_set

def make_epitope_fasta(epitope_list, fasta):
	''' Produces fasta file containing all neoepitope sequences for a sample

		epitope_list: list of neoepitopes to perform a blast search on
		fasta: path to output fasta file

		No return value
	'''
	# Write epitopes from list to fasta file
	with open(fasta, 'w') as fh:
		for epitope in epitope_list:
			print(''.join(['> seq=', epitope]), file=fh)
			print(epitope, file=fh)

def run_blast(fasta, db, blastp, output, remove_fasta=True):
	''' Runs blastp

		fasta: path fasta containing sequences to blast against database
		db: path to database to blast against
		output: path to write blast output

		No return value
	'''
	# Run blast
	subprocess.call([blastp, '-outfmt', '6 qseqid sseqid length qstart qend sseq evalue', '-db', db, 
					 '-query', fasta, '-matrix', 'BLOSUM62', '-evalue', '200000', '-ungapped', 
					 '-comp_based_stats', 'F', '-out', output])
	if remove_fasta:
		# Remove fasta file
		os.remove(fasta)

def score_match(pair, matrix):
	''' Gives a score from a matrix for a given pair of sequences

		pair: pair to score (tuple)
		matrix: scoring matrix (matrix)

		Return value: score
	''' 
	# If the pair is not in the matrix, reverse it and get score
	if pair not in matrix:
		return matrix[(tuple(reversed(pair)))]
	else:
		return matrix[pair]                                                                                                            

def score_pairwise(seq1, seq2, matrix):
	''' For two sequences of equal length, score them given a matrix
		Score returned is relative to score of first sequence against itself

		seq1: first sequence (string)
		seq2: second sequence, of same length as first sequence (string)
		matrix: scoring matrix (matrix)

		Return value: score
	''' 
	self_score = 0
	score = 0
	# Calculate score of first sequence against itself for baseline
	for i in range(len(seq1)):
		pair = (seq1[i], seq1[i])
		self_score += score_match(pair, matrix)
	# Calulate score of first and second sequence
	for i in range(len(seq1)):
		pair = (seq1[i], seq2[i])
		score += score_match(pair, matrix)
	adj_score = score/self_score
	return adj_score

def process_blast(blast_results, remove_blast_results=True):
	''' Processes blast output to return a dicionary of results

		blast_results: path to file contain blast results
		data_dictionary: dictionary linking peptide_IDs to genes/transcripts

		Return value: blast_dict - dictionary that links epitopes to a list of
					  [match E value, set of transcripts it comes from, set of 
					  genes it comes from, match pepetide sequence]
	'''
	blast_dict = {}
	# Process through all blast results
	with open(blast_results, 'r') as fh:
		for line in fh:
			# Grab info from blast line
			tokens = line.strip().split('\t')
			epitope = tokens[0].split('=')[1]
			evalue = float(tokens[6])
			match_seq = tokens[5]	
			match_transcript = tokens[1].split('|')[0]
			match_gene = tokens[1].split('|')[1]
			# Check for presence of invalid characters in match seq
			invalids = ['B', 'J', 'O', 'U', 'X', 'Z', '*']
			invalid_matches = []
			for char in invalids:
				if char in match_seq:
					invalid_matches.append(char)
			# Check for match criteria
			if len(match_seq) == len(epitope) and len(invalid_matches) == 0:
				# Match is the right length and doesn't contain invalid characters
				if epitope not in blast_dict:
					# Epitope has no valid matches yet, add this one
					blast_dict[epitope] = [evalue, set(), set(), match_seq]
					blast_dict[epitope][1].add(match_transcript)
					blast_dict[epitope][2].add(match_gene)
				elif evalue < blast_dict[epitope][0]:
					# This is a better match for the epitope, replace old one
					blast_dict[epitope] = [evalue, set(), set(), match_seq]
					blast_dict[epitope][1].add(match_transcript)
					blast_dict[epitope][2].add(match_gene)
				elif evalue == blast_dict[epitope][0]:
					# Equally good match - check more criteria
					if match_seq == blast_dict[epitope][3]:
						# This sequence matches previous match(es), store tx/gene
						blast_dict[epitope][1].add(match_transcript)
						blast_dict[epitope][2].add(match_gene)
					else:
						# This is a new sequence - check more criteria
						if match_seq == epitope:
							# This match is identical to the epitope, replace old one
							blast_dict[epitope] = [evalue, set(), set(), match_seq]
							blast_dict[epitope][1].add(match_transcript)
							blast_dict[epitope][2].add(match_gene)
						else:
							# Check whether old or new match is more similar by adjusted BLOSUM62
							orig_ps = score_pairwise(epitope, blast_dict[epitope][3], MatrixInfo.blosum62)
							match_ps = score_pairwise(epitope, match_seq, MatrixInfo.blosum62)
							if match_ps > orig_ps:
								# This match is more similar, replace the old one
								blast_dict[epitope] = [evalue, set(), set(), match_seq]
								blast_dict[epitope][1].add(match_transcript)
								blast_dict[epitope][2].add(match_gene)
	# Remove blast results
	if remove_blast_results:
		os.remove(blast_results)
	# Return dictionary of matches
	return blast_dict

def get_normal_binding_scores(blast_dict, alleles, available_alleles, output_dir, pat_id, remove_files=True):
	''' Creates dictionary linking matched normal epitopes to binding scores for 
		different HLA alleles

		blast_dict: blast_dict - dictionary that links epitopes to a list of
					[match E value, set of transcripts it comes from, set of 
					genes it comes from, match pepetide sequence] (from process_blast())
		alleles: list of HLA alleles to use for binding predictions
		available_alleles: path to pickled dictionary describing available HLA alleles
						   for different binding affinity predictors
		output_dir: path to output directory for writing temporary files
		pat_id: patient identifier

		Return value: nested dictionary, where keys are matched normal epitopes and 
					  values are dictionaries, where keys are HLA alleles and values
					  are binding scores for that epitope/allele combo
	'''
	# Create list of temporary files to remove
	files_to_remove = []
	# Extract matched normal peptide sequences
	normal_epitopes = set()
	for epitope in blast_dict:
		normal_epitopes.add(blast_dict[epitope][3])
	# Load available alleles
	with open(available_alleles, 'rb') as allele_stream:
		avail_alleles = pickle.load(allele_stream)
	# Initialize dictionary
	normal_dict = defaultdict(dict)
	for hla in alleles:
		# Determine if allele is valid for mhcnuggets
		if hla in avail_alleles["mhcnuggets_mhcI"]:
			# Class I allele
			allele_class = "I"
			max_length = 15
		elif hla in avail_alleles["mhcnuggets_mhcII"]:
			# Class II allele
			allele_class = "II"
			max_length = 30
		else:
			# Not a valid allele
			continue
		# Write relevant peptides to file
		peptide_file = os.path.join(output_dir, ''.join([pat_id, '.mhc.', hla, '.csv']))
		files_to_remove.append(peptide_file)
		with open(peptide_file, 'w') as f:
			for sequence in normal_epitopes:
				if len(sequence) <= max_length:
					print(sequence, file=f)
		# Run binding predictions
		mhc_out = os.path.join(output_dir, ''.join([pat_id, '.mhc.', hla, '.out']))
		files_to_remove.append(mhc_out)
		predict(class_=allele_class, peptides_path=peptide_file, mhc=hla, output=mhc_out)
		# Process mhcnuggets results
		score_dict = {}
		with open(mhc_out) as f:
			f.readline()
			for line in f:
				tokens = line.strip().split(',')
				score_dict[tokens[0]] = tokens[1]
		# Store score for each epitope if available
		for sequence in normal_epitopes:
			if sequence in score_dict:
				normal_dict[sequence][hla] = float(score_dict[sequence])
	# Remove temporary files
	if remove_files:
		for file_to_remove in files_to_remove:
			os.remove(file_to_remove)
	# Return dictionary
	return normal_dict

### Taken with modifications from https://www.biostars.org/p/306041/
def read_pair_generator(bam, contig):
	"""
	Generate read pairs in a BAM file or within a region string.
	Reads are added to read_dict until a pair is found.

	bam: path to BAM file
	contig: name of contig to process

	Return value: iterator of read1, read2 for each paired read
	"""
	read_dict = defaultdict(lambda: [None, None])
	for read in bam.fetch(contig=contig):
		if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
			# Ignore reads that aren't mapped as pairs or from primary alignments
			continue
		if read.query[-2:] in ['.1', '.2']:
			# Strip '.1' or '.2' from query name
			qname = read.query_name[:-2]
		else:
			qname = read.query_name
		if qname not in read_dict:
			# Establishing new read pair
			if read.is_read1:
				read_dict[qname][0] = read
			else:
				read_dict[qname][1] = read
		else:
			# Finish established read pair
			if read.is_read1:
				yield read, read_dict[qname][1]
			else:
				yield read_dict[qname][0], read
			del read_dict[qname]

def create_bed_file(bed_path, transcripts, cds_dict):
	""" Creates BED file of genome coordinates for relevant transcripts

		bed_path: path to write BED file
		transcripts: set of relevant transcript IDs
		cds_dict: dictionary linking transcript IDs to GTF file entry info

		Return value: None
	"""
	with open(bed_path, 'w') as f:
		for tx in transcripts:
			for cds in cds_dict[tx]:
				out_line = [cds[0], str(cds[2]), str(cds[3]+1)]
				print('\t'.join(out_line), file=f)

# Taken from Rail-RNA: https://github.com/nellore/rail
def parsed_md(md):
	""" Divides an MD string up by boundaries between ^, letters, and numbers
		md: an MD string (example: 33A^CC).
		Return value: MD string split by boundaries described above.
	"""
	md_to_parse = []
	md_group = [md[0]]
	for i, char in enumerate(md):
		if i == 0: continue
		if (re.match('[A-Za-z]', char) is not None) \
			!= (re.match('[A-Za-z]', md[i-1]) is not None) or \
			(re.match('[0-9]', char) is not None) \
			!= (re.match('[0-9]', md[i-1]) is not None):
			if md_group:
				md_to_parse.append(''.join(md_group))
			md_group = [char]
		else:
			md_group.append(char)
	if md_group:
		md_to_parse.append(''.join(md_group))
	return [char for char in md_to_parse if char != '0']

# Taken from Rail-RNA: https://github.com/nellore/rail
def indels_junctions_exons_mismatches(cigar, md, pos, seq, drop_deletions=False, junctions_only=False):
	""" Finds indels, junctions, exons, mismatches from CIGAR, MD string, POS

		cigar: CIGAR string
		md: MD:Z string
		pos: position of first aligned base
		seq: read sequence
		drop_deletions: drops deletions from coverage vectors iff True
		junctions_only: does not populate mismatch list
		Return value: tuple (insertions, deletions, junctions, exons, mismatches).
		Insertions is a list of tuples (last genomic position before insertion, 
			string of inserted bases). Deletions is a list of tuples (first genomic 
			position of deletion, string of deleted bases). Junctions is a list
			of tuples (intron start position (inclusive), intron end position (exclusive),
			left_diplacement, right_displacement). Exons is a list of tuples (exon start 
			position (inclusive), exon end position (exclusive)). Mismatches is a list
			of tuples (genomic position of mismatch, read base)
	"""
	insertions, deletions, junctions, exons, mismatches = [], [], [], [], []
	cigar = re.split(r'([MINDS])', cigar)[:-1]
	md = parsed_md(md)
	seq_size = len(seq)
	cigar_chars, cigar_sizes = [], []
	cigar_index, md_index, seq_index = 0, 0, 0
	max_cigar_index = len(cigar)
	while cigar_index != max_cigar_index:
		if cigar[cigar_index] == 0:
			cigar_index += 2
			continue
		if cigar[cigar_index+1] == 'M':
			aligned_base_cap = int(cigar[cigar_index])
			aligned_bases = 0
			while True:
				try:
					aligned_bases += int(md[md_index])
					if aligned_bases <= aligned_base_cap:
						md_index += 1
				except ValueError:
					# Not an int, but should not have reached a deletion
					assert md[md_index] != '^', '\n'.join(['cigar and md:', ''.join(cigar), ''.join(md)])
					if not junctions_only:
						mismatches.append((pos + aligned_bases, seq[seq_index + aligned_bases]))
					correction_length = len(md[md_index])
					m_length = aligned_base_cap - aligned_bases
					if correction_length > m_length:
						md[md_index] = md[md_index][:m_length]
						aligned_bases = aligned_base_cap
					else:
						aligned_bases += correction_length
						md_index += 1
				if aligned_bases > aligned_base_cap:
					md[md_index] = aligned_bases - aligned_base_cap
					break
				elif aligned_bases == aligned_base_cap:
					break
			# Add exon
			exons.append((pos, pos + aligned_base_cap))
			pos += aligned_base_cap
			seq_index += aligned_base_cap
		elif cigar[cigar_index+1] == 'N':
			skip_increment = int(cigar[cigar_index])
			# Add junction
			junctions.append((pos, pos + skip_increment, seq_index, seq_size - seq_index))
			# Skip region of reference
			pos += skip_increment
		elif cigar[cigar_index+1] == 'I':
			# Insertion
			insert_size = int(cigar[cigar_index])
			insertions.append((pos - 1, seq[seq_index:seq_index+insert_size]))
			seq_index += insert_size
		elif cigar[cigar_index+1] == 'D':
			assert md[md_index] == '^', '\n'.join(['cigar and md:', ''.join(cigar), ''.join(md)])
			# Deletion
			delete_size = int(cigar[cigar_index])
			md_delete_size = len(md[md_index+1])
			assert md_delete_size >= delete_size
			deletions.append((pos, md[md_index+1][:delete_size]))
			if not drop_deletions: exons.append((pos, pos + delete_size))
			if md_delete_size > delete_size:
				md[md_index+1] = md[md_index+1][delete_size:]
			else:
				md_index += 2
			# Skip deleted part of reference
			pos += delete_size
		else:
			assert cigar[cigar_index+1] == 'S'
			seq_index += int(cigar[cigar_index])
		cigar_index += 2
	# Merge exonic chunks/deletions; insertions/junctions could have chopped them up
	new_exons = []
	last_exon = exons[0]
	for exon in exons[1:]:
		if exon[0] == last_exon[1]:
			# Merge ECs
			last_exon = (last_exon[0], exon[1])
		else:
			# Push last exon to new exon list
			new_exons.append(last_exon)
			last_exon = exon
	new_exons.append(last_exon)
	return insertions, deletions, junctions, new_exons, mismatches

def get_expressed_transcripts(bam, reference_index, interval_dict, bed_path, remove_files=True):
	''' Gets transcripts that are expressed by at least one read pair

		bam: path to RNA-seq BAM file
		reference_index: bowtie reference index
		interval_dict: dictionary linking contigs to IntervalTree objects
					   connecting genomic intervals to overlapping transcripts
		bed_path: path to BED file describing coordinates of transcripts of interest

		Return value: set of transcript IDs that are expressed
	'''
	# Set up BAM reader/transcript set
	transcript_set = set()
	# Produce reduced BAM file
	new_bam = bam.replace('.bam', '.reduced.bam')
	subprocess.call(['samtools', 'view', '-b', '-L', bed_path, '-o', new_bam, bam])
	subprocess.call(['samtools', 'index', new_bam, ''.join([new_bam, '.bai'])])
	# Process BAM file
	bam_reader = pysam.AlignmentFile(new_bam, 'rb')
	for contig in reference_index.recs.keys():
		for read1, read2 in read_pair_generator(bam_reader, contig):
			# Get read info
			r1_tokens = str(read1).split('\t')
			r2_tokens = str(read2).split('\t')
			r1_tags = ast.literal_eval(r1_tokens[11])
			r2_tags = ast.literal_eval(r2_tokens[11])
			r1_md = [x[1] for x in r1_tags if x[0] == 'MD'][0]
			r2_md = [x[1] for x in r2_tags if x[0] == 'MD'][0]
			r1_cigar = r1_tokens[5]
			r2_cigar = r2_tokens[5]
			r1_pos = int(r1_tokens[3])+1
			r2_pos = int(r2_tokens[3])+1
			r1_seq = r1_tokens[9]
			r2_seq = r2_tokens[9]
			# Extract data from each read
			r1_insertions, r1_deletions, r1_junctions, r1_exons, r1_mismatches = indels_junctions_exons_mismatches(r1_cigar, r1_md, r1_pos, r1_seq)
			r2_insertions, r2_deletions, r2_junctions, r2_exons, r2_mismatches = indels_junctions_exons_mismatches(r2_cigar, r2_md, r2_pos, r2_seq)
			# Find genomic intervals covered, accounting for junctions
			intervals = []
			for e in r1_exons + r2_exons:
				intervals.extend([e[0], e[1]+1])
			# Find and store transcripts that overlap intervals
			for i in range(0, len(intervals), 2):
				overlapping_transcripts = interval_dict[contig].overlap(intervals[i], intervals[i+1])
				for tx in [x.data for x in overlapping_transcripts]:
					transcript_set.add(tx)
	# Remove temporary files
	if remove_files:
		for file_to_remove in [new_bam, ''.join([new_bam, '.bai']), bed_path]:
			os.remove(file_to_remove)
	# Return trancript IDs
	return transcript_set

def process_tcga_expression(tcga_dict, disease):
	''' Produces set of transcripts expressed among TCGA samples
		for a given cancer type

		tcga_dict: TCGA expression dictionary - nested dictionary where 
				   keys are TCGA disease types and values are sets of 
				   transcript IDs for transcripts that are expressed at 
				   a TPM of at least 1 for the 75th percentile expression 
				   value for that transcript among tumor samples
		disease: disease type, one of [melanoma, NSCLC, colon, 
				 endometrial, thyroid, prostate, RCC]

		Return value: set of expressed transcripts
	'''
	# ID TCGA disease
	if disease == 'melanoma':
		disease_list = ['SKCM']
	elif disease == 'NSCLC':
		disease_list = ['LUAD', 'LUSC']
	elif disease == 'colon':
		disease_list = ['COAD']
	elif disease == 'endometrial':
		disease_list = ['UCEC']
	elif disease == 'thyroid':
		disease_list = ['THCA']
	elif disease == 'prostate':
		disease_list = ['PRAD']
	elif disease == 'RCC':
		disease_list = ['KIRC']
	# Identify TCGA-expressed transcripts
	tx_set = set()
	for disease_id in disease_list:
		for transcript in tcga_dict[disease_id]:
			# transcript is expressed
			tx_set.add(transcript)
	# Return transcript set
	return tx_set

def extended_burden(patient, tumor, rna_id, epitope_dir, output_dir, blastp, db, rna_dir, 
					tcga_dict_path, disease, available_alleles, chain_file, liftover, binding_dir):
	''' Generates extended burden information for each neoepitope

		patient: patient identifier (string)
		tumor: tumor sample identifier (string)
		rna_id: tumor RNA sample identifier (string)
		epitope_dir: path to directory with neoepitope predictions from neoepiscope (string)
		output_dir: path to directory for script to write files to (string)
		blastp: path to blastp executable (string)
		db: path to human protein database (string)
		rna_dir: path to directory with with RNA alignment subdirectories by tumor (string)
		tcga_dict_path: path to dictionary with TCGA transcript expression values (string)
		disease: disease of sample (string)
		available alleles: path to neoepiscope's available alleles dictionary (string)
		chain_file: path to chain file for liftover (string)
		liftover: path to liftover executable (string)
		binding_dir: path to directory with binding affinity dictionaries (string)

		No return value, output file written.
	'''
	# Establish neoepitope directory/file
	epitope_file = os.path.join(epitope_dir, ''.join([patient, '.', tumor, '.neoepiscope.comprehensive.out']))
	# Find liftover coordinates
	hg19_bed = os.path.join(output_dir, ''.join([patient, '.', tumor, '.hg19.bed']))
	hg38_bed = os.path.join(output_dir, ''.join([patient, '.', tumor, '.hg38.bed']))
	hg38_coordinates = liftover_coordinates(epitope_file, hg19_bed, hg38_bed, liftover, chain_file)
	# Process neoepiscope data
	with open(os.path.join(paths.gencode_v29, 'intervals_to_transcript.pickle'), 'rb') as interval_stream:
		v29_intervals = pickle.load(interval_stream)
	wildcard_path = ''.join([binding_dir, patient, '.', tumor, '.mhcflurry.*.pickle'])
	epitope_dict, allele_list, relevant_transcripts = process_epitope_file(epitope_file, v29_intervals, hg38_coordinates, wildcard_path)
	# Produce FASTA file for blast
	fasta_path = os.path.join(output_dir, ''.join([patient, '.', tumor, '.epitopes.fasta']))
	make_epitope_fasta(epitope_dict.keys(), fasta_path)
	# Run blast
	blast_output = os.path.join(output_dir, ''.join([patient, '.', tumor, '.blast.out']))
	run_blast(fasta_path, db, blastp, blast_output)
	# Process blast results
	blast_data = process_blast(blast_output)
	# Get matched normal binding scores
	normal_binding = get_normal_binding_scores(blast_data, allele_list, available_alleles, 
											   output_dir, '.'.join([patient, tumor]))
	# Process RNA-seq data, if available
	rna_bam = os.path.join(rna_dir, rna_id, ''.join([rna_id, '.sorted.bam']))
	if os.path.isfile(rna_bam):
		# Patient specific expression data available - load annotation info
		with open(os.path.join(paths.gencode_v19, 'intervals_to_transcript.pickle'), 'rb') as interval_stream:
			interval_dict = pickle.load(interval_stream)
		with open(os.path.join(paths.gencode_v19, 'transcript_to_CDS.pickle'), 'rb') as cds_stream:
			cds_dict = pickle.load(cds_stream)
		reference_index = bowtie_index.BowtieIndexReference(paths.bowtie_hg19)
		bed_file = os.path.join(output_dir, ''.join([patient, '.', tumor, '.bed']))
		create_bed_file(bed_file, relevant_transcripts, cds_dict)
		# Identify expressed transcripts
		expressed_transcripts = get_expressed_transcripts(rna_bam, reference_index, interval_dict, bed_file)
	else:
		# No patient-specific expression data available
		expressed_transcripts = None
	# Process TCGA expression data
	with open(tcga_dict_path, 'rb') as f:
		tcga_dict = pickle.load(f)
	tcga_expressed_transcripts = process_tcga_expression(tcga_dict, disease)
	# Write output
	output_file = os.path.join(output_dir, ''.join([patient, '.', tumor, '.extended_epitope_burden.tsv']))
	with open(output_file, 'w') as o:
		header = ['Patient_ID', 'Tumor_ID', 'Epitope', 'Mutation', 'Mutation_type', 'Paired_normal', 'Paired_mismatches', 
				  'Blast_normal', 'Blast_mismatches', 'Anchor_mismatches', 'Transcript_IDs', 'Match_transcripts', 
				  'HLA_binders', 'Alleles_bound', 'Matched_binders', 'Match_alleles', 'Patient_expression', 
				  'Patient_transcripts', 'TCGA_expression', 'TCGA_transcripts']
		print('\t'.join(header), file=o)
		output_lines = set()
		for epitope in epitope_dict:
			# Determine matched normal data from blast, if available
			if epitope in blast_data:
				blast_normal = blast_data[epitope][3]
				mismatches = str(sum(x1!=x2 for x1,x2 in zip(epitope,blast_normal)))
				# Get anchor residue mismatch count
				if mismatches != '0':
					anchor_count = 0
					if len(epitope) <= 11:
						for i in [1, 4, 5, -1]:
							if epitope[i] != blast_normal[i]:
								anchor_count += 1
					else:
						for i in [0, 3, 5, 8]:
							if epitope[i] != blast_normal[i]:
								anchor_count += 1
					anchor_mismatches = str(anchor_count)
				else:
					anchor_mismatches = '0'
				normal_transcripts = ','.join(blast_data[epitope][1])
				normal_genes = blast_data[epitope][1]
				# Get matched normal binding scores
				normal_alleles = [x for x in normal_binding[blast_normal] if normal_binding[blast_normal][x] <= 500.0]
				if len(normal_alleles) > 0:
					matched_binders = str(len(normal_alleles))
					matched_alleles = ','.join(normal_alleles)
				else:
					matched_binders = '0'
					matched_alleles = 'NA'
			else:
				blast_normal = 'NA'
				mismatches = 'NA'
				normal_transcripts = 'NA'
				normal_genes = 'NA'
				matched_binders = 'NA'
				matched_alleles = 'NA'
				anchor_mismatches = 'NA'
			# Process through mutations associated with epitope
			for mutation in epitope_dict[epitope]:
				mut_id = ','.join(mutation[0:-1])
				mut_type = mutation[-1]
				for ep_data in epitope_dict[epitope][mutation]:
					paired_normal =  ep_data[0]
					if paired_normal != 'NA':
						paired_mismatches = str(sum(x1!=x2 for x1,x2 in zip(epitope,paired_normal)))
					else:
						paired_mismatches = 'NA'
					transcripts = ep_data[1]
					# Establish binding alleles
					binding_alleles = []
					for i in range(3, len(ep_data)):
						if ep_data[i] == 1:
							binding_alleles.append(allele_list[i-3])
					if len(binding_alleles) > 0:
						binders = str(len(binding_alleles))
						alleles = ','.join(binding_alleles)
					else:
						binders = '0'
						alleles = 'NA'
					# Establish transcript expression if RNA-seq is available
					patient_transcripts = []
					for tx in transcripts:
						if expressed_transcripts is not None:
							if tx in expressed_transcripts:
								patient_transcripts.append(tx) 
					if len(patient_transcripts) > 0:
						patient_expression = str(len(patient_transcripts))
						patient_tx_ids = ','.join(patient_transcripts) 
					elif expressed_transcripts is not None:
						patient_expression = '0'
						patient_tx_ids = 'NA'
					else:
						patient_expression = 'NA'
						patient_tx_ids = 'NA'
					# Establish TCGA expression if liftover was successful
					tcga_transcripts = []
					hg38_transcripts = ep_data[2]
					for tx in hg38_transcripts:
						if tx in tcga_expressed_transcripts:
							tcga_transcripts.append(tx)
					if len(tcga_transcripts) > 0:
						tcga_expression = str(len(tcga_transcripts))
						tcga_tx_ids = ','.join(tcga_transcripts)
					else:
						tcga_expression = '0'
						tcga_tx_ids = 'NA'
					out_line = [patient, tumor, epitope, mut_id, mut_type, paired_normal, paired_mismatches,
								blast_normal, mismatches, anchor_mismatches, ','.join(transcripts), normal_transcripts, 
								binders, alleles, matched_binders, matched_alleles, patient_expression, patient_tx_ids,
								tcga_expression, tcga_tx_ids]
					print('\t'.join(out_line), file=o)

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--patient', type=str, required=True,
					help='patient ID'
				)
parser.add_argument('-w', '--wes', type=str, required=True,
					help='tumor WES ID'
				)
parser.add_argument('-r', '--rna', type=str, required=True,
					help='tumor RNA-seq ID'
				)
parser.add_argument('-n', '--neoepitope-directory', type=str, required=True,
					help='path to directory containing neoepiscope output (comprehensive mode)'
				)
parser.add_argument('-o', '--output-directory', type=str, required=True,
					help='path to output directory'
				)
parser.add_argument('-e', '--blastp-executable', type=str, required=True,
					help='path to blastp executable'
				)
parser.add_argument('-b', '--blastp-database', type=str, required=True,
					help='path to blastp peptide database'
				)
parser.add_argument('-a', '--rna-alignments', type=str, required=True,
					help='path to directory containing RNA-seq alignment subdirectories'
				)
parser.add_argument('-t', '--tcga-expression', type=str, required=True,
					help='path to pickled dictionary with TCGA expression info'
				)
parser.add_argument('-s', '--disease-site', type=str, required=True,
					help='disease site (melanoma, NSCLC, colon, endometrial, thyroid, prostate, RCC)'
				)
parser.add_argument('-d', '--available-alleles-dict', type=str, required=True,
					help='path to pickled dictionary with available alleles for MHC binding predictors'
				)
parser.add_argument('-l', '--liftover', type=str, required=True,
					help='path to liftOver executable'
				)
parser.add_argument('-c', '--chain-file', type=str, required=True,
					help='path to hg19-to-hg38 chain file'
				)
parser.add_argument('-m', '--mhc-binding-score-dir', type=str, required=True,
					help='path to MHC binding score dict directory'
				)
args = parser.parse_args()

extended_burden(args.patient, args.wes, args.rna, args.neoepitope_directory, args.output_directory, 
				args.blastp_executable, args.blastp_database, args.rna_alignments,
				args.tcga_expression, args.disease_site, args.available_alleles_dict, 
				args.chain_file, args.liftover, args.mhc_binding_score_dir)

