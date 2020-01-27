#!/usr/bin/env python

from __future__ import print_function
import argparse
import bisect
import os
import pickle
from neoepiscope import paths, bowtie_index
from neoepiscope.transcript import Transcript, kmerize_peptide, seq_to_peptide
from neoepiscope.binding_scores import get_affinity_mhcnuggets
from collections import defaultdict
from string import maketrans
from bisect import bisect_left
from numpy import median


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
    parser.add_argument("--hla-type-dir", "-t",
                        type=str,
                        required=True,
                        help="path to directory with optitype and seq2hla output")
    parser.add_argument("--introns-to-transcripts", "-i",
                        type=str,
                        required=True,
                        help="path to file linking introns to transcripts")
    parser.add_argument("--filtered-outliers", "-f",
                        type=str,
                        required=True,
                        help="path to file with filtered outlier introns")
    parser.add_argument("--hla-reference-fasta", "-r",
    					type=str, required=True, 
    					help="path to HLA reference FASTA file from Optitype")
    args = parser.parse_args()

	# Load reference data
	with open(os.path.join(paths.gencode_v19, "intervals_to_transcript.pickle"), "rb") as interval_stream:
		interval_dict = pickle.load(interval_stream)
	with open(os.path.join(paths.gencode_v19, "transcript_to_CDS.pickle"), "rb") as cds_stream:
		cds_dict = pickle.load(cds_stream)
	reference_index = bowtie_index.BowtieIndexReference(paths.bowtie_hg19)
	revcomp_translation_table = maketrans("ATCG", "TAGC")

	# Get HLA decoder
	hla_id_linker = {}
	hla_fasta = os.path.abspath(args.hla_reference_fasta)
	for record in SeqIO.parse(hla_fasta, 'fasta'):
		header = str(record.description).split()
		hla_id_linker[header[0]] = ':'.join(header[1].split(':')[0:2])

	# Get patient info
	patients = defaultdict(set)
	allele_dict = {}
	with open(os.path.abspath(args.manifest)) as f:
		for line in f:
			tokens = line.strip().split('\t')
			if tokens[2] != 'NA' and tokens[3] != 'NA':
				allele_set = set()
				#Get class 1 alleles
				with open(os.path.join(os.path.abspath(args.hla_type_dir), ''.join([tokens[2], '_result.tsv']))) as h:
					h.readline()
					hla_tokens = h.readline().strip().split('\t')
					for i in range(1, 7):
						if 'HLA' in hla_tokens[i]:
							hla = hla_id_linker[hla_tokens[i]]
						else:
							hla = '-'.join(['HLA', hla_tokens[i]])
						allele_set.add(hla)
				# Get class 2 alleles
				with open(os.path.join(os.path.abspath(args.hla_type_dir), ''.join([tokens[2], '-ClassII.HLAgenotype4digits']))) as h:
					h.readline()
					for line in h:
						hla_tokens = line.strip('\n').split('\t')
						for i in [1, 3]:
							allele_opt = hla_tokens[1].split(',')
							for allele in allele_opt:
								if allele != 'no':
									hla = ''.join(['HLA-', allele.strip("'")])
									allele_set.add(hla)
				# Add info to dictionaries
				patients[tokens[0]].add((tokens[2], tokens[3]))
				allele_dict[tokens[3]] = allele_set

	# Process intron to transcript linker file
	intron_to_transcripts = {}
	with open(os.path.abspath(args.introns_to_transcripts)) as f:
		f.readline()
		for line in f:
			tokens = line.strip().split('\t')
			if tokens[0] not in intron_to_transcripts:
				intron_to_transcripts[tokens[3]] = [[tokens[1]], [tokens[2]], tokens[0], tokens[4]]
			else:
				intron_to_transcripts[tokens[3]][0].append(tokens[1])
				intron_to_transcripts[tokens[3]][1].append(tokens[2])

	# Establish dictionaries
	count_dict = defaultdict(int)
	epitope_dict = defaultdict(set)
	binding_epitope_dict = defaultdict(set)

	# Process intron retention data
	with open(os.path.abspath(args.filtered_outliers)) as f:	
		# Establish patient list
		patient_list = f.readline().strip().split('\t')[1:]
		for line in f:
			tokens = line.strip().split('\t')
			# ID intron chrom, start, end
			intron = intron_to_transcripts[tokens[0]][2].split(':')
			chromosome = intron[0]
			start_pos = int(intron[1].split('-')[0]) + 1
			end_pos = int(intron[1].split('-')[1])
			intron_length = end_pos - start_pos + 1
			# Eliminate introns of length 0
			if intron_length == 0:
				continue
			# Find relevant transcripts
			transcripts = intron_to_transcripts[tokens[0]][0]
			relevant_transcripts = [tx for tx in transcripts if tx in cds_dict]
			if len(relevant_transcripts) > 0:
				# Tally intron for relevant patients 
				for i in range(1, len(tokens)):
					if tokens[i] == '1':
						count_dict[patient_list[i-1]] += 1
				epitopes = set()
				for tx in relevant_transcripts:
					# Set up transcript object
					cds_lines = [[str(chrom), ".", seq_type, str(start), str(end), ".", strand]
					             for (chrom, seq_type, start, end, strand, tx_type) in cds_dict[tx]]
					tx_object = Transcript(reference_index, cds_lines, tx)
					strand = 1 - tx_object.rev_strand * 2
					# Reverse strand
					if tx_object.rev_strand:
						# Eliminate introns outside coding region
						if tx_object.stop_codon is not None:
							if start_pos < tx_object.stop_codon:
								continue
						if tx_object.start_codon is not None:
							if end_pos > tx_object.start_codon:
								continue
						# Check for full intron
						if (start_pos - 2) in tx_object.intervals and (end_pos - 1) in tx_object.intervals:
							# Find left genomic location
							left_frame = tx_object.reading_frame(end_pos+1)
							if left_frame is None:
								continue
							elif left_frame == 0:
								left_index = end_pos + 1
							elif left_frame == 1:
								left_index = end_pos + 2
							elif left_frame == 2:
								left_index = end_pos
							# Adjust starting postion to get sequence divisble by 3
							adj_length = (((left_index - end_pos) + intron_length) % 3)
							if adj_length == 0:
								start = start_pos - 1
								offset = (left_index - end_pos) + intron_length
							elif adj_length == 1:
								start = start_pos - 3
								offset = (left_index - end_pos) + intron_length + 2
							elif adj_length == 2:
								start = start_pos - 2
								offset = (left_index - end_pos) + intron_length + 1
							# Grab intron sequence
							intron_sequence = reference_index.get_stretch(chromosome, 
																		  start, 
																		  offset)[::-1].translate(revcomp_translation_table)
							assert (len(intron_sequence) % 3) == 0
						else:
							# Find position of last base in upstream exon
							index = bisect_left(tx_object.intervals, start_pos - 1)
							exon_pos = tx_object.intervals[index]+2
							exon_pos2 = tx_object.intervals[index-1]+1
							# Get left genomic position
							left_frame = tx_object.reading_frame(exon_pos)
							if left_frame is None:
								continue
							elif left_frame == 0:
								left_index = exon_pos
								offset = 1
							elif left_frame == 1:
								left_index = exon_pos
								offset = 2
							elif left_frame == 2:
								left_index = exon_pos
								offset = 0
							# Determine downstream exon sequence needed
							adj_length = (offset + intron_length) % 3
							if adj_length == 0:
								right_index = 0
							elif adj_length == 1:
								right_index = exon_pos2 - 1
								offset2 = 2
							elif adj_length == 2:
								right_index = exon_pos2
								offset2 = 1
							# Grab intron sequence
							intron_seqs = []
							# Get upstream exon sequence if necessary
							if offset > 0:
								pre_seq = reference_index.get_stretch(chromosome, 
																	  left_index - 1, 
																	  offset)[::-1].translate(revcomp_translation_table)
								intron_seqs.append(pre_seq)
							# Get main sequence
							main_seq = reference_index.get_stretch(chromosome, 
																   start_pos - 1, 
																   intron_length)[::-1].translate(revcomp_translation_table)
							intron_seqs.append(main_seq)
							# Grab downstream intron sequence if necessary
							if right_index > 0:
								post_seq = reference_index.get_stretch(chromosome, 
																	   right_index - 1, 
																	   offset2)[::-1].translate(revcomp_translation_table)
								intron_seqs.append(post_seq)
							# Combine sequences
							intron_sequence = ''.join(intron_seqs)
							assert (len(intron_sequence) % 3) == 0
					# Forward strand
					else:
						# Eliminate introns outside coding region
						if tx_object.stop_codon is not None:
							if start_pos > tx_object.stop_codon:
								continue
						if tx_object.start_codon is not None:
							if end_pos < tx_object.start_codon:
								continue
						# Check for full intron
						if (start_pos - 2) in tx_object.intervals and (end_pos - 1) in tx_object.intervals:
							# Find left genomic location
							left_frame = tx_object.reading_frame(start_pos-1)
							if left_frame is None:
								continue
							elif left_frame == 0:
								left_index = start_pos - 1
							elif left_frame == 1:
								left_index = start_pos - 2
							elif left_frame == 2:
								left_index = start_pos
							# Adjust offset to get sequence divisble by 3
							adj_length = (((start_pos - left_index) + intron_length) % 3)
							if adj_length == 0:
								offset = (start_pos - left_index) + intron_length
							elif adj_length == 1:
								offset = (start_pos - left_index) + intron_length + 2
							elif adj_length == 2:
								offset = (start_pos - left_index) + intron_length + 1
							# Grab intron sequence
							intron_sequence = reference_index.get_stretch(chromosome, 
																		  left_index - 1, 
																		  offset)
							assert (len(intron_sequence) % 3) == 0
						else:
							# Find position of last base in upstream exon
							index = bisect_left(tx_object.intervals, start_pos - 1)
							exon_pos = tx_object.intervals[index-1]+1
							exon_pos2 = tx_object.intervals[index]+2
							# Get left genomic position
							left_frame = tx_object.reading_frame(exon_pos)
							if left_frame is None:
								continue
							elif left_frame == 0:
								left_index = exon_pos
								grab_seq = True
							elif left_frame == 1:
								left_index = exon_pos - 1
								grab_seq = True
							elif left_frame == 2:
								left_index = exon_pos + 1
								grab_seq = False
							# Determine downstream exon sequence needed
							adj_length = ((exon_pos - left_index + 1) + intron_length) % 3
							if adj_length == 0:
								right_index = 0
							elif adj_length == 1:
								right_index = exon_pos2
								offset = 2
							elif adj_length == 2:
								right_index = exon_pos2
								offset = 1
							# Grab intron sequence
							intron_seqs = []
							# Get upstream exon sequence if necessary
							if grab_seq:
								pre_seq = reference_index.get_stretch(chromosome, 
																	  left_index - 1, 
																	  (exon_pos - left_index) + 1)
								intron_seqs.append(pre_seq)
							# Get main sequence
							main_seq = reference_index.get_stretch(chromosome, 
																   start_pos - 1, 
																   intron_length)
							intron_seqs.append(main_seq)
							# Grab downstream intron sequence if necessary
							if right_index > 0:
								post_seq = reference_index.get_stretch(chromosome, 
																	  right_index - 1, 
																	  offset)
								intron_seqs.append(post_seq)
							# Combine sequences
							intron_sequence = ''.join(intron_seqs)
							assert (len(intron_sequence) % 3) == 0
					# Get peptides + binding affinities	
					peptides = kmerize_peptide(seq_to_peptide(intron_sequence, reverse_strand=False), min_size=8, max_size=24)
					for i in range(1, len(tokens)):
						for pep in peptides:
							epitope_dict[patient_list[i-1]].add(pep)
		# Get binding affinities for each relevant patient
		for patient in patient_list:
			peptides = list(epitope_dict[patient])
			for allele in allele_dict[patient]:
				affinities = get_affinity_mhcnuggets(peptides, allele, '2')
				for j in range(len(affinities)):
					epitope_dict[patient].add(affinities[j][0])
					if affinities[j][1] != 'NA':
						if float(affinities[j][1]) <= 500.0:
							binding_epitope_dict[patient].add(affinities[j][0])

	# Write output file
	with open(os.path.join(os.path.abspath(args.output_dir), 'full_intron_retention_burden.tsv'), 'w') as o:
		header = ['Patient', 'Tumor_ID', 'Intron_burden', 'Intron_epitope_burden', 'Binding_intron_epitope_burden']
		print('\t'.join(header), file=o)
		for pat in patients:
			tumor = ';'.join(sorted([x[0] for x in patients[pat]]))
			rna_ids = [x[1] for x in patients[pat]]
			intron_counts = str(median([count_dict[x] for x in rna_ids]))
			ep_counts = str(median([len(epitope_dict[x]) for x in rna_ids]))
			binding_counts = str(median([len(binding_epitope_dict[x]) for x in rna_ids]))
			out_line = [pat, tumor, intron_counts, ep_counts, binding_counts]
			print('\t'.join(out_line), file=o)

