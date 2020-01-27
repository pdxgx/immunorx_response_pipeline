#!/usr/bin/env python

from __future__ import print_function
from astropy.stats import median_absolute_deviation
from collections import defaultdict
from numpy import median
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--read-counts", "-r",
                        type=str,
                        required=True,
                        help="path to input file with intron retention read counts")
    parser.add_argument("--output-file", "-o",
                        type=str,
                        required=True,
                        help="path to output file")
    parser.add_argument("--introns-to-transcripts", "-i",
                        type=str,
                        required=True,
                        help="path to file linking introns to transcripts")
    args = parser.parse_args()

	# Link introns to transcripts and vice versa, storing all possible introns and transcripts
	intron_to_transcripts = defaultdict(set)
	transcripts_to_introns = defaultdict(set)
	intron_set = set()
	tx_set = set()
	with open(os.path.abspath(args.introns_to_transcripts)) as f:
		f.readline()
		for line in f:
			tokens = line.strip().split('\t')
			intron_to_transcripts[tokens[3]].add(tokens[1])
			transcripts_to_introns[tokens[1]].add(tokens[3])
			intron_set.add(tokens[3])
			tx_set.add(tokens[1])

	# Get read counts for introns for each sample
	intron_count_dict = defaultdict(dict)
	tx_count_dict = defaultdict(dict)
	with open(os.path.abspath(args.read_counts)) as f:
		# Store all sample IDs from header
		samples = f.readline().strip().split('\t')[:-1]
		for line in f:
			tokens = line.strip().split('\t')
			intron = tokens[-1]
			if 'ENST' in intron:
				# Transcript entry
				tx = intron.split('|')[0]
				bp = intron.split('|')[6]
				if tx not in tx_set:
					print(tx)
				# Add exon read count info to dictionary
				for i in range(len(tokens)-1):
					tx_count_dict[samples[i]][tx] = (float(tokens[i]), float(bp))
			else:
				# Intron entry
				if intron not in intron_set:
					print(intron)
				# Add intron read count info to dictionary
				bases = intron.split(':')[1].split('-')
				bp = float(bases[1]) - float(bases[0]) + 1
				for i in range(len(tokens)-1):
					intron_count_dict[samples[i]][intron] = (float(tokens[i]), bp)

	outliers_to_samples = defaultdict(set)
	samples_to_outliers = defaultdict(set)
	# Identify outlier introns for each sample
	for sample in intron_count_dict:
		for transcript in transcripts_to_introns:
			if transcript in tx_count_dict[sample]:
				introns = list(transcripts_to_introns[transcript])
				intron_read_counts = [intron_count_dict[sample][x][0] for x in introns]
				intron_cutoff = (3*median_absolute_deviation(intron_read_counts))+median(intron_read_counts)
				tx_read_count = tx_count_dict[sample][transcript][0]/tx_count_dict[sample][transcript][1]
				for i in range(len(intron_read_counts)):
					av_count = intron_read_counts[i]/intron_count_dict[sample][introns[i]][1]
					if intron_read_counts[i] > intron_cutoff and av_count >= tx_read_count:
						outliers_to_samples[introns[i]].add(sample)
						samples_to_outliers[sample].add(introns[i])

	# Write output
	with open(os.path.abspath(args.output_file), 'w') as o:
		header = ['Intron'] + sorted(samples_to_outliers.keys())
		print('\t'.join(header), file=o)
		for outlier in outliers_to_samples:
			relevant_samples = outliers_to_samples[outlier]
			out_line = [outlier]
			for i in range(1, len(header)):
				if header[i] in relevant_samples:
					out_line.append('1')
				else:
					out_line.append('0')
			print('\t'.join(out_line), file=o)

