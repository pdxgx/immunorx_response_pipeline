#!/usr/bin/env python

import datetime
import argparse
from collections import defaultdict
from intervaltree import Interval, IntervalTree

def combineVariants(tumor_id, indel_size, *args):
	''' Parses a VCF to obtain variant entries and store as dictionary
		Keys in the dictionary are quadruples:
			(contig, position on contig, reference nucleotide, alternate nucleotide)
		Values are lists of the other data associated with that variant entry in the VCF,
			as well as the caller which produced the vcf
		If more than one caller called the variant, the associated data is added to the entry in the dictionary
		Unique header lines are also retained as a list for subsequent processing
		
		tumor_id: tumor identifier (string)
		indel_size: size cutoff for indels
		args: tuple of (caller, vcf)
			vcf: path vcf which has been normalized and decomposed using vt
			caller: name of the caller which produced the vcf
		
		Return values: variant dictionary, header info list
	'''
	# Establish dictionary to hold processed variants
	variants = defaultdict(set)
	headers = []
	# Parse each set
	for call_set in args:
		no_variants = False
		caller = call_set[0]
		vcf = call_set[1]
		with open(vcf, "r") as f:
			# Establish tumor sample as last column by default
			tumor_last = True
			# Filter calls to throw out - REJECT for MuTect, Tier 5 for MuSE
			reject_filters = ["REJECT", "Tier5", "DETP20", "IRC", "MMQS100", "MMQSD50", 
							  "MQD30", "MVC4", "MVF5", "NRC", "PB10", "RLD25", "SB1", "str10"]
			# Parse VCF
			for line in f:

				# Skip line if it's been established that there are not variants in VCF
				if no_variants:
					continue

				# Extract data from header lines
				if line[0:2] == "##":
					line = line.strip("\n")
					if "ID=" in line and line not in headers:
						headers.append(line)
			
				# Determine sample order from column headers
				elif line[0] == "#":
					tokens = line.strip("\n").split("\t")
					if len(tokens) < 10:
						no_variants = True
						continue
					if tokens[9] == tumor_id:
						tumor_last = False
					elif tokens[9] == "TUMOR" or tokens[9] == "PRIMARY":
						tumor_last = False
 					elif "NORM" in tokens[10] or "_N" in tokens[10]:
						tumor_last = False
			
				# Parse variant entries
				else:
					tokens = line.strip("\n").split("\t")
					filter_field = tokens[6]
					too_large = False
					# Check that variant passes filters before proceeding
					if len([item for item in reject_filters if (item in filter_field)]) == 0:
						contig = tokens[0]
						if contig.isdigit():
							contig = int(contig)
						position = int(tokens[1])
						id_field = tokens[2]
						ref = tokens[3]
						alt = tokens[4]
						# Eliminate large indels that may be structural variants
						if (len(ref) > len(alt) and len(ref) > indel_size) or (len(alt) > len(ref) and len(alt) > indel_size):
							too_large = True
						qual = tokens[5]
						info = tokens[7]
						if 'OLD_VARIANT' in info or 'OLD_CLUMPED' in info:
							new_info_list = []
							info_parts = info.split(';')
							for part in info_parts:
								if (part.startswith('OLD_VARIANT') or part.startswith('OLD_CLUMPED')) and len(part.split(':')[2]) > 1000:
									new_entry = ':'.join([':'.join(part.split(':')[0:2]), 'LONG_VARIANT'])
									new_info_list.append(new_entry)
								else:
									new_info_list.append(part)
							info = ';'.join(new_info_list)
						format_field = tokens[8]
						if tumor_last == True:
							normal = tokens[9]
							tumor = tokens[10]
						else:
							normal = tokens[10]
							tumor = tokens[9]
						entry = (contig, position, ref, alt)
						data = (id_field, qual, filter_field, info, format_field, normal, tumor, caller)
						if not too_large:
							variants[entry].add(data)
	
	return variants, headers


def parseHeaders(header_list):
	''' Parses a list of header lines to produce lists of format, filter, and info definitions
			
		header_list: list of header lines to process
		
		Return values: format, filter, and info lists
	'''
	# Establish output dictionaries
	format_list = []
	filter_list = []
	info_list = []
	
	# Iterate through header lines
	for item in header_list:
		
		# Store format header line	
		if 'FORMAT=<' in item:
			format_list.append(item)
		
		# Store filter header line	
		elif 'FILTER=<' in item:
			filter_list.append(item)
		
		# Store info header line		
		elif 'INFO=<' in item:
			info_list.append(item)
	
	return format_list, filter_list, info_list

def write_entry(variant_entry, variant_dict, filehandle):
	# Initialize lists to store each data type
	set_callers = []
	ID_list = []
	qual_list = []
	filter_list = []
	info_list = []
	tumor_format_dict = {'GT': set(), 'FREQ': set()}
	normal_format_dict = {'GT': set(), 'FREQ': set()}
	tumor_data = []
	normal_data = []
	
	# Loop through data from each caller that called the variant
	for callset in variant_dict[variant_entry]:
		# Store data on callers
		caller = callset[7].upper()
		set_callers.append(caller)
		
		# Store data on IDs
		if callset[0] != "." and callset[0] not in ID_list:
			ID_list.append(callset[0])
	
		# Store data on quality scores
		if callset[1] != "." and callset[1] != '0':
			qual_list.append(float(callset[1]))
	
		# Store data on info
		if callset[3] != "" and callset[3] != ".":
			info_items = callset[3].split(";")
			for item in info_items:
				if item not in info_list:
					info_list.append(item)
				
		# Adjust and store genotype fields
		format_field = callset[4].split(':')
		normal_entry = callset[5].split(':')
		tumor_entry = callset[6].split(':')
		for j in range(0, len(format_field)):
			# Genotype data
			if format_field[j] == 'GT':
				# Normal entry - adjust 0 or 1 to 0/0 or 1/1
				if (normal_entry[j] != '.' and normal_entry[j] != '0' and normal_entry[j] != '1'):
					normal_format_dict['GT'].add(normal_entry[j])
				elif normal_entry[j] == '0':
					normal_format_dict['GT'].add('0/0')
				elif normal_entry[j] == '1':
					normal_format_dict['GT'].add('1/1')
				# Tumor entry - adjust 0 or 1 to 0/0 or 1/1
				if (tumor_entry[j] != '.' and tumor_entry[j] != '0' and tumor_entry[j] != '1'):
					tumor_format_dict['GT'].add(tumor_entry[j])
				elif tumor_entry[j] == '0':
					tumor_format_dict['GT'].add('0/0')
				elif tumor_entry[j] == '1':
					tumor_format_dict['GT'].add('1/1')
			# VAF data
			elif format_field[j] == 'FREQ':
				if normal_entry[j] != '.':
					normal_format_dict['FREQ'].add(float(normal_entry[j].replace('%', '')))
				if tumor_entry[j] != '.':
					tumor_format_dict['FREQ'].add(float(tumor_entry[j].replace('%', '')))
		if len(normal_format_dict['GT']) > 0:
			if '0/1' in normal_format_dict['GT']:
				norm_gt = '0/1'
			elif '1/1' in normal_format_dict['GT']:
				norm_gt = '1/1'
			else:
				norm_gt = '0/0'
			if '0/1' in tumor_format_dict['GT']:
				tum_gt = '0/1'
			elif '1/1' in tumor_format_dict['GT']:
				tum_gt = '1/1'
			else:
				tum_gt = '0/0'
		else:
			norm_gt = '0/0'
			tum_gt = '0/0'
		if len(tumor_format_dict['FREQ']) > 0:
			norm_freq = str(sum(normal_format_dict['FREQ'])/len(normal_format_dict['FREQ'])) + '%'
			tum_freq = str(sum(tumor_format_dict['FREQ'])/len(tumor_format_dict['FREQ'])) + '%'
		else:
			norm_freq = '.'
			tum_freq = '.'

	# Finalize ID data
	if len(ID_list) == 0:
		ID = "."
	else:
		ID = ";".join(ID_list)
	
	# Finalize quality score data
	if len(qual_list) != 0:
		QUAL = sum(qual_list)/len(qual_list)
	else:
		QUAL = "."

	# Finalize info data
	set_callers = list(set(set_callers))
	if len(info_list) != 0:
		INFO = ";".join(info_list) + ";" + ";".join(set_callers)
	else:
		INFO = ";".join(set_callers)
	
	# Finalize genotype fields
	FORMAT = "GT:FREQ"
	NORMAL = ":".join([norm_gt, norm_freq])
	TUMOR = ":".join([tum_gt, tum_freq])
	
	# Write data to VCF
	outline = "\t".join([str(variant_entry[0]), str(variant_entry[1]), ID, variant_entry[2], variant_entry[3], 
						 str(QUAL), "PASS", INFO, FORMAT, TUMOR, NORMAL]) + "\n"
	filehandle.write(outline)

def writeVCF(number, variants, outdir, sample, callers, format, filter, info):
	''' Write a consensus VCF from a refined set of variants
		Also produces output file with descriptive data about the call sets
			Number of callers used and which callers were used
		
		variants = dictionary of filtered variants
		outdiir = path to output directory
		sample = sample name to attach to output files
		callers = list of callers used to generate call sets
		format = header format definition dictionary
		filter = header filter definition dictionary
		info = header info definition dictionary
		
		Return value: none
	'''
	# Establish output files
	outputVCF = os.path.join('.'.join([sample, "consensus.vcf"]))

	# Establish header info
	today = datetime.date.today()
	caller_headers = []
	for caller in callers:
		line = '##INFO=<ID=' + caller.upper() + ',Number=.,Type=String,Description="Variant called by ' + caller.lower() + '">\n'
		caller_headers.append(line)
	
	with open(outputVCF, "w") as f:
		# Write header lines
		f.write("##fileformat=VCFv4.2\n")
		f.write("##fileDate=" + str(today) + "\n")
		# Write format header lines
		for ID in format:
			f.write(ID+'\n')
		f.write('##FORMAT=<ID=CALLER,Number=.,Type=String,Description="Caller used to produce data"\n')
		# Write info header lines
		for ID in info:
			f.write(ID+'\n')
		# Write filter header lines
		for ID in filter:
			f.write(ID+'\n')
		for header in caller_headers:
			f.write(header)
		# Write column header lines
		f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n")
		
		# Loop through each variant to collect data about pindel overlap
		pindel_intervals = defaultdict(IntervalTree)
		for entry in variants.keys():
			for callset in variants[entry]:
				caller = callset[7].upper()
				if caller == "PINDEL":
					# Pindel variant - add to interval tree and set of pindel variants
					chrom = entry[0]
					if len(entry[2]) > len(entry[3]):
						indel_length = len(entry[2])
					else:
						indel_length = len(entry[3])
					pindel_intervals[chrom][int(entry[1]):int(entry[1])+indel_length] = entry

		# Loop through variants and write to VCF if relevant
		for entry in sorted(variants.keys()):
			
			# Collect data on which caller(s) called the variant and store data on pindel structural events
			caller_list = []
			for callset in variants[entry]:
				caller = callset[7].upper()
				caller_list.append(caller)

			caller_list.sort()
			callers_used = ",".join(caller_list)
			
			# Check whether non-pindel variants are overlapped by pindel
			pindel_overlap = False
			if "PINDEL" not in caller_list:
				overlap = pindel_intervals[entry[0]][int(entry[1]):int(entry[1])+max(len(entry[2]), len(entry[3]))]
				if len(overlap) > 0:
					pindel_overlap = True

			# If at least the required number of callers called the variant, write to output VCF
			if len(caller_list) >= number and not pindel_overlap:
				write_entry(entry, variants, f)
			elif caller_list == ["PINDEL"]:
				write_entry(entry, variants, f)

			
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcfs', type=str, required=True,
            help='comma-separated list of paths to input vcf'
        )
    parser.add_argument('-c', '--callers', type=str, required=True,
            help='comma-separated list of callers which generated the input vcfs (in same order)'
        )
    parser.add_argument('-o', '--outputdir', type=str, required=True,
            help='path to output consensus vcf'
        )
    parser.add_argument('-n', '--number', type=int, required=False,
            help='number of callers required to maintain a variant'
        )
    parser.add_argument('-s', '--sample', type=str, required=False,
            help='unique sample name for output'
        )
    parser.add_argument('-t', '--tumor-id', type=str, required=False,
            help='tumor ID in mutect VCF'
        )
    parser.add_argument('-i', '--indel-size', type=int, required=False,
            help='tumor ID in mutect VCF'
        )
    args = parser.parse_args()
    
    # Turn comma separated VCF and caller lists into list objects
    vcf_list = args.vcfs.split(",")
    caller_list = args.callers.split(",")
    
    # Turn VCF and callers lists into a list of tuples for combineVariants()
    set_list = []
    for i in range(0, len(caller_list)):
    	this_set = (caller_list[i], vcf_list[i])
    	set_list.append(this_set)
    
    # Produce a variant dictionary and list of header data	
    variant_dict, header_data = combineVariants(args.tumor_id, args.indel_size, *set_list)
    
    # Produce a dictionary of unique header data
    format_dict, filter_dict, info_dict = parseHeaders(header_data)
    
    # Write consensus VCF and associated metadata
    writeVCF(args.number, variant_dict, args.outputdir, args.sample, caller_list, format_dict, filter_dict, info_dict)
