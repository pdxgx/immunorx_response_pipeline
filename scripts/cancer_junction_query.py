#!/usr/bin/env python3

"""
_______.py
Python code for querying single sample RNA-seq results against TCGA/GTEx
junction database to identify cancer-specific junctions.

"""

import argparse
from collections import defaultdict
import csv
from datetime import datetime
import glob
from intervaltree import IntervalTree
import logging
import os
import pandas as pd
import re
import sqlite3 as sql
import sys


_ID = '_IDs'
_COV = '_Covs'
_PER = '_Sample_Percents'
_MED_COV = '_Median_Coverage'
_STG = '_Stages'
_JX_SAMP_TABLE = 'jx_sample_map'
_JX_ANN_TABLE = 'jx_annotation_map'
_PHEN_TABLE = 'sample_phenotype_map'

_TCGA_CANCER_TYPES = [
    'Acute_Myeloid_Leukemia', 'Adrenocortical_Carcinoma',
    'Bladder_Urothelial_Carcinoma', 'Brain_Lower_Grade_Glioma',
    'Breast_Invasive_Carcinoma',
    'Cervical_Squamous_Cell_Carcinoma_and_Endocervical_Adenocarcinoma',
    'Cholangiocarcinoma', 'Colon_Adenocarcinoma', 'Esophageal_Carcinoma',
    'Glioblastoma_Multiforme', 'Head_and_Neck_Squamous_Cell_Carcinoma',
    'Kidney_Chromophobe', 'Kidney_Renal_Clear_Cell_Carcinoma',
    'Kidney_Renal_Papillary_Cell_Carcinoma', 'Liver_Hepatocellular_Carcinoma',
    'Lung_Adenocarcinoma', 'Lung_Squamous_Cell_Carcinoma',
    'Lymphoid_Neoplasm_Diffuse_Large_B_cell_Lymphoma', 'Mesothelioma',
    'Ovarian_Serous_Cystadenocarcinoma', 'Pancreatic_Adenocarcinoma',
    'Pheochromocytoma_and_Paraganglioma', 'Prostate_Adenocarcinoma',
    'Rectum_Adenocarcinoma', 'Sarcoma', 'Skin_Cutaneous_Melanoma',
    'Stomach_Adenocarcinoma', 'Testicular_Germ_Cell_Tumors', 'Thymoma',
    'Thyroid_Carcinoma', 'Uterine_Carcinosarcoma',
    'Uterine_Corpus_Endometrial_Carcinoma', 'Uveal_Melanoma'
]

_CANCER_TYPES_PRIMARY = [
    'Acute_Myeloid_Leukemia', 'Adrenocortical_Carcinoma', 'Astrocytoma',
    'Bladder_Urothelial_Carcinoma', 'Breast_Invasive_Carcinoma',
    'Cervical_Adenosquamous', 'Cervical_Squamous_Cell_Carcinoma',
    'Cholangiocarcinoma', 'Colon_Adenocarcinoma',
    'Dedifferentiated_Liposarcoma', 'Desmoid_Tumor',
    'Endocervical_Adenocarcinoma', 'Esophagus_Adenocarcinoma',
    'Esophagus_Squamous_Cell_Carcinoma', 'Glioblastoma_Multiforme',
    'Head_and_Neck_Squamous_Cell_Carcinoma', 'Kidney_Chromophobe',
    'Kidney_Renal_Clear_Cell_Carcinoma',
    'Kidney_Renal_Papillary_Cell_Carcinoma', 'Leiomyosarcoma',
    'Liver_Hepatocellular_Carcinoma', 'Lung_Adenocarcinoma',
    'Lung_Squamous_Cell_Carcinoma',
    'Lymphoid_Neoplasm_Diffuse_Large_B_cell_Lymphoma',
    'Malignant_Peripheral_Nerve_Sheath_Tumors', 'Mesothelioma',
    'Myxofibrosarcoma', 'Oligoastrocytoma', 'Oligodendroglioma',
    'Ovarian_Serous_Cystadenocarcinoma', 'Pancreatic_Adenocarcinoma',
    'Paraganglioma', 'Pheochromocytoma', 'Prostate_Adenocarcinoma',
    'Rectum_Adenocarcinoma', 'Skin_Cutaneous_Melanoma',
    'Stomach_Adenocarcinoma', 'Synovial_Sarcoma',
    'Testicular_Germ_Cell_Tumors', 'Thymoma', 'Thyroid_Carcinoma',
    'Undifferentiated_Pleomorphic_Sarcoma', 'Uterine_Carcinosarcoma',
    'Uterine_Corpus_Endometrial_Carcinoma', 'Uveal_Melanoma'
    ]

_ALL_CANCERS = list(
    set(_CANCER_TYPES_PRIMARY).union(set(_TCGA_CANCER_TYPES))
)


def gtf_to_cds(gtf_file):
    """ References cds_dict to get cds bounds for later Bowtie query
    Keys in the dictionary are transcript IDs, while entries are lists of
        relevant CDS/stop codon data
        Data: [chromosome, left, right, +/- strand]
    Writes cds_dict as a pickled dictionary
    gtf_file: input gtf file to process
    dictdir: path to directory to store pickled dicts

    NOTE: from neoepiscope transcript.py

    Return value: dictionary
    """
    cds_dict = defaultdict(list)
    # Parse GTF to obtain CDS/stop codon info
    with open(gtf_file) as f:
        for line in f:
            if line[0] != '#':
                tokens = line.strip().split('\t')
                if tokens[2] == 'gene' and 'protein_coding' in line:
                    gene_id = re.sub(r'.*gene_id \"([A-Z0-9._]+)\"[;].*',
                                     r'\1', tokens[8])
                    gene_type = re.sub(r'.*gene_type \"([a-z_]+)\"[;].*',
                                       r'\1', tokens[8])
                    if gene_type == 'protein_coding':
                        # Create new dictionary entry for new each gene
                        cds_dict[gene_id].append(
                            [tokens[0], int(tokens[3]), int(tokens[4]),
                             tokens[6]]
                        )
    return cds_dict


def cds_to_tree(cds_dict):
    """ Creates searchable tree of chromosome intervals from CDS dictionary

    Each chromosome is stored in the dictionary as an interval tree object
        Intervals are added for each CDS, with the associated transcript ID
        Assumes transcript is all on one chromosome - does not work for
            gene fusions
    Writes the searchable tree as a pickled dictionary
    cds_dict: CDS dictionary produced by gtf_to_cds()

    NOTE: from neoepiscope transcript.py

    Return value: searchable tree
    """
    searchable_tree = {}
    # Add genomic intervals to the tree for each transcript
    for gene_id in cds_dict:
        gene = cds_dict[gene_id]
        chrom = gene[0][0]
        # Add new entry for chromosome if not already encountered
        if chrom not in searchable_tree:
            searchable_tree[chrom] = {}
        # Add CDS interval to tree with transcript ID
        for cds in gene:
            left = cds[1]
            right = cds[2]
            strand = cds[3]
            # Interval coordinates are inclusive of left, exclusive of right
            if strand not in searchable_tree[chrom]:
                searchable_tree[chrom][strand] = IntervalTree()
            if right > left:
                searchable_tree[chrom][strand][left:right+1] = gene_id
    return searchable_tree


def jx_ends_in_cds(junction, cds_tree, id_name_dict):
    """Check found junctions for coding region overlap

    Input:
        junction information: chromosome, left and right boundaries of the
            junction, and strand
        CDS tree containing coding regions, created by gtf_to_cds
        dictionary mapping gene ids from .gtf file to gene names

    Checks to see whether either junction end overlaps coding regions.  If
    either or both do, collects gene ids and names for the overlapping genes.
    If both sides overlap, checks to see if any of the genes are the same.

    Returns eight entries comprising new column information for the junction
    database; columns contain the following information:
    1) binary: whether both ends of the junction overlap no known gene regions
    2) binary: whether one end of the junction only overlaps gene regions
    3) binary: whether two ends of the junction overlap different genes
    4) binary: whether two ends of the junction overlap the same gene
    5) comma-separated list of gene ids overlapped by the 5' junction end
    6) comma-separated list of gene ids overlapped by the 3' junction end
    7) comma-separated list of gene names overlapped by the 5' junction end
    8) comma-separated list of gene names overlapped by the 3' junction end
    """
    no_overlap = 0
    one_overlap = 0
    both_same = 0
    both_diff = 0
    left_genes = []
    right_genes = []
    left_names = []
    right_names = []
    chrom, left, right, strand = junction.split(';')
    left = int(left)
    right = int(right)
    try:
        jx_start = cds_tree[chrom][strand].overlaps(left)
        jx_stop = cds_tree[chrom][strand].overlaps(right)
    except KeyError:
        return (no_overlap, one_overlap, both_diff, both_same, left_genes,
                right_genes, left_names, right_names)

    if jx_start or jx_stop:
        for start_set in list(cds_tree[chrom][strand][left]):
            if start_set[2] not in left_genes:
                left_genes.append(start_set[2])
            name = id_name_dict[start_set[2]]
            if name not in left_names:
                left_names.append(name)
        for stop_set in list(cds_tree[chrom][strand][right]):
            if stop_set[2] not in right_genes:
                right_genes.append(stop_set[2])
            name = id_name_dict[stop_set[2]]
            if name not in right_names:
                right_names.append(name)
        if jx_start and jx_stop:
            num_same_genes = len(set(left_genes) & set(right_genes))
            if num_same_genes > 0:
                both_same = 1
            if ((len(right_genes) - num_same_genes > 0) or
                    (len(left_genes) - num_same_genes > 0)):
                both_diff = 1
        else:
            one_overlap = 1

    if strand == '+':
        fivepr_genes = ','.join(left_genes)
        threepr_genes = ','.join(right_genes)
        fivepr_names = ','.join(left_names)
        threepr_names = ','.join(right_names)
    else:
        fivepr_genes = ','.join(right_genes)
        threepr_genes = ','.join(left_genes)
        fivepr_names = ','.join(right_names)
        threepr_names = ','.join(left_names)

    return (no_overlap, one_overlap, both_diff, both_same, fivepr_genes,
            threepr_genes, fivepr_names, threepr_names)


def check_annotations(junction, junction_dict):
    """Adds annotated splice junctions from .gtf file to the junction list.

    Junction column key:
        0 = neither junction side is annotated
        1 = one junction side is annotated
        2 = both junction sides are annotated, but not together
        3 = junction is fully annotated

    NOTE: right end of junction is 2 locations lower than what is annotated in
          GENCODE. To accommodate this, all junction right ends have 2 added
          before being checked against the junction annotation dictionary.

    """
    tokens = junction.split(';')
    chrom, left, right, strand = tokens[0], tokens[1], tokens[2], tokens[3]
    right = str(int(right) + 2)
    junction = chrom + ';' + left + ';' + right + ';' + strand
    try:
        if junction in junction_dict[chrom][strand]['full']:
            annotated_col = 3
        else:
            annotated_col = 0
            if strand == '+':
                five_site = left
                three_site = right
            else:
                five_site = right
                three_site = left
            if five_site in junction_dict[chrom][strand]['fivepr']:
                annotated_col += 1
            if three_site in junction_dict[chrom][strand]['threepr']:
                annotated_col += 1
    except KeyError:
        annotated_col = 0
    return annotated_col


def extract_splice_sites(gtf_file):
    """Extracts splice site anns_same_same from .gtf file

    This function is a modified version of one that is part of HISAT.
    Copyright 2014, Daehwan Kim <infphilo@gmail.com>

    HISAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HISAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
    """
    genes = defaultdict(list)
    trans = {}

    annotations = {}

    with open(gtf_file) as gtf:
        # Parse valid exon lines from the annotation file into a dict by
        # transcript_id
        for line in gtf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()

            try:
                item = line.split('\t')
                chrom, left, right, strand = item[0], item[3], item[4], item[6]
                feature, values = item[2], item[8]
            except ValueError:
                continue
            left, right = int(left), int(right)

            if feature != 'exon' or left >= right:
                continue

            values_dict = {}
            for attr in values.split(';')[:-1]:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

            if 'gene_id' not in values_dict or \
                    'transcript_id' not in values_dict:
                continue

            transcript_id = values_dict['transcript_id']
            if transcript_id not in trans:
                trans[transcript_id] = [chrom, strand, [[left, right]]]
                genes[values_dict['gene_id']].append(transcript_id)
            else:
                trans[transcript_id][2].append([left, right])

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]

    # Calculate and print the unique junctions
    junctions = set()
    for chrom, strand, exons in trans.values():
        for i in range(1, len(exons)):
            junctions.add((chrom, exons[i - 1][1], exons[i][0], strand))
    junctions = sorted(junctions)

    for chrom, left, right, strand in junctions:
        if chrom not in annotations:
            annotations[chrom] = {}
            annotations[chrom]['+'] = {'full': [], 'fivepr': [], 'threepr': []}
            annotations[chrom]['-'] = {'full': [], 'fivepr': [], 'threepr': []}
        annotations[chrom][strand]['full'].append(
            ';'.join([chrom, str(left), str(right), strand]))
        if strand == '+':
            five_site = str(left)
            three_site = str(right)
        else:
            five_site = str(right)
            three_site = str(left)
        if five_site not in annotations[chrom][strand]['fivepr']:
            annotations[chrom][strand]['fivepr'].append(five_site)
        if three_site not in annotations[chrom][strand]['threepr']:
            annotations[chrom][strand]['threepr'].append(three_site)
    return annotations


def make_id_name_dict(gtf_file):
    """Creates dictionary mapping gene IDs to gene names

    Input: gtf_file containing gene names and gene IDs

    Parses each line of the gtf file, and adds gene name-ID pairs to the dict.

    Returns the created dictionary.
    """
    id_name_dict = {}
    with open(gtf_file) as gtf:
        for line in gtf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()
            try:
                values = line.split('\t')[-1]
            except ValueError:
                continue
            values_dict = {}
            for attr in values.split(';')[:-1]:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')
            if 'gene_id' not in values_dict or 'gene_name' not in values_dict:
                continue
            if values_dict['gene_id'] not in id_name_dict:
                id_name_dict[values_dict['gene_id']] = values_dict['gene_name']
    return id_name_dict


def sj_out_to_jxs(sj_out):
    """Accepts STAR-generated SJ.out.tab file & returns list of its junctions.

    Input: sj_out (str) string pointing to STAR output __SJ.out file containing
        junction calls from a STAR alignment.

    Returns a list of the file's unique junctions in 0-based closed coordinates
    """
    all_junctions = []
    strand_mapper = ['?', '+', '-']
    with open(sj_out) as sj:
        jx_file = csv.reader(sj, delimiter='\t')
        for line in jx_file:
            chrom, left, right, strand = line[0], line[1], line[2], line[3]
            strand = strand_mapper[int(strand)]
            if strand == '?':
                continue
            if 'chr' not in chrom:
                chrom = 'chr' + chrom
            left = str(int(left) - 1)
            right = str(int(right) - 1)
            jx = ';'.join([chrom, left, right, strand])
            all_junctions.append(jx)
    return list(set(all_junctions))


def jx_input_to_files(jx_input, recursive_glob=False):
    """Parses jx file input and returns list of actual junction files.

    Input: jx_input (str) a file containing junctions, or a directory with
        several of these, or with several subdirectories with these.

    Parses the generic junction input string and determines actual junction
    file names and locations.

    Returns a list of junction files.
    """
    file_list = []
    if not (jx_input.endswith('SJ.out.tab')):
        wd = os.path.abspath(jx_input)
        if recursive_glob:
            logging.info('using recursive glob...')
            logging.info('base path is {}'.format(jx_input))
            SJ_path = os.path.join(jx_input, '**/*SJ.out.tab')
            file_list.extend(glob.glob(SJ_path, recursive=True))
        else:
            SJ_path = os.path.join(wd, jx_input, '*SJ.out.tab')
            file_list.extend(glob.glob(SJ_path))
    else:
        file_list.append(jx_input)
    return file_list


def jx_input_to_jx_list(jx_input, separate=False):
    """
    Accepts command line file input; parses file(s) and returns list of jxs.

    :param jx_file:
    :return:
    """
    jx_list = []
    files = jx_input_to_files(jx_input)
    for file in files:
        if file.endswith('SJ.out.tab'):
            jxs = sj_out_to_jxs(file)
        if separate:
            jx_list.append(jxs)
        else:
            jx_list.extend(jxs)

    if jx_list == []:
        print('jx input is:', jx_input)
        print('files are:', files)
        print('junction input must be in SJ.out.tab format.')
        exit()
    return jx_list


def query_single_expt(jx_input, gtf_file, out_path, now, db_conn, prev_cans=[],
                      recursive=False):
    """Extracts provided junctions and annotates them with TCGA/GTEx info.

    Input:
        jx_input (str):
        gtf_file (str):
        out_path (str): directory where output should be stored.
        now (str): timestamp for labeling output
        db_conn (sqlite3 database connection): an open connection to the db
        prev_cans (list): list of TCGA cancer types to collect junction
            prevalences from
        recursive (bool): whether or not to look recursively for junction files
            inside the jx_input directory.

    Determines a list of single-sample jx files. Extracts junctions and cohort
    prevalences for select (prev_cans) TCGA cancer types.  Extracts normal
    junctions from GTEx and TCGA. For each set of single-sample junctions,
    retains only those that are not in TCGA/GTEx normal samples. Adds TCGA
    cohort prevalences to sample-specific, cancer-specific junctions. Adds
    annotation values (fully, partially, or not GENCODE-annotated).

    Returns None
    """
    logging.info('collecing junctions for single experiment files...')
    all_jx_sets = {}

    files = jx_input_to_files(jx_input, recursive)
    for file in files:
        logging.info('collecting junctions for {}...'.format(file))
        name_tag = os.path.basename(file).split('.')[0]
        all_jx_sets[name_tag] = jx_input_to_jx_list(file)

    logging.info('{} single experiments to compare:'.format(len(all_jx_sets)))
    for basename, junctions in all_jx_sets.items():
        logging.info('{}'.format(basename))

    if len(all_jx_sets) == 0:
        logging.info('no junction sets to compare: exiting')
        print('no junction sets to compare: exiting')
        exit()

    if gtf_file:
        logging.info('starting gtf file parsing...')
        annotations = extract_splice_sites(gtf_file)
        logging.info('gtf parsing complete.')

    logging.info('collecting normal junctions...')
    norm_command = (
        'SELECT jx FROM (SELECT DISTINCT jx_id id FROM {js} '
        'INNER JOIN {sp} ON {js}.recount_id == {sp}.recount_id '
        'WHERE tumor_normal == 1) INNER JOIN {ja} ON id == {ja}.jx_id;'
        ''.format(js=_JX_SAMP_TABLE, sp=_PHEN_TABLE, ja=_JX_ANN_TABLE)
    )
    norm_jxs = pd.read_sql_query(norm_command, db_conn)['jx'].tolist()
    logging.info('norm junction collection complete.')

    logging.info('collecting tcga junctions...')
    tcga_command = (
        'SELECT jx, annotation FROM '
        '(SELECT DISTINCT jx_id id FROM {js} '
        'INNER JOIN {sp} ON {js}.recount_id == {sp}.recount_id '
        'WHERE tumor_normal == 0) INNER JOIN {ja} ON id == {ja}.jx_id;'
        ''.format(js=_JX_SAMP_TABLE, sp=_PHEN_TABLE, ja=_JX_ANN_TABLE)
    )
    tcga_df = pd.read_sql_query(tcga_command, db_conn)
    tcga_jxs = tcga_df['jx'].tolist()
    logging.info('cancer junction collection complete.')

    full_df = pd.DataFrame({'jx': [], 'annotation': []})
    for cancer in prev_cans:
        logging.info('starting prevalence collection for {}...'.format(cancer))
        if cancer in _TCGA_CANCER_TYPES:
            cancer_col_name = 'project_type_label'
        else:
            cancer_col_name = 'primary_type'

        per_col = cancer + _PER

        sample_count_command = (
            'SELECT COUNT (*) FROM {} WHERE tumor_normal == 0 '
            'AND {} == "{}";'.format(_PHEN_TABLE, cancer_col_name, cancer)
        )
        count = pd.read_sql_query(sample_count_command, db_conn)
        count = count['COUNT (*)'][0]

        select_command = (
            'SELECT {ja}.jx, {ja}.annotation, COUNT (phen_recount) '
            'FROM (SELECT phen_recount, nonnorm_jxs '
            'FROM (SELECT {js}.jx_id nonnorm_jxs, phen_recount '
            'FROM (SELECT recount_id phen_recount '
            'FROM {sp} '
            'WHERE {sp}.{can_col} == "{can}" AND {sp}.tumor_normal == 0) '
            'INNER JOIN {js} ON phen_recount == {js}.recount_id) '
            'LEFT OUTER JOIN '
            '(SELECT DISTINCT jx_id nor_id FROM {js} '
            'INNER JOIN {sp} ON {js}.recount_id == {sp}.recount_id '
            'WHERE tumor_normal == 1)'
            'ON nonnorm_jxs == nor_id '
            'WHERE nor_id IS NULL) '
            'INNER JOIN {ja} ON {ja}.jx_id==nonnorm_jxs '
            'GROUP BY ({ja}.jx_id);'.format(
                js=_JX_SAMP_TABLE, ja=_JX_ANN_TABLE, sp=_PHEN_TABLE,
                can=cancer, can_col=cancer_col_name
            )
        )
        query_result = pd.read_sql_query(select_command, db_conn)
        col_rename = {'COUNT (phen_recount)': per_col}
        query_result.rename(columns=col_rename, inplace=True)
        query_result[per_col] = query_result[per_col] / count
        query_result = query_result.sort_values(by=[per_col], ascending=False)

        full_df = pd.merge(
            full_df, query_result, on=['jx', 'annotation'], how='outer'
        ).fillna(0)
        logging.info('prevalences for {} complete.'.format(cancer))
        logging.info('dataframe length is {}.'.format(len(full_df)))

    for name_tag, junctions in all_jx_sets.items():
        logging.info('starting junction collection for {}...'.format(name_tag))
        spf_jxs = list(set(junctions) - set(norm_jxs))
        tcga_all_overlaps = set(spf_jxs).intersection(set(tcga_jxs))
        tcga_prev_jxs = set(full_df['jx'].tolist()).intersection(set(spf_jxs))
        unique_jxs = list(set(spf_jxs) - tcga_all_overlaps)

        logging.info('For junction set {}'.format(name_tag))
        logging.info('total number of junctions is: {}'.format(len(junctions)))
        logging.info(
            'total number not found in GTEx is: {}'.format(len(spf_jxs))
        )
        logging.info(
            'number of unique junctions is: {}'.format(len(unique_jxs))
        )
        logging.info(
            'number of junctions found in TCGA is: {}'
            ''.format(len(tcga_all_overlaps))
        )
        logging.info(
            'number of junctions found in TCGA cancers of interest is {}'
            ''.format(len(tcga_prev_jxs))
        )

        logging.info('merging dataframes...')
        in_prev_overlaps = list(
            set(spf_jxs).intersection(set(full_df['jx'].tolist()))
        )
        prev_df = full_df[full_df.jx.isin(in_prev_overlaps)]
        logging.info('length of new prevalence df is {}'.format(len(prev_df)))
        tcga_not_prevs = list(set(spf_jxs) - set(in_prev_overlaps))

        overlap_df = pd.DataFrame({'jx': tcga_not_prevs})
        if gtf_file:
            overlap_df['annotation'] = overlap_df.jx.apply(
                lambda x: check_annotations(x, annotations)
            )
        else:
            overlap_df['annotation'] = ''

        logging.info('length of overlap df is {}'.format(len(overlap_df)))
        overlap_df = pd.concat([prev_df, overlap_df], sort=False)
        logging.info('after merging, length is {}'.format(len(overlap_df)))

        out_file = os.path.join(
            out_path, '{}_queryresults_{}.csv'.format(name_tag, now)
        )
        logging.info('writing output to {}'.format(out_file))
        print('writing output to {}'.format(out_file))
        with open(out_file, 'w') as output:
            overlap_df.to_csv(output, index=False)

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Query RNA-seq expts against GTEx/TCGA junction database.'
    )

    parser.add_argument(
        '--db-path', '-d', default='./',
        help='give the path for storing the created sql database.'
    )
    parser.add_argument(
        '--log-level', default='INFO', help='INFO is the only level supported.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give the path to store junction scoring output.'
    )
    parser.add_argument(
        '--junction-input', '-j',
        help='For querying external RNA-seq experiment(s) against the GTEX '
             'junction database: provide a .bed, *SJ.out.tab, or .csv file w/ '
             'junctions, or a directory containing multiple junction files.'
    )
    parser.add_argument(
        '--gtf-file', '-g',
        help='If "annotate" is selected: gtf file containing CDS annotation '
             'is also required to determine whether junctions occur in '
             'protein coding regions of a gene.'
    )
    parser.add_argument(
        '--tumor-prevalences', nargs='*', default=['Blood'],
        choices=_ALL_CANCERS,
        help='These are the cancers for which TCGA prevalence values will be '
             'given for junctions in a single experiment query.'
    )
    parser.add_argument(
        '--recursive-glob', action='store_true',
        help='for single sample queries, use this option to collect junction '
             'files recursively in the junction-input directory provided.'
    )

    args = parser.parse_args()
    db_path = args.db_path
    log_mode = args.log_level
    out_path = args.output_path
    jx_input = args.junction_input
    gtf_path = args.gtf_file
    prevs_to_print = args.tumor_prevalences
    rec_glob = args.recursive_glob

    try:
        db_name = os.path.join(db_path, 'new_jx_index.db')
        conn = sql.connect(db_name)
        index_db = conn.cursor()
    except sql.OperationalError:
        print('If OperationalError is "unable to open database file": ')
        print('make sure -d gives the PATH to the database directory, ')
        print('not the database itself.')
        raise sql.OperationalError

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'query_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))

    query_single_expt(
        jx_input, gtf_path, out_path, now, conn, prev_cans=prevs_to_print,
        recursive=rec_glob
    )
