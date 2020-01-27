#!/usr/bin/env python3

"""
singleexpt_jxs_SRA_filter.py
Python 3 code for checking experimental junctions against previously collected
snaptron junctions from SRA experiments.

"""

import argparse
from datetime import datetime
import glob
import logging
import os
import pandas as pd
import sys


def snaptron_results_to_jxs(snaptron_result_lines, total_count=0,
                            min_sample_count=1, min_read_count=1):
    """Collects junctions found in select snaptron-queried SRA experiments.

    Input:
        snaptron_result_lines (opened file or incoming snaptron query results):
            line-by-line junction data
        total_count (int): a value by which to scale the sample count if passed
        min_sample_count (int): required minimum number of samples for a
            junction to be valid
        min_read_count (int): required minimum number of reads for a junction
            to be valid

    Parses incoming lines.

    Returns a dictionary with each key being one junction and its value being
    the number of samples in which the junction is found.
    """
    expt_jxs = {}
    for line in snaptron_result_lines:
        items = line.split('\t')
        if not items[0].endswith('I'):
            continue

        samp_count = int(items[13])
        if samp_count < min_sample_count:
            continue

        read_count = int(items[14])
        if read_count < min_read_count:
            continue

        chrom, left, right, strand = items[2], items[3], items[4], items[6]
        left = str(int(left) - 1)
        right = str(int(right) - 1)
        jx = ';'.join([chrom, left, right, strand])
        if total_count:
            expt_jxs[jx] = samp_count / total_count
        else:
            expt_jxs[jx] = samp_count

    return expt_jxs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Check junctions for SRA experiment evidence.'
    )
    parser.add_argument(
        '--snaptron-results',
        help='.txt file containing results from a previous snaptron search, '
             'or directory containing multiple of these.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output files: sampled spots, aligned junction '
             'reads, and SRA numbers with their p-values.'
    )
    parser.add_argument(
        '--log-level', default='INFO', help='INFO is the only level supported.'
    )
    parser.add_argument(
        '--junction-input', '-j',
        help='For querying external RNA-seq experiment(s) against the GTEX '
             'junction database: provide a .bed, *SJ.out.tab, or .csv file w/ '
             'junctions, or a directory containing multiple junction files.'
    )
    parser.add_argument(
        '--sra-filter', nargs='*',
        choices=[
            'embryo_stemcells', 'melanocyte_cellline',
            'lateembryo_tissue', 'embryo_primarycell',
            'lateembryo_primarycell', 'embryo_cellline',
            'lateembryo_cellline', 'embryo_tissue',
            'placenta_tissue', 'placenta_cellline',
            'aorta_cellline', 'aorta_tissue',
            'astrocyte_cellline', 'astrocyte_primarycell',
            'astrocyte_stemcells', 'biliarytree_all',
            'bone_cellline', 'bone_primarycell',
            'bone_tissue', 'ectoderm_cellline',
            'epithelialcell_cellline', 'epithelialcell_primarycell',
            'epithelialcell_stemcells', 'eye_cellline',
            'eye_primarycell', 'eye_tissue',
            'fallopiantube_cellline', 'fallopiantube_tissue',
            'fibroblast_cellline', 'fibroblast_primarycell',
            'fibroblast_stemcells', 'gallbladder_primarycell',
            'glialcell_cellline''gc_cl', 'glialcell_primarycell',
            'glialcell_stemcells', 'hematopoietic_cellline',
            'hematopoietic_primarycell', 'hematopoietic_stemcells',
            'hepatocyte_cellline', 'hepatocyte_primarycell',
            'inducedpluripotentstemcell_cellline',
            'isletoflangerhans_cellline',
            'isletoflangerhans_primarycell',
            'leukocyte_cellline', 'leukocyte_primarycell',
            'lymphocyte_cellline', 'lymphocyte_primarycell',
            'macrophage_cellline', 'macrophage_primarycell',
            'melanocyte_primarycell', 'melanocyte_cellline',
            'myeloidcell_cellline',
            'myeloidcell_primarycell', 'myoblast_cellline',
            'myoblast_primarycell', 'neonate_cellline',
            'neonate_primarycell', 'neonate_tissue',
            'oligodendrocyte_primarycell', 'oocyte_primarycell',
            'placenta_primarycell', 'platelet_cellline',
            'platelet_primarycell', 'pluripotentstemcell_cellline',
            'pluripotentstemcell_stemcells',
            'somaticstemcell_stemcells', 'thymus_primarycell',
            'thymus_tissue', 'zygote_primarycell',
            'mesenchymalstemcell_stemcells',
            'mesenchyme_primarycells', 'mesenchyme_stemcells'
        ]
    )
    parser.add_argument(
        '--batch-number', type=int,
        help='For downloading snaptron junctions for multiple cell types in '
             'parallel via slurm: batch numbers correspond to item number in '
             'sra can_file(s) above.'
    )

    args = parser.parse_args()
    ont_in = args.ontology_df
    out_path = args.output_path
    log_mode = args.log_level
    jx_input = args.junction_input
    batch_num = args.batch_number
    sra_sets = args.sra_filter

    now = str(datetime.now())
    now = now.split('.')[0].replace(' ', '_').replace(':', '.')

    csv_path = os.path.join(jx_input, '**/*.csv')
    file_list = glob.glob(csv_path, recursive=True)
    file_to_collect = file_list[batch_num]

    log_file = os.path.join(
        out_path, '{}_sra_filter_{}_log.txt'.format(batch_num, now)
    )
    logging.basicConfig(filename=log_file, level=log_mode)
    logging.info('input is: {}'.format(' '.join(sys.argv)))
    logging.info('filtering can_file {}'.format(file_to_collect))
    snaptron_jxs = {}
    all_snaptron_jxs = set()

    if not ont_in.endswith('.txt'):
        txt_path = os.path.join(ont_in, '*rawresults*.txt')
        ont_files = glob.glob(txt_path)
    else:
        ont_files = [ont_in]

    for ont in ont_files:
        name_tag = os.path.basename(ont).split('.')[0]
        name_tag = name_tag.split('_rawresults')[0]
        try:
            name_tag = name_tag.split('metaSRA-runs_')[1]
        except IndexError:
            pass
        if name_tag not in sra_sets:
            continue
        logging.info('collecting {} sra junctions'.format(name_tag))
        with open(ont) as lines:
            snaptron_jxs[name_tag] = snaptron_results_to_jxs(lines)

        all_snaptron_jxs = all_snaptron_jxs.union(set(snaptron_jxs[name_tag]))
        logging.info('all jxs updated: {} jxs'.format(len(all_snaptron_jxs)))

    logging.info('loading single expt:')
    logging.info('{}'.format(file_to_collect))
    single_jxs = pd.read_csv(file_to_collect)
    logging.info('df loaded; filtering on SRA:')
    single_jxs = single_jxs[~single_jxs['jx'].isin(all_snaptron_jxs)]
    logging.info('df filtered; writing to can_file:')

    new_filename = os.path.basename(file_to_collect).split('.')[0]
    new_filename += '_sra_filtered.csv'
    with open(os.path.join(out_path, new_filename), 'w') as output:
        single_jxs.to_csv(output, index=False)
