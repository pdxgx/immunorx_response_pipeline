#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import pickle
import subprocess
import tempfile
import warnings
from neoepiscope import paths

def get_affinity_netMHCpan(
    peptides, allele, netmhcpan, version, scores, dic, remove_files=True
):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: allele to use for binding affinity
                    (string, format HLA-A02:01)
        netmhcpan: path to netMHCpan executable
        version: version of netMHCpan software
        scores: list of scoring methods
        dic: path to available alleles dictionary
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities
                        as strings)
    """
    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(dic, "rb") as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        allele = allele.replace("*", "")
        if allele not in avail_alleles["".join(["netMHCpan", str(version)])]:
            warnings.warn(
                " ".join([allele, "is not a valid allele for netMHCpan"]), Warning
            )
            score_form = tuple(["NA" for i in range(0, len(scores))])
            return [(peptides[i],) + score_form for i in range(0, len(peptides))]
        # Establish return list and sample id
        sample_id = ".".join(
            [peptides[0], str(len(peptides)), allele, "netmhcpan", version]
        )
        affinities = []
        # Write one peptide per line to a temporary file for input
        peptide_file = tempfile.mkstemp(
            suffix=".peptides", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(peptide_file)
        with open(peptide_file, "w") as f:
            for sequence in peptides:
                print(sequence, file=f)
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(
            suffix=".netMHCpan.out", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(mhc_out)
        err_file = tempfile.mkstemp(
            suffix=".netMHCpan.err", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(err_file)
        with open(err_file, "w") as e:
            # Run netMHCpan
            if version == "3":
                subprocess.check_call(
                    [
                        netmhcpan,
                        "-a",
                        allele,
                        "-inptype",
                        "1",
                        "-p",
                        "-xls",
                        "-xlsfile",
                        mhc_out,
                        peptide_file,
                    ],
                    stderr=e,
                )
            elif version == "4":
                subprocess.check_call(
                    [
                        netmhcpan,
                        "-BA",
                        "-a",
                        allele,
                        "-inptype",
                        "1",
                        "-p",
                        "-xls",
                        "-xlsfile",
                        mhc_out,
                        peptide_file,
                    ],
                    stderr=e,
                )
        with open(mhc_out, "r") as f:
            f.readline()
            f.readline()
            for i in range(0, len(peptides)):
                tokens = f.readline().strip("\n").split("\t")
                # for v3, tokens[5] is affinity, tokens[6] is rank
                # for v4, tokens[6] is affinity, tokens[7] is rank
                if version == "3":
                    result_dict = {"affinity": tokens[5],
                                   "rank": tokens[6]}
                elif version == "4":
                    result_dict = {"affinity": tokens[6],
                                   "rank": tokens[7]}
                nM = [peptides[i]]
                for value in sorted(scores):
                    nM.append(result_dict[value])
                affinities.append(tuple(nM))
        return affinities
    finally:
        # Remove temporary files
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)

def get_affinity_mhcnuggets(peptides, allele, version, dic, remove_files=True):
    """ Obtains binding affinities from list of peptides

        peptides: peptides of interest (list of strings)
        allele: Allele to use for binding affinity (string)
        scores: list of scoring methods
        version: version of mhcnuggets
        dic: path to available alleles dictionary
        remove_files: option to remove intermediate files

        Return value: affinities (a list of binding affinities
                        as strings)
    """
    from mhcnuggets.src.predict import predict
    files_to_remove = []
    try:
        # Check that allele is valid for method
        with open(dic, "rb") as allele_stream:
            avail_alleles = pickle.load(allele_stream)
        # Check that allele is valid for method
        allele = allele.replace("*", "")
        if allele in avail_alleles["mhcnuggets_mhcI"]:
            allele_class = "I"
            max_length = 15
        elif allele in avail_alleles["mhcnuggets_mhcII"]:
            allele_class = "II"
            max_length = 30
        else:
            warnings.warn(
                " ".join([allele, "is not a valid allele for mhcnuggets"]), Warning
            )
            return [(peptides[i], "NA") for i in range(0, len(peptides))]
        # Establish return list and sample id
        sample_id = ".".join(
            [peptides[0], str(len(peptides)), allele, "mhcnuggets", version]
        )
        affinities = []
        # Write one peptide per line to a temporary file for
        #   input if peptide length is at least 9
        # Count instances of smaller peptides
        # Establish temporary file to hold output
        peptide_file = tempfile.mkstemp(
            suffix=".txt", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(peptide_file)
        na_count = 0
        with open(peptide_file, "w") as f:
            for sequence in peptides:
                if len(sequence) > max_length or '?' in sequence:
                    na_count += 1
                else:
                    print(sequence, file=f)
        if na_count > 0:
            warnings.warn(
                " ".join(
                    [
                        str(na_count),
                        "peptides not compatible with",
                        "mhcnuggets will not receive score",
                    ]
                ),
                Warning,
            )
        # Establish temporary file to hold output
        mhc_out = tempfile.mkstemp(
            suffix=".mhcnuggets.out", prefix="".join([sample_id, "."]), text=True
        )[1]
        files_to_remove.append(mhc_out)
        # Run mhcnuggets
        predict(
            class_=allele_class, peptides_path=peptide_file, mhc=allele, output=mhc_out
        )
        # Retrieve scores for valid peptides
        score_dict = {}
        with open(mhc_out, "r") as f:
            # Skip headers
            f.readline()
            for line in f:
                tokens = line.strip("\n").split(",")
                score_dict[tokens[0]] = tokens[1]
        # Produce list of scores for valid peptides
        # Invalid peptides receive "NA" score
        for sequence in peptides:
            if sequence in score_dict:
                nM = (sequence, score_dict[sequence])
            else:
                nM = (sequence, "NA")
            affinities.append(nM)
        return affinities
    finally:
        if remove_files:
            for file_to_remove in files_to_remove:
                os.remove(file_to_remove)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tool', type=str, required=True,
            help='which tool to use for predictions'
        )
    parser.add_argument('-n', '--neoepiscope-file', type=str, required=True,
            help='Path to patient neoepiscope file'
        )
    parser.add_argument('-d', '--allele-dict', type=str, required=True,
            help='Path to pickled dictionary with available alleles from neoepiscope'
        )
    parser.add_argument('-a', '--allele', type=str, required=True,
            help='MHC allele'
        )
    parser.add_argument('-o', '--output-dir', required=True,
            help='Path to output directory'
        )
    args = parser.parse_args()

    # Get peptides
    peptides = set()
    with open(args.neoepiscope_file) as f:
    	f.readline()
    	f.readline()
    	for line in f:
    		tokens = line.strip().split('\t')
    		peptides.add(tokens[0])
    peptides = list(peptides)

    # Select tool and run
    if args.tool == 'mhcnuggets':
    	affinities = get_affinity_mhcnuggets(peptides, args.allele, "2", remove_files=True)
    elif args.tool == 'netMHCpan':
    	affinities = get_affinity_netMHCpan(peptides, args.allele, paths.netMHCpan4, "4", ["affinity"], remove_files=True)
    else:
    	raise NotImplementedError(''.join([args.tool, ' is not a valid tool']))

    # Store affinities as dict
    affinity_dict = {}
    for score in affinities:
    	affinity_dict[score[0]] = score[1]

    # Get patient ID
    patient_id = os.path.basename(args.neoepiscope_file).replace('.neoepiscope.comprehensive.out', '')

    # Pickle dictionary
    dict_file = os.path.join(args.output_dir, '.'.join([patient_id, args.tool, args.allele, 'pickle']))

    with open(dict_file, 'wb') as p:
    	pickle.dump(affinity_dict, p)


