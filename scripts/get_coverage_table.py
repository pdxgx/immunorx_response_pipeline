#!/usr/bin/env python

from __future__ import print_function
import os
import glob
from collections import defaultdict

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("--graph-dir", "-g",
                            type=str, required=True,
                            help="path to directory containing bedgraph files")
        parser.add_argument("--output-file", "-o",
                            type=str, required=True,
                            help="path to write output file")
        parser.add_argument("--manifest", "-m",
                            type=str, required=True,
                            help="path to tumor-normal pair manifest file")
        args = parser.parse_args()

        graph_dir = args.graph_dir
        out_file = args.output_file

        patients = []
        with open(out_file, 'w') as o:
                headers = ['Patient', 'Tumor', '6']
                print('\t'.join(headers), file=o)
                with open(args.manifest) as f:
                        for line in f:
                                tokens = line.strip().split('\t')
                                bg = os.path.join(graph_dir, '.'.join([tokens[0], tokens[2], 'bg']))
                                coverage = 0
                                with open(bg) as g:
                                        for line in g:
                                                gtokens = line.strip().split('\t')
                                                bp = float(gtokens[2]) - float(gtokens[1])
                                                reads = float(gtokens[3])
                                                if reads >= 6:
                                                        coverage += bp
                                coverage = coverage/1000000.0
                                out_list = [tokens[0], tokens[2], str(coverage)]
                                print('\t'.join(out_list), file=o)
