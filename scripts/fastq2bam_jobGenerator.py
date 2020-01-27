#!/usr/bin/env python

from __future__ import print_function
import argparse
import json
import os
import re
import subprocess


def main(manifest, cwl, output_prefix, output_dir, cmd_output_dir, fastq_dir, reference_dir):
    i = 0
    fh = open(manifest, 'r')
    for line in fh:
        s = line.strip().split("\t")
        normal_sample = s[1]
        normal_fastq = [os.path.join(fastq_dir, ''.join([normal_sample, '_1.fastq.gz'])), os.path.join(fastq_dir, ''.join([normal_sample, '_2.fastq.gz']))]
        tumor_sample = s[2]
        tumor_fastq = [os.path.join(fastq_dir, ''.join([tumor_sample, '_1.fastq.gz'])), os.path.join(fastq_dir, ''.join([tumor_sample, '_2.fastq.gz']))]

        if normal_sample == "NA" or tumor_sample == "NA":
            print("warning: no tumor / normal pair found. tumor: %s normal: %s" % (tumor_sample, normal_sample))
            continue

        r = Run(tumor_fastq, tumor_sample, normal_fastq, normal_sample, reference_dir)

        # write cwl inputs
        input_json = r.generate_inputs_json()
        inputs_filename = os.path.join(
            output_dir, "".join([output_prefix, ".", str(i), ".", tumor_sample, ".json"])
        )
        write_file(json.dumps(input_json), inputs_filename)

        exec_file_contents = r.generate_cwltool_submission(tumor_sample, cwl, inputs_filename, cmd_output_dir)

        exec_filename = os.path.join(
            output_dir, "".join([output_prefix, ".", str(i), ".sh"])
        )
        write_file(exec_file_contents, exec_filename)
        os.chmod(exec_filename, 0755)
        i += 1
    fh.close()


def write_file(file_contents, output_filename):
    with open(output_filename, 'w') as ofh:
        ofh.write(file_contents)


class Run(object):

    def __init__(self, tumor_fastq, tumor_sample, normal_fastq, normal_sample, reference_dir):
        self.tumor_fastq = tumor_fastq
        self.tumor_sample = tumor_sample
        self.normal_fastq = normal_fastq
        self.normal_sample = normal_sample
        self.phase1 = os.path.join(reference_dir, "1000G_phase1.indels.hg19.sites.fixed.vcf.gz")
        self.dbsnp = os.path.join(reference_dir, "dbsnp_132_b37.leftAligned.vcf.gz")
        self.bwa_idx = os.path.join(reference_dir, "bwa_idx_GRCh37d5.tar.gz")
        self.reference = os.path.join(reference_dir, "core_ref_GRCh37d5/genome.fa")
        self.referenceTar = os.path.join(reference_dir, "core_ref_GRCh37d5.tar.gz")
        self.mills_and_1000g = os.path.join(reference_dir, "Mills_and_1000G_gold_standard.indels.hg19.sites.fixed.vcf.gz")


    def generate_inputs_json(self):
        json_template = {
            "tumor_fastqs": [],
            "tumor_sample_name": self.tumor_sample,
            "normal_fastqs": [],
            "normal_sample_name": self.normal_sample,
            "dbsnp": {
                "path": self.dbsnp,
                "class": "File"
            },
            "1000g_phase1": {
                "path": self.phase1,
                "class": "File"
            },
            "mills_and_1000g": {
                "path": self.mills_and_1000g,
                "class": "File"
            },
            "reference": {
                "path": self.reference,
                "class": "File"
            },
            "referenceTar": {
                "path": self.referenceTar,
                "class": "File"
            },
            "bwa_idx": {
                "path": self.bwa_idx,
                "class": "File"
            }
        }
        for s in self.tumor_fastq:
            json_template["tumor_fastqs"].append(
                {
                    "path": s,
                    "class": "File"
                }
            )
        for s in self.normal_fastq:
            json_template["normal_fastqs"].append(
                {
                    "path": s,
                    "class": "File"
                }
            )
        return json_template

    def generate_script_contents(self, command):
        file_parts = [
            "#!/bin/bash",
            "set -e",
            "\n".join(command),
        ]
        return "\n\n".join(file_parts)

    def generate_cwltool_submission(self, name, cwl, cwl_inputs, outdir):
        command = [
            "mkdir -p {}".format(
                os.path.join(outdir, name, "tmp")
            ),
            "mkdir -p {}".format(
                os.path.join(outdir, name, "cache")
            ),
            "\n",
            "cwltool \\",
            "--outdir %s \\" % (os.path.join(outdir, name)),
            "--tmpdir-prefix %s \\" % (os.path.join(outdir, name, "tmp")),
            "--cachedir %s \\" % (os.path.join(outdir, name, "cache")),
            "%s \\" % (cwl),
            "%s" % (cwl_inputs)
        ]
        return self.generate_script_contents(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_manifest",
                        type=str,
                        help="tsv file containing the fields: Subject_ID, Normal_sample, Tumor_sample, RNA_sample")
    parser.add_argument("--cwl", "-c",
                        type=str,
                        default="fastq2bam.cwl.yaml",
                        help="cwl workflow")
    parser.add_argument("--output-prefix", "-p",
                        type=str,
                        default="fastq2bam",
                        help="prefix to name output files with")
    parser.add_argument("--output-dir", "-o",
                        type=str,
                        default=".",
                        help="output directory to write generated files to")
    parser.add_argument("--cmd-output-dir", "-O",
                        type=str,
                        default=".",
                        help="output directory for workflow to write files to")
    parser.add_argument("--fastq-dir", "-f",
                        type=str,
                        help="path to directory containing gzipped fastq files")
    parser.add_argument("--reference-dir", "-r",
                        type=str,
                        help="path to directory containing reference data")
    args = parser.parse_args()

    if args.output_dir is not None:
        if not os.path.isdir(os.path.abspath(args.output_dir)):
            subprocess.Popen(
                "mkdir -p {0}".format(os.path.abspath(args.output_dir)),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

    if args.cmd_output_dir is not None:
        if not os.path.isdir(os.path.abspath(args.cmd_output_dir)):
            subprocess.Popen(
                "mkdir -p {0}".format(os.path.abspath(args.cmd_output_dir)),
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

    if args.output_prefix is not None:
        if re.search("/", args.output_prefix):
            raise ValueError("[ERROR] The output prefix cannot contain '/'")
        
    main(args.input_manifest, 
         os.path.abspath(args.cwl),
         args.output_prefix, 
         os.path.abspath(args.output_dir),
         os.path.abspath(args.cmd_output_dir), 
         os.path.abspath(args.fastq_dir)
         os.path.abspath(args.reference_dir))
