#!/usr/bin/env python

from __future__ import print_function
import argparse
import json
import os
import re
import subprocess
import uuid
import logging

def main(manifest, cwl, output_prefix, output_dir, cmd_output_dir, bam_dir, reference_dir):
    i = 0
    fh = open(manifest, 'r')
    for line in fh:
        s = line.strip().split("\t")
        normal_sample = s[1]
        tumor_sample = s[2]
        normal_bam = os.path.join(bam_dir, tumor_sample, ''.join([normal_sample, '.reheadered.realigned.cleaned.bam']))
        tumor_bam = os.path.join(bam_dir, tumor_sample, ''.join([tumor_sample, '.reheadered.realigned.cleaned.bam']))
        
        if normal_sample == "NA" or tumor_sample == "NA":
            print("warning: no tumor / normal pair found. tumor: %s normal: %s" % (tumor_sample, normal_sample))
            continue

        r = Run(tumor_bam, tumor_sample, normal_bam, normal_sample)

        # write cwl inputs
        input_json = r.generate_inputs_json()
        inputs_filename = os.path.join(
            output_dir, "".join([output_prefix, ".", str(i), ".", tumor_sample, ".json"])
        )
        write_file(json.dumps(input_json), inputs_filename)

        # Generate a bash script
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

    def __init__(self, tumor_bam, tumor_sample, normal_bam, normal_sample):
        self.tumor_bam = tumor_bam
        self.tumor_sample = tumor_sample
        self.normal_bam = normal_bam
        self.normal_sample = normal_sample
        self.dbsnp = os.path.join(reference_dir, "dbsnp_132_b37.leftAligned.vcf.gz")
        self.cosmic = os.path.join(reference_dir, "b37_cosmic_v54_120711.vcf.gz")
        self.centromere = os.path.join(reference_dir, "centromere_hg19.bed")
        self.reference = os.path.join(reference_dir,"genome.fa.gz")
        self.referenceTar = os.path.join(reference_dir, "core_ref_GRCh37d5.tar.gz")
        self.bed_file = os.path.join(reference_dir, "gaf_20111020+broad_wex_1.1_hg19.bed")
        if os.path.isfile(self.tumor_bam) and os.path.isfile(self.normal_bam):
            subprocess.call(['cp', self.tumor_bam.replace('.bam', '.bai'), ''.join([self.tumor_bam, '.bai'])])
            subprocess.call(['cp', self.normal_bam.replace('.bam', '.bai'), ''.join([self.normal_bam, '.bai'])])

    def generate_inputs_json(self):
        json_template = {
            "tumor": {
                "path": self.tumor_bam,
                "class": "File"
           },
            "tumor_aliquot_uuid": str(uuid.uuid4()),
            "tumor_analysis_uuid": str(uuid.uuid4()),
            "tumor_aliquot_name": self.tumor_sample,
            "normal": {
                "path": self.normal_bam,
                "class": "File"
            },
            "normal_aliquot_uuid": str(uuid.uuid4()),
            "normal_analysis_uuid": str(uuid.uuid4()),
            "normal_aliquot_name": self.tumor_sample,
            "dbsnp": {
                "path": self.dbsnp,
                "class": "File"
            },
            "cosmic": {
                "path": self.cosmic,
                "class": "File"
            },
            "centromere": {
                "path": self.centromere,
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
            "bed_file": {
                "path": self.bed_file,
                "class": "File"
            },
            "platform": "illumina",
            "center": "OHSU"
        }
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
                        default="mc3_variant.cwl",
                        help="cwl workflow")
    parser.add_argument("--output-prefix", "-p",
                        type=str,
                        default="bam2variants_mc3",
                        help="prefix to name output files with")
    parser.add_argument("--output-dir", "-o",
                        type=str,
                        default=".",
                        help="output directory to write generated files to")
    parser.add_argument("--cmd-output-dir", "-O",
                        type=str,
                        default=".",
                        help="output directory for workflow to write files to")
    parser.add_argument("--bam-dir", "-b",
                        type=str,
                        default=".",
                        help="path to directory containing BAM files")
    parser.add_argument("--reference-dir", "-r",
                        type=str,
                        default=".",
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
         os.path.abspath(args.bam_dir),
         os.path.abspath(args.reference_dir))
