# Immunotherapy response pipeline

Pipeline for reproducing the workflow used to identify predictors of immunotherapy response as described in our [manuscript](https://www.biorxiv.org/content/10.1101/665026v2). 

## Software requirements

*Python and python packages:*

We ran python scripts with python 2.7.14, installed using [Anaconda](https://www.anaconda.com/distribution/). We made use of the following python packages, installable with `pip`:

* [`astropy` v3.1.2](https://github.com/astropy/astropy)
* [`biopython` v1.75](https://github.com/biopython/biopython)
* [`cwltool` v1.0.20180211183944](https://github.com/common-workflow-language/cwltool)
* [`neoepiscope` v0.3.5](https://github.com/pdxgx/neoepiscope)
* [`pyliftover` v0.4](https://github.com/konstantint/pyliftover)

You will also need `Open-CRAVAT`, which is installable via `pip`, but will require a virtual environment with python 3.6 or higher. Use the following commands to create and set up the environment (including installing `intervaltree` v2.1.0 for another analysis):

```
conda create -n cravat_env python=3.6

conda activate cravat_env==1.6.1

pip install intervaltree==2.1.0

pip install open-cravat

cravat-admin install-base

cravat-admin ls -a -t clinvar

conda deactivate
```

*R and R packages*

We ran `R` scripts with [`R` v3.6.1](https://ftp.osuosl.org/pub/cran/src/base/R-3/R-3.6.1.tar.gz). We used the following packages installable from CRAN:

* `car` (v3.0-3)
* `data.table` (v1.12.8)
* `devtools` (v2.2.1)
* `ggplot2` (v3.2.1)
* `gridExtra` (v2.3)
* `magritter` (v1.5) 
* `pROC` (v1.15.3)
* `RColorBrewer` (v1.1-2)
* `scales` (v1.0.0)
* `survival` (v2.44-1.1)
* `survminer` (v0.4.4)

And the following package installable through `devtools`:

* [`kma`](https://github.com/adamtongji/kma) (commit 01ec25f575e4f4c8b1c90cd6223d2d8e9bc8e825)

We also used [`R studio`(https://rstudio.com/products/rstudio/download/#download) for interactively generating figures and analyzing results.

*CWL workflows and their dependencies:*

We used the following workflows:

* [dockstore-cgpmap](https://github.com/cancerit/dockstore-cgpmap) (commit 0bacb0bee2e5c04b268c629d589ff1c551d34745)
* [gatk-cocleaning-tool](https://github.com/OpenGenomics/gatk-cocleaning-tool) (commit d2bafc23221f6a8dceedd45a534163e0e1bf5c68)
* [mc3](https://github.com/OpenGenomics/mc3) (commit 72a24b55544e3011ede1c46b13d531a7d05ef4e0)

These workflows require [`docker`](https://docs.docker.com/v17.12/install/) to be installed.

They also require the following perl modules, installable with [cpan minus](http://www.cpan.org/modules/INSTALL.html):

* `Getopt::Long`
* `File::Path`
* `File::Temp`
* `File::Spec`
* `IO::File`
* `Pod::Usage`
* `Data::Dumper`

*Repository to clone:*

You will need to clone [this repository](https://github.com/JulianneDavid/shared-cancer-splicing) for use in tumor-specific junction detection.

*Additional required software:*

We also installed the following tools, and any dependencies they listed in their installation instructions:

* [`bedtools` v2.23.0](https://github.com/arq5x/bedtools2/releases/tag/v2.23.0)
* [`blast` v2.6.0](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/)
* [`bowtie2` v2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
* [`express` v1.5.1](http://bio.math.berkeley.edu/eXpress/)
* [`GATK` v3.7](https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.7-0-gcfedb67)
* [`HapCUT2`](https://github.com/vibansal/HapCUT2) (commit 1c90cbb207f0834acc1ba778515251df9aa0eab2)
* [`liftOver`](https://genome-store.ucsc.edu/) (add to cart and download, free for academic use; our version was downloaded in March 2019)
* [`mSINGS`](https://bitbucket.org/uwlabmed/msings/src/master/) (commit 030289381f3b7aee24d8eccbb69b3e66711f5bb0)
* [`NetCTLpan` v1.1](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netCTLpan)
* [`OptiType`](https://github.com/FRED-2/OptiType) (commit 9de4a0f82e2e6e472797e7cec114f60a57e61d2a)
* [`seq2hla`](https://bitbucket.org/sebastian_boegel/seq2hla/src/default/) (commit c37148c10e8569c37bc846f618fc5f6690c7990c)
* [`samtools` v1.6](https://github.com/samtools/samtools/releases/tag/1.6)
* [`STAR` vv2.6.1c](https://github.com/alexdobin/STAR/releases/tag/2.6.1c)
* [`tabix` v1.4](https://github.com/samtools/htslib/releases/tag/1.4)
* [`vt` v0.5772](https://github.com/atks/vt/releases/tag/0.5772)


## Data

### Patient genomic data

We downloaded the paired-end whole exome sequencing (WES) and RNA sequencing (RNA-seq) FASTQ files from the following BioProjects from the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra):

* PRJNA278450
* PRJNA293912
* PRJNA305077
* PRJNA306070
* PRJNA307199
* PRJNA312948
* PRJNA324705
* PRJNA343789
* PRJNA357321
* PRJNA369259
* PRJNA414014
* PRJNA420786
* PRJNA82745

We also downloaded from the SRA RNA-seq reads from melanocytes:

* PRJNA421623

Note that some of these BioProjects require dbGAP access, so you will need to request dbGAP access to obtain them if you do not already have access.

We also downloaded the paired-end WES FASTQ files from the European Genome-phenome Archive project [EGAD00001004352](https://www.ebi.ac.uk/ega/datasets/EGAD00001004352), which also requires special access.

Finally, we obtained paired-end WES FASTQ data from [Graff et al.](http://www.oncotarget.com/index.php?journal=oncotarget&page=article&op=view&path[]=10547&pubmed-linkout=1) and [Le et al.](https://science.sciencemag.org/content/357/6349/409) by reaching out to the authors.

Place all the FASTQ files in the same directory. Each FASTQ will need to be processed for compatibility with the alignment workflow (see below):

```zcat [INPUT_FASTQ] | sed -E 's~^(@.+[.][0-9]+)[.]([12]) .+( length=[0-9]+)~\\1/\\2\\3~' | gzip -f > [OUTPUT_FASTQ]```

`[INPUT_FASTQ]` is the original FASTQ file, and `[OUTPUT_FASTQ]` is the processed FASTQ, named with the original prefix and with file ending `_1.fastq.gz` or `_2.fastq.gz`.

A [manifest file](data/immunorx_response_pairs.txt) summarizing tumor-normal pairs for samples used in the study is available and used in many of the analysis steps below. Additionally, a list of [melanocyte sample accession numbers](data/melanocyte_accessions.txt) is available for retained intron analysis (see "Tumor-specific retained intron identification" below).


### Genomic annotation data

*For the WES alignment/bam processing and somatic variant calling workflows:*

We downloaded the [Human GRCh37d5 reference bundle](ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/core_ref_GRCh37d5.tar.gz) and [`bwa` index](ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/bwa_idx_GRCh37d5.tar.gz) for use with the alignment workflows. You will need to both retain a copy of the `core_ref_GRCh37d5.tar.gz` file, as well as decompress it and retain a copy of the `genome.fa.gz` it contains. Store all of these in the same directory.

We also downloaded some annotation data from the Broad Institute [b37 resource bundle](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/): the [1000 genomes indels](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz), the [Mills and 1000 genomes indels](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz), and the [dbSNP variants](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz). All of these compressed VCF files should be indexed with tabix.

We also downloaded some of the [FireCloud reference files from the TCGA](https://console.cloud.google.com/storage/browser/firecloud-tcga-open-access): `hg19_cosmic_v54_120711.vcf` and `gaf_20111020+broad_wex_1.1_hg19.bed`. You will need permissions to access these files

Store the GRCh37d5 reference bundle files, `bwa` index files, Broad reference files, and TCGA FireCloud reference files in one directory, and copy the [centromere bed file](data/centromere_hg19.bed) in this repo to that directory as well.

*For RNA-seq alignments:*

We downloaded the GENCODE [GRCh37](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz) and [GRCh38](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz) reference FASTA files, and the [GRCh37](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz) and [GRCh38](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz) GTF files.

*For DNA neoepitope prediction:*

We ran `neoepiscope`'s `download` functionality, answering yes to downloading/indexing the GENCODE v19 annotation, and yes to downloading the hg19 bowtie index.

*For tumor-specific splice junction identification:*

We used some data described in the "Data" section of [this repository](github.com/JulianneDavid/shared-cancer-splicing). We downloaded all of the exon-exon junction BEDs for GTEx and TCGA, as well as the GENCODE v28 GTF file.

*For extended neoepitope burden analysis:*

We downloaded the [hg19-to-hg38 chain file](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz) from UCSC and the [GRCh37 protein coding translations FASTA file](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_translations.fa.gz) from GENCODE.

We also downloaded TPM expression data for TCGA generated by the National Cancer Institute, including [COAD data](https://osf.io/95wnv/), [KIRC data](https://osf.io/egxyw/), [LUAD data](https://osf.io/3yngu/), [LUSC data](https://osf.io/tyfha/), [PRAD data](https://osf.io/m5nh6/), [SKCM data](https://osf.io/cxj8h/), [THCA data](https://osf.io/7pdnr/), and [UCEC data](https://osf.io/h2tu6/) cancer types, as well as the [ID map](https://osf.io/7qpsg/) for linking CGHubAnalysis UUIDs to TCGA aliquot barcodes. Unzip these files and store them all in one directory.

*For survival analysis:*

We obtained data from TCGA for the SKCM and KIRC cancer types. We downloaded [SKCM clinical data](http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/SKCM/20160128/gdac.broadinstitute.org_SKCM.Merge_Clinical.Level_1.2016012800.0.0.tar.gz) and [KIRC clinical data](http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KIRC/20160128/gdac.broadinstitute.org_KIRC.Merge_Clinical.Level_1.2016012800.0.0.tar.gz). We also downloaded [SKCM MAF files](http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/SKCM/20160128/gdac.broadinstitute.org_SKCM.Mutation_Packager_Calls.Level_3.2016012800.0.0.tar.gz) and [KIRC MAF files](http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KIRC/20160128/gdac.broadinstitute.org_KIRC.Mutation_Packager_Calls.Level_3.2016012800.0.0.tar.gz). Unzip all the files.

*For other data processing:*

We downloaded the [HLA reference FASTA file](https://raw.githubusercontent.com/FRED-2/OptiType/master/data/hla_reference_dna.fasta) from Optitype.


## Execution

### WES alignment/BAM processing

We used the dockstore-cgpmap workflow to align reads to the GRCh37d5 reference, generating genome-coordinate sorted alignments with duplicates marked, and the gatk-cocleaning-tool workflow to realign around indels and perform base recalibration for paired tumor and normal sequence read data. The [fastq2bam_jobGenerator.py](scripts/fastq2bam_jobGenerator.py) script generates the required input json files and shell scripts to execute the workflows:

```python fastq2bam_jobGenerator.py -c [WORKFLOW_FILE] -o [SCRIPT_DIR] -O [BAM_DIR] -f [FASTQ_DIR] -r [REFERENCE_DIR] [MANIFEST]```

The `WORKFLOW_FILE` is the path to the [fastq2bam.cwl.yaml](scripts/fastq2bam.cwl.yaml) file in this repository - you will need to copy this file to the same directory that includes your dockstore-cgpmap and gatk-cocleaning-tool repos. You will also need to copy the [rename.cwl.yaml](scripts/rename.cwl.yaml) file to the same directory. The `SCRIPT_DIR` is the path to the directory this python script will generate input json files and shell scripts to, and the `BAM_DIR` is the file that processed tumor and normal alignments will be written to (within subdirectories named for each tumor sample) after the shell scripts in `SCRIPT_DIR` are run. `FASTQ_DIR` is the path to the directory containing your FASTQ files, and `REFERENCE_DIR` is the path to the directory which contains your reference data for the workflow (see the "For the WES alignment/bam processing and somatic variant calling workflows" subsection of "Genomic annotation data" above). The `MANIFEST` is the path to the tumor-normal pair manifest file described in "Patient genomic data" above.

We ran all of the generated shell scripts in `SCRIPT_DIR`.


### Somatic variant calling and consensus calling

We used the mc3 workflow to call somatic variants using MuSE, MuTect, Pindel, RADIA, SomaticSniper, and VarScan 2. The [bam2variants_jobGenerator.py](scripts/bam2variants_jobGenerator.py) script generates the required input json files and shell scripts to execute the workflows.

```python bam2variants_jobGenerator.py  -c [MC3_WORKFLOW_FILE] -o [SCRIPT_DIR] -O [VARIANT_DIR] -b [BAM_DIR] -r [REFERENCE_DIR] [MANIFEST]```

The `MC3_WORKFLOW_FILE` is the path to your `mc3_variant.cwl` file in your mc3 repository. The `SCRIPT_DIR` is the path to the directory this python script will generate input json files and shell scripts to, and the `VARIANT_DIR` is the path to the directory that VCFs will be written to (within subdirectories for each tumor sample) after the shell scripts in `SCRIPT_DIR` are run. `BAM_DIR` is the path to the directory containing your alignment files (which live within subdirectories for each tumor, see "WES alignment/BAM processing" above), and `REFERENCE_DIR` is the path to the directory which contains your reference data (see the "For the WES alignment/bam processing and somatic variant calling workflows" subsection of "Genomic annotation data" above). The `MANIFEST` is the path to the tumor-normal pair manifest file described in "Patient genomic data" above.

We ran all of the generated shell scripts in `SCRIPT_DIR`.

After running mc3 then used the software `vt` to process the VCF files from each caller, *EXCEPT* those from `Indelocator` and `VarScan 2` indels:

First, we normalized each VCF:

```vt normalize [INPUT_VCF] -n -r [REFERENCE_FASTA] -o [OUTPUT_VCF]```

`INPUT_VCF` is the path to the VCF for that caller from the mc3 workflow, `REFERENCE_FASTA` is the path to the FASTA file in your reference data directory (see the "For the WES alignment/bam processing and somatic variant calling workflows" subsection of "Genomic Annotation Data" above), and `OUTPUT_VCF` is the path to which `vt normalize` will write the normalized VCF.

Then we decomposed block substitutions for each normalized VCF:

```vt decompose_blocksub -a -p [INPUT_VCF] -o [OUTPUT_VCF]```

`INPUT_VCF` is the path to the VCF for that caller after running `vt normalize`, and `OUTPUT_VCF` is the path to which `vt decompose_blocksub` will write the normalized/partially decomposed VCF.

Then we decomposed multi-allelic variants for each VCF from the previous step:

```vt decompose -s [INPUT_VCF] -o [OUTPUT_VCF]```

`INPUT_VCF` is the path to the VCF for that caller after running `vt decompose_blocksub`, and `OUTPUT_VCF` is the path to which `vt decompose` will write the normalized/fully decomposed VCF.

Then, for *`MuTect` variants ONLY*, we sorted variants:

```vt sort [INPUT_VCF] -o [OUTPUT_VCF]```

`INPUT_VCF` is the pat to the VCF for that caller from `vt decompose`, and `OUTPUT_VCF` is the path to which `vt sort` will write the normalized/fully decomposed/sorted VCF.

Then we took unique variants for all callers:

```vt uniq [INPUT_VCF] -o [OUTPUT_VCF]```

`INPUT_VCF` is the path to the VCF for that caller from `vt decompose` (or `vt sort` for `MuTect` VCFs), and `OUTPUT_VCF` is the path to which `vt uniq` will write the final VCF for that caller.

Finally, we used a [python script](scripts/VCF_parse.py) to produce a consensus call set of variants produced by at least 2 callers and not overlapped by a Pindel variant, or called by solely Pindel:

```python VCF_parse.py -v [MUSE_VCF],[MUTECT_VCF],[PINDEL_VCF],[RADIA_VCF],[SOMATICSNIPER_VCF],[VARSCAN_FPF_VCF] -c muse,mutect,pindel,radia,somaticsniper,varscan -o [OUTPUT_DIR] -n 2 -s [PATIENT_ID].[TUMOR_ID], -t [TUMOR_ID] -i 1000```

The VCFs listed for the `-v` option are the paths to those VCFs processed by `vt uniq` for each caller, the `OUTPUT_DIR` is the path to the directory to which the script will write consensus VCFs, the `PATIENT_ID` is the patient identifier for the sample (first column of the tumor-normal pair manifest file described in "Patient genomic data" above), and `TUMOR_ID` is the tumor identifier for the sample (third column of the tumor-normal pair manifest file described in "Patient genomic data" above).


### Germline variant calling

We used `GATK` to run `HaplotypeCaller` for calling germline variants, and to run `VariantFiltration` for filtering the results of `HaplotypeCaller`.

HaplotypeCaller was run using default options:

```java -jar GenomeAnalysisTK.jar -R [REFERENCE_FASTA] -T HaplotypeCaller -I [NORMAL_BAM] -o [RAW_GERMLINE_VCF]```

The `REFERENCE_FASTA` is the path to the FASTA file in your reference data directory (see the "For the WES alignment/bam processing and somatic variant calling workflows" subsection of "Genomic Annotation Data" above). The `NORMAL_BAM` file is the path to the processed normal sample alignment file (see "WES alignment/BAM processing" above). `[RAW_GERMLINE_VCF]` is the path to which `HaplotypeCaller` will write the raw VCF.

`GATK`'s `VariantFiltration` was run as below:

```java -jar GenomeAnalysisTK.jar -R [REFERENCE_FASTA] -T VariantFiltration --variant [RAW_GERMLINE_VCF] -o [FILTERED_GERMLINE_VCF] --clusterSize 3 --clusterWindowSize 15 --missingValuesInExpressionsShouldEvaluateAsFailing --filterName 'QDFilter' --filterExpression 'QD < 2.0' --filterName 'QUALFilter' --filterExpression 'QUAL < 100.0' --filterName DPFilter --filterExpression 'DP < 10.0'```

The `REFERENCE_FASTA` is the path to the FASTA file in your reference data directory (see the "For the WES alignment/bam processing and somatic variant calling workflows" subsection of "Genomic Annotation Data" above). `[RAW_GERMLINE_VCF]` is the path to the raw VCF produced by `HaplotypeCaller`. `FILTERED_GERMLINE_VCF` is the path to which `VariantFiltration` will write the VCF with filtering information.

Only variants passing all filters were retained and saved to a final for downstream analyses:

```grep '^#' [FILTERED_GERMLINE_VCF] > [POSTFILTER_GERMLINE_VCF]```

```grep -v '^#' [FILTERED_GERMLINE_VCF] | grep PASS >> [POSTFILTER_GERMLINE_VCF]```

`FILTERED_GERMLINE_VCF` is the path to the VCF with filtering information from `VariantFiltration`, and `POSTFILTER_GERMLINE_VCF` is the path to the cleaned VCF containing only variants which passed all filters.

 We also ran `vt` to normalize, decompose (both block substitutions and multi-allelic variants), sort, and obtain a unique list of variants from the filtered/cleanred germline VCFs (as described above for the VCFs from individual somatic variant callers in "Somatic variant calling and consensus calling"). On the VCF resulting from `vt`'s `uniq`, we ran one final filtering step to correct genotyping formats:

```grep -v '^#' [UNIQ_GERMLINE_VCF] | sed -E 's/1[/][.][:]/0\/1:/' | sed -E 's/[.][/]1[:]/0\/1:/' | sed -E 's/[.]\/[.]/0\/0/'  >> [FINAL_GERMLINE_VCF]```

`UNIQ_GERMLINE_VCF` is the path to the the VCF from `vt uniq`, and `FINAL_GERMLINE_VCF` is the path to which to write the final germline VCF file for downstream analyses.


### Genome coverage

We used `bedtools genomecov` to determine the Mbp of genome covered in each tumor sample:

```bedtools genomecov -ibam [TUMOR_BAM] -bg > [BEDGRAPH_DIR][PATIENT_ID].[TUMOR_ID].bg```

The `TUMOR_BAM` is the processed tumor sample alignment file (see "WES alignment/BAM processing" above). `BEDGRAPH_DIR` is the path to the directory to which raw coverage output will be written. `PATIENT_ID` is the patient identifier for the sample (first column of the tumor-normal pair manifest file described in "Patient genomic data" above), and `TUMOR_ID` is the tumor identifier for the sample (third column f the tumor-normal pair manifest file described in "Patient genomic data" above). 

We processed the output BedGraph files for all samples to determine the number of bp pairs per sample that were covered by at least 6 reads. A coverage summary file for all patients was generated with a [python script](scripts/get_coverage_table.py):

```python get_coverage_table.py --graph-dir [BEDGRAPH_DIR] --output-file [OUTPUT_FILE] --manifest [MANIFEST]```

`BEDGRAPH_DIR` is the path to the directory containing the bedgraph files produced by `bedtools genomecov` in the previous step, `OUTPUT_FILE` is the path to which the script will write the coverage table, and `MANIFEST` is the path to the tumor-normal pair manifest file described in "Patient genomic data" above.


### Haplotype phasing

Before performing haplotype phasing, we merged our final germline variants and our consensus somatic variants using `neoepiscope`:

```neoepiscope merge -g [FINAL_GERMLINE_VCF] -s [SOMATIC_VCF] -o [MERGED_VCF]```

`FINAL_GERMLINE_VCF` is the path to the final germline VCF (see "Germline variant calling" above), `SOMATIC_VCF` is the path to the consensus somatic VCF produced by the [VCF_parse.py](scripts/VCF_parse.py) script (see "Somatic variant calling and consensus calling" above), and `MERGED_VCF` is the path to which `neoepiscope merge` will write the merged VCF.

Then, we predicted tumor haplotypes using `HapCUT2`'s `extractHAIRS` and `HAPCUT2` functionalities:

```extractHAIRS --indels 1 --bam [TUMOR_BAM] --VCF [MERGED_VCF] --out [FRAGMENT_FILE]```

```HAPCUT2 --fragments [FRAGMENT_FILE] --vcf [MERGED_VCF] --output [HAPLOTYPES]```

`TUMOR_BAM` is the path to the processed tumor sample alignment file (see "WES alignment/BAM processing" above). `MERGED_VCF` is the path to the merged VCF produced by `neoepiscope merge`. `FRAGMENT_FILE` is the path to the output of `HapCUT2`'s `extractHAIRS`, and `HAPLOTYPES` is the path to the output of `HAPCUT2`.

Finally, we prepared our haplotype predictions for neoepitope prediction using `neoepiscope`:

```neoepiscope prep -v [MERGED_VCF] -c [HAPLOTYPES] -o [PREPPED_HAPLOTYPES]```

`MERGED_VCF` is the path to the merged VCF produced by `neoepiscope merge`, `HAPLOTYPES` is the path to the output of `HAPCUT2`, and `PREPPED_HAPLOTYPES` is the path to which `neoepiscope prep` will write the prepared haplotype data.


### HLA typing

We performed MHC Class I HLA typing with `OptiType` for each tumor WES FASTQ pair as below:

```python OptiTypePipeline.py -i [FASTQ1] [FASTQ2] --dna --outdir [OUTPUT_DIR] --prefix [TUMOR_ID]```

`[FASTQ1]` and `[FASTQ2]` are the paths to the forward and reverse FASTQ files for the tumor sample, and `TUMOR_ID` is the tumor identifier for the sample (third column f the tumor-normal pair manifest file described in "Patient genomic data" above). `OUTPUT_DIR` is the path to the directory where `OptiType` will write the output files.

We performed MHC Class II HLA typing with `seq2hla` for each tumor WES FASTQ pair:

```python seq2HLA.py  -1 [FASTQ1] -2 [FASTQ2] -r [OUTPUT_DIR]/[TUMOR_ID]```

`[FASTQ1]` and `[FASTQ2]` are the paths to the forward and reverse FASTQ files for the tumor sample, and `TUMOR_ID` is the tumor identifier for the sample (third column f the tumor-normal pair manifest file described in "Patient genomic data" above). `OUTPUT_DIR` is the path to the same directory where you had `OptiType` write output files.


### MSI calculations

We determined tumor MSI status using `mSINGS`. You will first need to edit the `run_msings.sh` file in the `scripts` subdirectory of the `mSINGS` repository to include paths to your `VarScan 2` jar file (the `VARSCAN` variable), `samtools` executable (the `samtools` variable), and output directory (`SAVEPATH` variable). You will also need to generate a list of paths to processed tumor sample alignments, one per line. Then, we ran this command to determine MSI status across all patients:

```sh run_msings.sh [BAM_LIST] [INTERVALS_FILE] [BEDFILE] [REFERENCE_FASTA] [BASELINE_FILE]```

`BAM_LIST` is the path to your list of tumor sample alignments for each patient. `INTERVALS_FILE` is path to the `mSINGS_TCGA.msi_intervals` file in the `doc` subdirectory of the `mSINGS` repository. `BEDFILE` is the path to the `mSINGS_TCGA.bed` file in the `doc` subdirectory of the `mSINGS` repository. `REFERENCE_FASTA` is the path to the FASTA file in your reference data directory (see the "For the WES alignment/bam processing and somatic variant calling workflows" subsection of "Genomic Annotation Data" above). `BASELINE_FILE` is the path to the `mSINGS_TCGA.baseline` file in the `doc` subdirectory of the `mSINGS` repository. Data will be written into a subdirectory named for each alignment in the `SAVEPATH` directory.


### DNA neoepitope prediction

We predicted neoepitopes of 8-24 amino acids in length for each tumor sample from our consensus somatic variants with `neoepiscope call`, both accounting for phasing and germline and somatic variants with a comprehensive suite of transcript types, and not accounting for phasing with only protein-coding transcripts:

```neoepiscope call -b hg19 -c [PREPPED_HAPLOTYPES] -o [OUTPUT_DIR]/[PATIENT_ID].[TUMOR_ID].neoepiscope.comprehensive.out -k 8,24 --nmd --pp --igv --trv --no-affinity --fasta```

```neoepiscope call -b hg19 -c [PREPPED_HAPLOTYPES] -o [OUTPUT_DIR]/[PATIENT_ID].[TUMOR_ID].neoepiscope.somatic.out -k 8,24 --isolate --germline exclude --no-affinity --fasta```

`PREPPED_HAPLOTYPES` is the path to output of `neoepiscope prep` (see "Haplotype phasing" above). `OUTPUT_DIR` is the path the directory to which `neoepiscope` will write the output files. `PATIENT_ID` is the patient identifier for the sample (first column of the tumor-normal pair manifest file described in "Patient genomic data" above), and `TUMOR_ID` is the tumor identifier for the sample (third column f the tumor-normal pair manifest file described in "Patient genomic data" above). This will create output in your working directory, so you may add a directory as a prefix to the `PATIENT_ID` to direct output to a specific location.

Additionally, we predicted neopitopes of 8-11 amino acids in length from individual variant calling tools (`MuSE`, `MuTect`, `Pindel`, `RADIA`, `SomaticSniper`, and `VarScan 2` SNV) with `neoepiscope`, not accounting for phasing with only protein-coding transcripts. For each caller, we first removed variants that did not pass filters - for `MuTect`, `SomaticSniper`, and `VarScan 2` we retained variants with "PASS" in the `FILTER` field; for `MuSE` we retained variants that did NOT have "Tier5" in the `FILTER` field; and we retained all `Pindel` and `RADIA` variants, as they are already filtered. Then we ran `neoepiscope prep` and `neoepiscope call` for each caller/sample:

```neoepiscope prep -v [CLEAN_VCF] -o [PATIENT_ID].[TUMOR_ID].[CALLER].haplotypes```

```neoepiscope call -b hg19 -c [PATIENT_ID].[TUMOR_ID].caller.haplotypes -o [BY_CALLER_EPITOPE_DIR]/[PATIENT_ID].[TUMOR_ID].neoepiscope.[CALLER].out --no-affinity```

`CLEAN_VCF` is the path to the filtered VCF for each caller. `BY_CALLER_EPITOPE_DIR` is the path to the directory to which neoepiscope will write the output `PATIENT_ID` is the patient identifier for the sample (first column of the tumor-normal pair manifest file described in "Patient genomic data" above), and `TUMOR_ID` is the tumor identifier for the sample (third column f the tumor-normal pair manifest file described in "Patient genomic data" above). `CALLER` is the variant caller used (formatted as "muse", "mutect", "pindel", "radia", "somatic_sniper", or "varscan"). 

We then processed the neoepitope predictions from each caller to identify differences/similarities between tools using a [python script](process_epitopes_by_caller.py):

```python process_epitopes_by_caller.py -p [MANIFEST] -o [OUTPUT_DIR] -e [BY_CALLER_EPITOPE_DIR] -c [CONSENSUS_EPITOPE_DIR] ```

`MANIFEST` is the path to the tumor-normal pair manifest file described in "Patient genomic data" above, `OUTPUT_DIR` is the path to the directory to which the script will write files, `BY_CALLER_EPITOPE_DIR` is the directory containing the neoepitope predictions for each caller (see above), and `CONSENSUS_EPITOPE_DIR` is the path to the directory containing unphased neoepitope predictions based on consensus variant calls (see above).

Finally, to prepare data for plotting, we used another [python script](scripts/epitope_upset_combos.py) to simplify the data:

```python epitope_upset_combos.py -i [OUTPUT_DIR]/epitope_upset_table.tsv -o [OUTPUT_DIR]```

`OUTPUT_DIR` is the path to the same output directory used in the previous step. The script will write two files, `epitope_upset.tsv` and `epitope_upset_counts.tsv`, to the `OUTPUT_DIR`.


### Binding affinity analysis

To identify binding neoepitopes, we used `MHCnuggets` to predict binding affinities of phased neoepitopes. For all patient MHC Class I or II alleles for each patient as predicted by `Optitype` or `seq2hla`, we used a [python script](scripts/get_binding_scores.py) to predict binding affinities:

```python get_binding_scores.py -t mhcnuggets -n [PATIENT_ID].[TUMOR_ID].neoepiscope.comprehensive.out -d [NEOEPISCOPE_REPO]/neoepiscope/availableAlleles.pickle -a [HLA_ALLELE] -o [OUTPUT_DIR]```

`PATIENT_ID` is the patient identifier for the sample (first column of the tumor-normal pair manifest file described in "Patient genomic data" above), and `TUMOR_ID` is the tumor identifier for the sample (third column f the tumor-normal pair manifest file described in "Patient genomic data" above). `NEOEPISCOPE_REPO` is the path to the `neoepiscope` git repository. `HLA_ALLELE` is the HLA allele to use. `OUTPUT_DIR` is the path the directory where the script will write a pickled dictionary for each allele, linking each neoepitope as a key to it's binding affinity for the HLA allele used as a value (or 'NA' if binding affinity predictions for that allele/peptide are not possible).

We also used `netMHCpan` to predicting binding affinities of unphased neoepitopes:

```python get_binding_scores.py -t netMHCpan -n [PATIENT_ID].[TUMOR_ID].neoepiscope.somatic.out -d [NEOEPISCOPE_REPO]/neoepiscope/availableAlleles.pickle -a [HLA_ALLELE] -o [OUTPUT_DIR]```

Note that a different `OUTPUT_DIR` should be used in this case to separate phased vs. unphased binding affinity predictions, otherwise all inputs are the same as above.


### Compiling DNA/clinical data table

To compile the main data table summarizing DNA genomic data and clinical information, we used a [python script](scripts/immunorx_data_table.py):

```python immunorx_data_table.py -o [OUTPUT_DIR] -c [COVERAGE_SUMMARY] -s [SAMPLE_DATA] -m [MANIFEST] -v [CONSENSUS_VCF_DIR] -r [RAW_VCF_DIR] -t [HLA_TYPE_DIR] -i [MSI_DIR] -n [NEOEPITOPE_DIR] -b [BINDING_DIR] -u [UNPHASED_BINDING_DIR] -f [HLA_REFERENCE_FASTA]```

`OUTPUT_DIR` is the path to the directory to which the script will write the summary file, `immunotherapy_data_table.tsv`. `COVERAGE_SUMMARY` is the path to the output file you specified for the genomic coverage table (see "Genome coverage" above). `SAMPLE_DATA` is the path to the [sample information file](data/immunotherapy_sample_key.tsv) containing clinical information. `MANIFEST` is the path to the tumor-normal pair manifest described in "Patient genomic data" above. `CONSENSUS_VCF_DIR` is the path to the output directory containing consensus somatic VCFs you specified for the [consensus calling script](scripts/VCF_parse.py) (see "Somatic variant calling and consensus calling" above). `RAW_VCF_DIR` is the path to the directory you specified for mc3 output of raw VCFs in per-tumor subdirectories (see "Somatic variant calling and consensus calling" above).  `HLA_TYPE_DIR` is the path to the directory containing `Optitype` and `seq2hla` data (see "HLA typing" above). `MSI_DIR` is the path to the directory you specified for the `SAVEPATH` variable for mSINGS (see "MSI calculations" above). `NEOEPITOPE_DIR` is the path to the directory containing the output from neoepiscope (see "DNA neoepitope prediction" above). `BINDING_DIR` is the path to the output directory you specified for binding affinity predictions for phased epitopes (see "Binding affinity analysis" above). `UNPHASED_BINDING_DIR` is the path to the output directory you specified for binding affinity predictions for unphased epitopes (see "Binding affinity analysis" above). `HLA_REFERENCE_FASTA` is the path to the HLA reference file from `Optitype` (see the "For other data processing" subsection of "Genomic annotation data" above). 


### Processed epitope analysis

For each tumor sample, we ran netCTLpan to predict naturally processed neoepitopes of sizes 8, 9, 10, 11 for each patient MHC Class I allele using `NetCTLpan`. We used a [python script](scripts/generate_netctlpan_data.py) to facilitate this:

```python generate_netctlpan_data.py -i [NEOEPISCOPE_FILE] -o [OUTPUT_DIR] -e [NETCTLPAN] -a [ALLELE] -s [SIZE]```

`NEOEPISCOPE_FILE` is the path to the output for that sample from neoepiscope, run with phasing (the `.neoepiscope.comprehensive.out` file ending, see "DNA neoepitope prediction" above). `OUTPUT_DIR` is the path to the directory where FASTAs for input to `NetCTLpan` and the `NetCTLpan` output files will be written. `NETCTLPAN` is the path to your `NetCTLpan` executable. `ALLELE` is the HLA allele for the run, and `SIZE` is the peptide size for the run (8, 9, 10, or 11). Run for every peptide size/MHC Class I allele combination for each sample.

To process the results, we used a [python script](scripts/process_netctlpan_results.py):

```python process_netctlpan_results.py -m [MANIFEST] -i [NETCTL_DIR] -o [OUTPUT_DIR]/netctlpan_burdens.tsv```

`MANIFEST` is the path to the tumor-normal pair manifest described in "Patient genomic data" above. `NETCTL_DIR` is the `OUTPUT_DIR` from the previous step, containing the pickled dictionaries generated by the script. `OUTPUT_DIR` is the directory to which the script will write the `netctlpan_burdens.tsv` file.


### RNA-seq alignment

For use in identifying cancer-specific splicing junctions and for determining epitope expression, we aligned RNA-seq reads using `STAR` to both the GRCh37 and GRCh38 genome builds. We first created a set of `STAR` indexes from the GENCODE reference FASTA files (see the "For RNA-seq alignments" subsection of "Genomic annotation data" above) for each build:

```STAR --runMode genomeGenerate --genomeDir [INDEX_DIR] --genomeFastaFiles [REFERENCE_FASTA] --sjdbGTFfile [GTF]```

`INDEX_DIR` is the path to the directory to write the index to, `REFERENCE_FASTA` is the path to the GENCODE reference FASTA file, and `GTF` is the path to the GENCODE GTF file. Be sure to do this for both GRCh37 and GRCh38 references, storing each in a separate `INDEX_DIR`.

Then, we aligned each pair of RNA-seq FASTQ files to each reference genome:

```STAR --runMode alignReads --outSAMattributes NH HI AS nM MD --outSAMstrandField intronMotif --outFileNamePrefix [OUTPUT_DIR]/[RNA_SAMPLE_ID] --genomeDir [INDEX_DIR] --readFilesCommand zcat --readFilesIn [FASTQ1] [FASTQ2]```

`OUTPUT_DIR` is the path to the directory to which `STAR` will write the output files. `RNA_SAMPLE_ID` is the RNA sample identifier for the sample (fourth column of the tumor-normal pair manifest file described in "Patient genomic data" above). This will create a directory named for that identifier in your `OUTPUT_DIR`. `INDEX_DIR` is the path to the directory containing the `STAR` index to use. `FASTQ1` and `FASTQ2` are the paths to the forward and reverse FASTQ files for the sample. Again, be sure to do this for each reference build, specifying a different `OUTPUT_DIR` for each build.

For use in identifying cancer-specific retained introns, we also aligned RNA-seq reads using `bowtie2`. First, we generated transcript and intron reference FASTA files using a python script from `kma`: 

```python [KMA LIBRARY]/generate_introns.py --genome [REFERENCE_FASTA] --gtf [GTF] --out [OUTPUT_DIR]```

`KMA LIBRARY` is the path to the `kma` package library for `R`. The `REFERENCE_FASTA` is the path to the FASTA file in your reference data directory (see the "For the WES alignment/bam processing and somatic variant calling workflows" subsection of "Genomic Annotation Data" above). `GTF` is the path to the GENCODE GRCh37 GTF file (see the "For RNA-seq alignments" subsection of "Genomic annotation data" above). `OUTPUT_DIR` is the path to the directory to which the script will write the reference FASTA files.

Then, we combined the transcript and reference FASTA files into a single FASTA file:

```cat [OUTPUT_DIR]/trans.fa [OUTPUT DIR]/introns.fa > [OUTPUT_DIR]/trans_and_introns.fa```

`OUTPUT_DIR` is the path to the output directory from the previous step.

Next, we indexed the transcript/intron reference with `bowtie2`:

```bowtie2-build --offrate 1 [OUTPUT_DIR]/trans_and_introns.fa [OUTPUT_DIR]/trans_and_introns```

Once again, `OUTPUT_DIR` is the path to the output directory from the previous step.

Finally, we aligned the RNA-seq reads for each tumor sample, *as well as each melanocyte sample*, to the indexed reference transcript/intron FASTA:

```bowtie2 -p 12 -k 200 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 -x [OUTPUT_DIR]/trans_and_introns -1 [FASTQ1] -2 [FASTQ2] | samtools view -Sb - > [BAM_DIR]/[RNA_SAMPLE_ID]_hits.bam```

Once again, `OUTPUT_DIR` is the path to the output directory from the previous step. `FASTQ1` and `FASTQ2` are the forward and reverse FASTQ files for the sample. `BAM_DIR` is the path to the directory to which the alignments will be written. `RNA_SAMPLE_ID` is the RNA sample identifier for the sample (fourth column of the tumor-normal pair manifest file described in "Patient genomic data" above, or the melanocyte sample accession number in the melanocyte accession file).


### Tumor-specific splice junctions

We first created TCGA and GTEx junction indexes following the instructions in step 1 of the "Execute" section of [this repository](github.com/JulianneDavid/shared-cancer-splicing), using a virtual environment:

```
conda activate cravat_env

python3 jx_indexer.py -d [DB_DIR] index -c [GTEx_JUNCTION_COVERAGE] -C [TCGA_JUNCTION_COVERAGE] -b [GTEx_JUNCTION_BED] -B [TCGA_JUNCTION_BED] -p [GTEx_PHEN] -P [TCGA_PHEN] -s [RECOUNT_SAMPLE_IDS] -g [GENCODE_ANNOTATION_GTF]

conda deactivate
```

See the "Data" section of the respository for descriptions of the input files. The `DB_DIR` now contains junction indexes. (Note that this process may take around 24 hours to complete.)

Next, we called junctions present in our tumor samples but not GTEx or SRA using a [python script](scripts/cancer_junction_query.py):

```python cancer_junction_query.py -d [DB_DIR] -o [OUTPUT_DIR] -j [RNA_ALIGNMENT_DIR] -g [GENCODE_GTF] --recursive-glob --tumor-prevalences Skin_Cutaneous_Melanoma Kidney_Renal_Clear_Cell_Carcinoma```

`DB_DIR` is the path to the directory containing the junction indexes from the previous step. `OUTPUT_DIR` is the path to the directory to which the script will write output. `RNA_ALIGNMENT_DIR` is that path to your GRCh38 RNA-seq alignments (stored in per-sample subdirectories). `GENCODE_GTF` is the path to the GENCODE v28 GTF file (see the "For tumor-specific splice junction identification" subsection of "Genomic Annotation Data" above).

We then filtered out junctions that were present in any melanocyte samples in the SRA using a [python script](scripts/singleexpt_jxs_SRA_filter.py):

```python singleexpt_jxs_SRA_filter.py --snaptron-results [THIS_REPO]/data/ -o [JUNCTION_DIR] --sra-filter melanocyte_primarycell melanocyte_cellline```

`THIS_REPO` is the path to this repository, which contains necessary `Snaptron` data in the `data` subdirectory. (These data represent the results of `Snaptron` queries to the SRA for identifying melanocyte junctions.) `JUNCTION_DIR` is the path to the directory to which this script will write the filtered junctions.

Finally, to tally tumor-specific junction burden, we used a [python script](scripts/process_neojunctions.py):

```python process_neojunctions.py -m [MANIFEST] -o [OUTPUT_DIR] -j [JUNCTION_DIR]```

`MANIFEST` is the path to the tumor-normal pair manifest described in "Patient genomic data" above. `OUTPUT_DIR` is the path to the directory to which the script will write the output file, `patient_jx_burdens.tsv`. `JUNCTION_DIR` is the path to the directory containing the predicted tumor-specific junctions.


### Tumor-specific retained intron and retained intron neoepitope identification

To identify retained introns, we used `kma`. Per the tool's recommendation, we first quantified reads for each tumor sample *and each melanocyte sample* using `express`:

```mkdir [OUTPUT_DIR]/RNA_SAMPLE_ID]/```

```express -o [OUTPUT_DIR]/RNA_SAMPLE_ID]/ [TRASCRIPTOME_INTRON_REFERENCE_DIR]/trans_and_introns.fa [RNA_SAMPLE_ID]_hits.bam```

`OUTPUT_DIR` is the path to the directory you would like `express` to write output to - use a different directory for tumor samples vs. melanocyte samples. `RNA_SAMPLE_ID` is the RNA sample identifier for the sample (fourth column of the tumor-normal pair manifest file described in "Patient genomic data" above or the melanocyte sample accession number in the melanocyte accession file). `TRASCRIPTOME_INTRON_REFERENCE_DIR` is the path to the directory containing the merged transcript and intron reference FASTA (described in "RNA-seq alignment" above). 

Next, edit the [`kma` analysis R script](scripts/kma_analysis.R) at any line with a comment that says "## UPDATE THIS PATH". On line 6, update the directory to be the `OUTPUT_DIR` you used for storing express output for your tumor samples; on line 32, update the directory to be the `OUTPUT_DIR` you used for storing express output for your melanocyte samples. On line 17, update the directory to be the one containing your `intron_to_transcripts.txt` file created by `kma`'s `generate_introns.py` script (described in "RNA-seq alignment" above). On lines 27 and 51, update the directory to wherever you'd like to store information on retained intron read counts (may be the same or different directories). Then, run the script:

```Rscript kma_analysis.R```

This will generate two output files in the directory/directories you specified on lines 27 and 51, `intron_retention_read_counts.tsv` and `melanocyte_intron_retention_read_counts.tsv`. 

Next, we used a [python script](scripts/find_intron_retention_outliers.py) to find outliers in the data:

```python find_intron_retention_outliers.py -r intron_retention_read_counts.tsv -i [TRASCRIPTOME_INTRON_REFERENCE_DIR]/intron_to_transcripts.txt -o intron_retention_outliers.tsv```

```python find_intron_retention_outliers.py -r melanocyte_intron_retention_read_counts.tsv -i [TRASCRIPTOME_INTRON_REFERENCE_DIR]/intron_to_transcripts.txt -o melanocyte_intron_retention_outliers.tsv```

`TRASCRIPTOME_INTRON_REFERENCE_DIR` is the path to the directory containing the merged transcript and intron reference FASTA (described in "RNA-seq alignment" above). Adjust commands to include paths to directories containing the input and output files as relevant. 

Next, we filtered out outlier introns that are also in melanocyte data using a [python script](scripts/intron_retention_summary.py):

```python intron_retention_summary.py -m [MANIFEST] -o melanocyte_intron_retention_outliers.tsv -t intron_retention_outliers.tsv -f [FILTERED_OUTLIERS]```

The `MANIFEST` is the path to the tumor-normal pair manifest described in "Patient genomic data" above. `FILTERED_OUTLIERS` is the TSV file filtered file that the script will write. Add directory paths in front of the melanocyte and tumor outlier files as necessary. 

Finally, we used a [python script](scripts/retained_intron_epitopes.py) to summarize the burdens of retained introns and neoepitopes derived from retained introns for each patient:

```python retained_intron_epitopes.py -m [MANIFEST] -o [OUTPUT_DIR] -t [HLA_TYPE_DIR] -i [TRASCRIPTOME_INTRON_REFERENCE_DIR]/intron_to_transcripts.txt -f [FILTERED_OUTLIERS] -r [HLA_REFERENCE_FASTA]```

The `MANIFEST` is the path to the tumor-normal pair manifest described in "Patient genomic data" above. `OUTPUT_DIR` is the path to the directory to which the script will write the burden summary file. `HLA_TYPE_DIR` is the path to the directory containing `Optitype` and `seq2hla` data (see "HLA typing" above). `TRASCRIPTOME_INTRON_REFERENCE_DIR` is the path to the directory containing the merged transcript and intron reference FASTA (described in "RNA-seq alignment" above). `FILTERED_OUTLIERS` is the path to the output from the previous step. `HLA_REFERENCE_FASTA` is the path to the HLA reference file from `Optitype` (see "Genomic annotation data" above).


### Extended neoepitope burden

First, we made a blast protein database the GRCh37 protein annotation:

```makeblastdb -in [GRCH37_PROTEIN_FASTA] -input_type fasta  -title [DB_TITLE] -parse_seqids -out [DB_TITLE]```

`GRCH37_PROTEIN_FASTA` is the path to the GRCh37 protein reference FASTA file from GENCODE (see the "For RNA-seq alignments" subsection of "Genomic annotation data" above). `DB_TITLE` is the name for the database, including path information.

Next, we created a dictionary linking TCGA disease types to transcripts that are expressed in that disease type using a [python script](scripts/TCGA_expression_map.py):

```python TCGA_expression_map.py -i [TCGA_EXPRESSION_DIR]```

`TCGA_EXPRESSION_DIR` is the path to the directory containing the TPM expression data/mapping file (see the "For extended neoepitope burden analysis" subsection of "Genomic annotation data" above). The script will write a pickled dictionary, `TCGA_expression.pickle`, into that directory.

For each tumor sample, we used a [python script](scripts/extended_epitope_burden.py) to gather information about each neopeptide:

```python extended_epitope_burden.py -p [PATIENT_ID] -w [TUMOR_ID] -r [RNA_SAMPLE_ID] -n [NEOEPITOPE_DIR] -o [OUTPUT_DIR] -e [BLASTP] -b [DB_TITLE] -a [RNA_ALIGNMENTS] -t [TCGA_EXPRESSION_DIR]/TCGA_expression.pickle -s [DISEASE_SITE] -d [ALLELE_DICT] -l [LIFTOVER] -c [CHAIN_FILE] -m [BINDING_DIR]```

`PATIENT_ID` is the patient identifier for the sample (first column of the tumor-normal pair manifest file described in "Patient genomic data" above), `TUMOR_ID` is the tumor identifier for the sample (third column f the tumor-normal pair manifest file described in "Patient genomic data" above), and `RNA_SAMPLE_ID` is the RNA sample identifier for the sample (fourth column of the tumor-normal pair manifest file described in "Patient genomic data" above). `NEOEPITOPE_DIR` is the path to the directory containing the output from neoepiscope (see "DNA neoepitope prediction" above). `OUTPUT_DIR` is the path to the directory the script will write output to. `BLASTP` is the path to your `blastp` executable. `DB_TITLE` is the name for the protein database you created. `RNA_ALIGNMENTS` is the path to the directory containing your alignments of RNA-seq data to GRCh37 using `STAR` (see "RNA-seq alignment" above). `TCGA_EXPRESSION_DIR` is the path to your TCGA expression directory. `DISEASE_SITE` is the cancer type for the patient (written as melanoma, NSCLC, colon, endometrial, thyroid, prostate, or RCC). `ALLELE_DICT` is the path to the `availableAlleles.pickle` dictionary from `neoepiscope`. `LIFTOVER` is the path to your `liftOver` executable. `CHAIN_FILE` is the path to your hg19-to-hg38 chain file (see "Genomic annotation data" above). `BINDING_DIR` is the path to the output directory you specified for binding affinity predictions (see "Binding affinity analysis" above).

Finally, we processed this data for all patients using a [python script](scripts/process_extended_burden.py):

```python process_extended_burden.py -b [BURDEN_DIR] -r [RNA_ALIGNMENTS] -m [MANIFEST]```

`BURDEN_DIR` is the path to the directory you used as the `OUTPUT_DIR` in the previous step. `RNA_ALIGNMENTS` is the path to the directory containing your alignments of RNA-seq data to GRCh37 using `STAR` (see "RNA-seq alignment" above). `MANIFEST` is the path to the tumor-normal pair manifest described in "Patient genomic data" above. This script will produce two output files, `summarized_epitope_burden_updated.tsv` and `raw_epitope_table_updated.tsv`.


### Driver variant analysis

First, we annotated each consensus somatic VCF using `Open-CRAVAT`, selecting ClinVar annotations:

```conda activate cravat_env```

```cravat [SOMATIC_VCF] -a clinvar -v -l hg19```

```conda deactivate```

`SOMATIC_VCF` is the path to the consensus somatic VCF produced by the [VCF_parse.py script](scripts/VCF_parse.py) (see "Somatic variant calling and consensus calling" above). `Open-CRAVAT` will write its output into the same directory as the VCF file.

Next, we used a [python script](scripts/get_driver_mutations.py) to identify driver variants and their resulting neoepitopes:

```python get_driver_mutations.py -m [MANIFEST] -o [OUTPUT_DIR] -v [VCF_DIR] -e [NEOEPITOPE_DIR] -b [BINDING_DIR]```

`MANIFEST` is the path to the tumor-normal pair manifest described in "Patient genomic data" above. `OUTPUT_DIR` is the path to the directory to which the script will write its output files, `drivers.tsv`, `mut_binders.tsv`, `driver_mut_binders.tsv`, `ep_binders.tsv`, and `driver_ep_binders.tsv`. `VCF_DIR` is the path to the directory you used as the `OUTPUT_DIR` when producing consensus somatic VCFs using the [VCF_parse.py script](scripts/VCF_parse.py) (see "Somatic variant calling and consensus calling" above). `NEOEPITOPE_DIR` is the path to the directory containing the output from neoepiscope (see "DNA neoepitope prediction" above). `BINDING_DIR` is the path to the output directory you specified for binding affinity predictions (see "Binding affinity analysis" above).


### Processing TCGA survival data

To process TCGA mutation/clinical data for survival analysis, we used a [python script](scripts/process_TCGA_clinical.py):

```python process_TCGA_clinical.py -m [SKCM_MAF_DIR] -c [SKCM_CLINICAL_DATA] -o [OUTPUT_DIR]/skcm_tmb_survival_data.tsv```

```python process_TCGA_clinical.py -m [KIRC_MAF_DIR] -c [KIRC_CLINICAL_DATA] -o [OUTPUT_DIR]/kirc_tmb_survival_data.tsv```

`SKCM_MAF_DIR` and `KIRC_MAF_DIR` are the paths to the directories containing the SKCM and KIRC MAF files, respectively (see "Genomic annotation data" above). `SKCM_CLINICAL_DATA` and `KIRC_CLINICAL_DATA` are the `SKCM.clin.merged.txt` file from the SKCM clinical data and `KIRC.clin.merged.txt` file from the KIRC clinical data, respectively (see the "For survival analysis" subsection of "Genomic annotation data" above). `OUTPUT_DIR` is the path to the directory where output will be written.


### Statistical analysis and figure generation

To create files in a suitable format for bar plotting of variant/neoepitope burdens, we used a [python script](scripts/enumerate_variants_and_epitopes.py):

```python enumerate_variants_and_epitopes.py -i [DNA_SUMMARY_FILE] -j [JX_SUMMARY_FILE] -r [RI_SUMMARY_FILE] -o [OUTPUT_DIR]```

`DNA_SUMMARY_FILE` is the path to the `immunotherapy_data_table.tsv` file generated by the [`immunorx_data_table.py` script](scripts/immunorx_data_table.py) (see "Compiling DNA/clinical data table" above). `JX_SUMMARY_FILE` is the path to the `patient_jx_burdens.tsv` file generated by the[`process_neojunctions.py` script](scripts/process_neojunctions.py) (see "Tumor-specific splice junctions" above). `RI_SUMMARY_FILE` is the path to the `full_intron_retention_burden.tsv` file generated by the [`retained_intron_epitopes.py` script](scripts/retained_intron_epitopes.py) `OUTPUT_DIR` is the path to the directory where the output files, `rna_enumerated_mutations.tsv` and `rna_enumerated_epitopes.tsv`, will be written.

To perform statistical analysis and generate figures, we used [3](scripts/immunotherapy_analysis.R) [R](scripts/survival-plots.R) [scripts](scripts/immuno-rx-models.R). First, move or copy the following files to a working directory of your choice:

* `immunotherapy_data_table.tsv` (from [`immunorx_data_table.py`](scripts/immunorx_data_table.py), see "Compiling DNA/clinical data table" above)
* `epitope_upset.tsv` (from [`epitope_upset_combos.py`](scripts/epitope_upset_combos.py), see "DNA neoepitope prediction", above)
* `epitope_upset_counts.tsv` (from [`epitope_upset_combos.py`](scripts/epitope_upset_combos.py), see "DNA neoepitope prediction", above)
* `drivers.tsv` (from [`get_driver_mutations.py`](scripts/get_driver_mutations.py), see "Driver variant analysis" above)
* `ep_binders.tsv` (from [`get_driver_mutations.py`](scripts/get_driver_mutations.py), see "Driver variant analysis" above)
* `mut_binders.tsv` (from [`get_driver_mutations.py`](scripts/get_driver_mutations.py), see "Driver variant analysis" above)
* `driver_ep_binders.tsv` (from [`get_driver_mutations.py`](scripts/get_driver_mutations.py), see "Driver variant analysis" above)
* `driver_mut_binders.tsv` (from [`get_driver_mutations.py`](scripts/get_driver_mutations.py), see "Driver variant analysis" above)
* `summarized_epitope_burden.tsv` (from [`process_extended_burden.py`](scripts/process_extended_burden.py), see "Extended neoepitope burden" above)
* `patient_jx_burdens.tsv` (from [`process_neojunctions.py`](scripts/process_neojunctions.py), see "Tumor-specific splice junctions" above)
* `full_intron_retention_burden.tsv` (from [`retained_intron_epitopes.py`](scripts/retained_intron_epitopes.py), see "Tumor-specific retained intron and retained intron neoepitope identification" above)
* `rna_enumerated_mutations.tsv` (from [`enumerate_variants_and_epitopes.py`](scripts/enumerate_variants_and_epitopes.py), see above)
* `rna_enumerated_neoepitopes.tsv` (from [`enumerate_variants_and_epitopes.py`](scripts/enumerate_variants_and_epitopes.py), see above)
* `netctlpan_burdens.tsv` (from [`process_netctlpan_results.py`](scripts/process_netctlpan_results.py), see "Processed epitope analysis" above)
* `skcm_tmb_survival_data.tsv` (from [`process_TCGA_clinical.py`](scripts/process_TCGA_clinical.py), see "Processing TCGA survival data" above)
* `kirc_tmb_survival_data.tsv` (from [`process_TCGA_clinical.py`](scripts/process_TCGA_clinical.py), see "Processing TCGA survival data" above)

Then, run the R scripts interactively within R studio, allowing you to adjust figure dimensions and see results as desired. First, run [`immunotherapy_analysis.R`](scripts/immunotherapy_analysis.R), changing the working directory on line 11 to the directory containing the above files. Then, run [`survival-plots.R`](scripts/survival-plots.R), changing the working directory on line 1 to the directory containing the above files. Finally, run [`immuno-rx-models.R`](scripts/immuno-rx-models.R), changing the working directory on line 2 to the directory containing the above files.


## Support

For any questions or concerns regarding the analysis process, please raise an issue on this repository or email hellopdxgx@gmail.com.


## References

Graff et al. Early evidence of anti-PD-1 activity in enzalutamide-resistant prostate cancer. Oncotarget. 2016; 7:52810-52817.

Le et al. Mismatch repair deficiency predicts response of solid tumors to PD-1 blockade. Science. 2017; 357:409-413.

Wood et al. Burden of tumor mutations, neoepitopes, and other variants are dubious predictors of cancer immunotherapy response and overall survival. Preprint.


