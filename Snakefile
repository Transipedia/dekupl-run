#!/bin/python3

#######################################################################
# The MIT License
#
# Copyright (c) 2017, Jérôme Audoux (jerome.audoux@inserm.fr)
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files
# (the “Software”), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# The Software is provided “as is”, without warranty of any kind,
# express or implied, including but not limited to the warranties of
# merchantability, fitness for a particular purpose and
# noninfringement. In no event shall the authors or copyright holders
# be liable for any claim, damages or other liability, whether in an
# action of contract, tort or otherwise, arising from, out of or in
# connection with the software or the use or other dealings in the
# Software.
#######################################################################

import os
import gzip
import datetime
from sys import platform

__author__ = "Jérôme Audoux (jerome.audoux@inserm.fr)"


def getbasename(fileName):
    fileName = os.path.basename(fileName)
    *name, extension, compression = fileName.split(os.path.extsep)
    baseName = '.'.join(name)
    return(baseName)


configfile: "config.json"

# COMMON VARIABLES
SAMPLE_NAMES    = [i['name'] for i in config["samples"]]
CONDITION_COL   = "condition"
CONDITION_A     = config['diff_analysis']['condition']['A']
CONDITION_B     = config['diff_analysis']['condition']['B']
PVALUE_MAX      = config['diff_analysis']['pvalue_threshold']
LOG2FC_MIN      = config['diff_analysis']['log2fc_threshold']
MIN_REC         = config['dekupl_counter']['min_recurrence']
MIN_REC_AB      = config['dekupl_counter']['min_recurrence_abundance']
LIB_TYPE        = config['lib_type']    if 'lib_type'     in config else "rf"
R1_SUFFIX       = config['r1_suffix']   if 'r1_suffix'    in config else "_1.fastq.gz"
R2_SUFFIX       = config['r2_suffix']   if 'r2_suffix'    in config else "_2.fastq.gz"
CHUNK_SIZE      = config['chunk_size']  if 'chunk_size'   in config else 1000000
TMP_DIR         = config['tmp_dir']     if 'tmp_dir'      in config else os.getcwd()
KMER_LENGTH     = config['kmer_length'] if 'kmer_length'  in config else 31
DIFF_METHOD     = config['diff_method'] if 'diff_method'  in config else 'DESeq2'
DATA_TYPE       = config['data_type']   if 'data_type'    in config else 'RNA-Seq'
FRAG_LENGTH     = config['fragment_length'] if 'fragment_length' in config else 200
FRAG_STD_DEV    = config['fragment_standard_deviation'] if 'fragment_standard_deviation' in config else 30
OUTPUT_DIR      = config['output_dir']
FASTQ_DIR       = config['fastq_dir']

# DIRECTORIES
BIN_DIR         = workflow.basedir + "/bin"
TMP_DIR         = temp(TMP_DIR + "/dekupl_tmp")
GENE_EXP_DIR    = OUTPUT_DIR + "/gene_expression"
KALLISTO_DIR    = GENE_EXP_DIR + "/kallisto"
COUNTS_DIR      = OUTPUT_DIR + "/kmer_counts"
KMER_DE_DIR     = OUTPUT_DIR + "/" + CONDITION_A + "_vs_" + CONDITION_B + "_kmer_counts"
METADATA_DIR    = OUTPUT_DIR + "/metadata"
REFERENCE_DIR   = OUTPUT_DIR + "/references"
LOGS            = OUTPUT_DIR + "/Logs"

# FILES
RAW_COUNTS                  = COUNTS_DIR    + "/raw-counts.tsv.gz"
MASKED_COUNTS               = COUNTS_DIR    + "/masked-counts.tsv.gz"
NORMALIZATION_FACTORS       = COUNTS_DIR  + "/normalization_factors.tsv"
DIFF_COUNTS                 = KMER_DE_DIR   + "/diff-counts.tsv.gz"
PVALUE_ALL                  = KMER_DE_DIR   + "/raw_pvals.txt.gz"
MERGED_DIFF_COUNTS          = KMER_DE_DIR   + "/merged-diff-counts.tsv.gz"
ASSEMBLIES_FASTA            = KMER_DE_DIR   + "/merged-diff-counts.fa.gz"
ASSEMBLIES_BAM              = KMER_DE_DIR   + "/merged-diff-counts.bam"
SAMPLE_CONDITIONS           = METADATA_DIR  + "/sample_conditions.tsv"
SAMPLE_CONDITIONS_FULL      = METADATA_DIR  + "/sample_conditions_full.tsv"
DEFAULT_TRANSCRIPTS         = "".join(REFERENCE_DIR + "/gencode.v24.transcripts.fa.gz")
REF_TRANSCRIPT_FASTA        = config['transcript_fasta'] if 'transcript_fasta' in config else DEFAULT_TRANSCRIPTS
REF_TRANSCRIPT_COUNTS       = REFERENCE_DIR + "/" + getbasename(REF_TRANSCRIPT_FASTA) + ".tsv.gz"
TRANSCRIPT_TO_GENE_MAPPING  = config['transcript_to_gene'] if 'transcript_to_gene' in config else REFERENCE_DIR + "/transcript_to_gene_mapping.tsv"
KALLISTO_INDEX              = REFERENCE_DIR + "/" + getbasename(REF_TRANSCRIPT_FASTA) + "-kallisto.idx"
TRANSCRIPT_COUNTS           = KALLISTO_DIR  + "/transcript_counts.tsv.gz"
GENE_COUNTS                 = KALLISTO_DIR  + "/gene_counts.tsv.gz"
DEGS                        = GENE_EXP_DIR  + "/" + CONDITION_A + "vs" + CONDITION_B + "-DEGs.tsv"
CHECKING_PLOTS              = KMER_DE_DIR   + "/checking_plots.pdf"
DIST_MATRIX                 = GENE_EXP_DIR  + "/clustering_of_samples.pdf"
NORMALIZED_COUNTS           = GENE_EXP_DIR  + "/normalized_counts.tsv"
PCA_DESIGN                  = GENE_EXP_DIR  + "/pca_design.tsv"

# binaries
REVCOMP                 = BIN_DIR + "/revCompFastq.pl"
DEKUPL_COUNTER          = BIN_DIR + "/dekupl-counter"
DIFF_FILTER             = BIN_DIR + "/diffFilter.pl"
TTEST_FILTER            = BIN_DIR + "/TtestFilter"
KALLISTO                = BIN_DIR + "/kallisto"
JOIN_COUNTS             = BIN_DIR + "/joinCounts"
MERGE_COUNTS            = BIN_DIR + "/mergeCounts.pl"
MERGE_TAGS              = BIN_DIR + "/mergeTags"
COMPUTE_NF              = BIN_DIR + "/computeNF"
DESEQ2_REF_TRANSCRIPTS  = BIN_DIR + "/DESeq2_ref_transcripts.R"
JELLYFISH               = "jellyfish"
JELLYFISH_COUNT         = JELLYFISH + " count"
JELLYFISH_DUMP          = JELLYFISH + " dump"
PIGZ                    = "pigz"
ZCAT                    = "gunzip -c"
SORT                    = "sort"
JOIN                    = "join"

if platform == "darwin":
    SORT = "gsort"
    JOIN = "gjoin"

# SET MEMORY/THREAD USAGE FOR EACH RULE
MAX_MEM_KALLISTO  = 4000
MAX_MEM_JELLYFISH = 8000
MAX_MEM_SORT      = 3000

MAX_CPU           = 1000
MAX_CPU_JELLYFISH = 10
MAX_CPU_SORT      = 10

# GET THE METHOD USED FOR DETECT DE KMERS
if DIFF_METHOD == "DESeq2":
    TEST_DIFF_SCRIPT   = BIN_DIR + "/DESeq2_diff_method.R"
elif DIFF_METHOD == "Ttest":
    TEST_DIFF_SCRIPT   = BIN_DIR + "/Ttest_diff_method.R"
else:
    sys.exit("Invalid value for 'diff_method', possible choices are: 'DESeq2' and 'Ttest'")

# VERIFY LIB_TYPE VALUE
if LIB_TYPE not in ['rf', 'fr', 'unstranded', 'single']:
    sys.exit("Invalid value for 'lib_type', possible choices are: 'rf', 'rf' and 'unstranded'")

# VALIDATE sample names, because they will be used with R and cause errors if malformed
for name in SAMPLE_NAMES:
    if not re.match(r"^[a-zA-Z][1-9a-zA-Z_]*$", name):
        sys.exit("Invalid sample name '" + name + "'.\nSample names must start with at least one letter and then only letters, numbers and underscore characters are allowed")

# Print Variables in use
onstart:
    sys.stderr.write(
    """                ___  __ _                _
               /   \/__\ | ___   _ _ __ | |
              / /\ /_\ | |/ / | | | '_ \| |
             / /_///__ |   <| |_| | |_) | |
            /___,'\__/ |_|\_\\\__,_| .__/|_|
                                  |_|

    """)

    sys.stderr.write("***************** PARAMETERS ******************\n")

    sys.stderr.write("\n* General\n")
    sys.stderr.write("KMER_LENGTH                   = " + str(KMER_LENGTH) + "\n")
    sys.stderr.write("REF_TRANSCRIPT_FASTA          = " + str(REF_TRANSCRIPT_FASTA) + "\n")
    sys.stderr.write("TRANSCRIPT_TO_GENE_MAPPING    = " + str(TRANSCRIPT_TO_GENE_MAPPING) + "\n")

    sys.stderr.write("\n* K-mer counting\n")
    sys.stderr.write("MIN_REC     = " + str(MIN_REC) + "\n")
    sys.stderr.write("MIN_REC_AB  = " + str(MIN_REC_AB) + "\n")

    sys.stderr.write("\n* Diff analysis\n")
    sys.stderr.write("CONDITION_A = " + CONDITION_A + "\n")
    sys.stderr.write("CONDITION_B = " + CONDITION_B + "\n")
    sys.stderr.write("PVALUE_MAX  = " + str(PVALUE_MAX) + "\n")
    sys.stderr.write("LOG2FC_MIN  = " + str(LOG2FC_MIN) + "\n")
    sys.stderr.write("DIFF_METHOD = " + DIFF_METHOD + "\n")
    return []

if DATA_TYPE == "RNA-Seq":
    rule all:
      input: MERGED_DIFF_COUNTS, DEGS
else:
    rule all:
      input: MERGED_DIFF_COUNTS


# LOG FUNCTIONS
def current_date():
    return datetime.datetime.now().strftime("%d %B %Y, %H:%M:%S")

def start_log(log_file, rule_name):
    with open(log_file, "w") as f:
        f.write("******\n")
        f.write("start of rule " + rule_name + " : " + current_date() + "\n")

def end_log(log_file, rule_name):
    with open(log_file, "a") as f:
        f.write("\nend of rule " + rule_name + " : " + current_date() + "\n")
        f.write("******\n");

###############################################################################
#
# SOFTWARE INSTALLATION

rule compile_joinCounts:
  output: JOIN_COUNTS
  run:
    if not os.path.isfile(BIN_DIR + JOIN_COUNTS):
      shell("cd share/joinCounts && make")
      shell("ln -s -f ../share/joinCounts/joinCounts " + BIN_DIR)

rule compile_mergeTags:
  output: MERGE_TAGS
  input: "share/mergeTags/mergeTags.c"
  run:
    if not os.path.isfile(BIN_DIR + MERGE_TAGS):
      shell("cd share/mergeTags && make")
      shell("ln -s -f ../share/mergeTags/mergeTags " + BIN_DIR)


rule compile_computeNF:
  output: COMPUTE_NF
  input: "share/computeNF/computeNF.c"
  run:
    if not os.path.isfile(BIN_DIR + COMPUTE_NF):
      shell("cd share/computeNF && make")
      shell("ln -s -f ../share/computeNF/computeNF " + BIN_DIR)


rule compile_TtestFilter:
  input: "share/TtestFilter/TtestFilter.c"
  output: TTEST_FILTER
  run:
    if not os.path.isfile(BIN_DIR + TTEST_FILTER):
      shell("cd share/TtestFilter && make")
      shell("ln -s -f ../share/TtestFilter/TtestFilter " + BIN_DIR)


rule download_kallisto:
  output:
    kallisto_symlink = KALLISTO,
    kallisto_tarball = temp("share/kallisto.tar.gz")
  run:
    if not os.path.isfile(BIN_DIR + KALLISTO):
      shell("wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz -O {output.kallisto_tarball}")
      shell("tar -xzf {output.kallisto_tarball} -C share")
      shell("ln -s ../share/kallisto_linux-v0.43.0/kallisto " + KALLISTO)



###############################################################################
#
# DOWNLOAD REFERENCE FILES
#
# Download the gencode transcripts in fasta format (if no input transcriptome)
rule gencode_download:
  output: REF_TRANSCRIPT_FASTA
  shell: "wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.transcripts.fa.gz -O {output}"


###############################################################################
#
# BUILD INDEXES FROM REFERENCE FILES
#
# Create a Kallisto index of the reference transrciptome
rule kallisto_index:
  input:
    transcripts   = REF_TRANSCRIPT_FASTA,
    kallisto_bin  = KALLISTO
  resources: ram = MAX_MEM_KALLISTO
  output:
    KALLISTO_INDEX
  shell: "{KALLISTO} index -i {output} {input.transcripts}"


###############################################################################
#
# UTILS
# Creates :
#   1. A tabulated file with the sample names and conditions
#   2. A tabulated file with the sample names and normalization factors
#   3. A tabulated file with the sample names, condition and normalization factors

rule sample_conditions:
  output: SAMPLE_CONDITIONS
  run:
    with open(output[0], "w") as f:
      f.write("\t".join(["sample",CONDITION_COL]) + "\n")
      for sample in config["samples"]:
        f.write("\t".join([sample["name"],sample[CONDITION_COL]]) + "\n")

rule compute_normalization_factors:
  input:
    raw_counts = RAW_COUNTS,
    binary = COMPUTE_NF
  output:
    nf      = NORMALIZATION_FACTORS
  log: LOGS + "/compute_norm_factors.log"
  shell: "{COMPUTE_NF} {input.raw_counts} > {output.nf} 2> {log}"

rule sample_conditions_full:
  output:
    SAMPLE_CONDITIONS_FULL
  input:
    sample_conditions     = SAMPLE_CONDITIONS,
    normalization_factors  = NORMALIZATION_FACTORS
  shell: "{JOIN} --header -t $'\t' {input.sample_conditions} {input.normalization_factors} > {output}"


##############################################################################

# STEP 1: DIFFERENTIAL GENE EXPRESSION
#         Download kallisto, and quantify gene expression for aLogsll
#         the samples

# 1.3 Generic rule to quantify a sample with kallisto
rule kallisto_quantif:
  input:
    r1 = FASTQ_DIR + "/{sample}" + R1_SUFFIX,
    r2 = FASTQ_DIR + "/{sample}" + R2_SUFFIX,
    index = KALLISTO_INDEX
  resources: ram = MAX_MEM_KALLISTO
  params:
    output_dir = KALLISTO_DIR + "/{sample}",
  output:
    abundance_h5  = KALLISTO_DIR + "/{sample}/abundance.h5",
    abundance_tsv = KALLISTO_DIR + "/{sample}/abundance.tsv",
    run_info      = KALLISTO_DIR + "/{sample}/run_info.json"
  log : LOGS + "/{sample}_kallisto.log"
  threads: 1
  run:
    start_log(log[0],"kallisto_quantif")
    shell("{KALLISTO} quant -i {input.index} -o {params.output_dir} {input.r1} {input.r2} 2>> {log}")
    end_log(log[0],"kallisto_quantif")

rule kallisto_quantif_single_end:
  input:
    reads = FASTQ_DIR + "/{sample}.fastq.gz",
    index = KALLISTO_INDEX
  resources: ram = MAX_MEM_KALLISTO
  params:
    fragment_length     = FRAG_LENGTH,  #-l, --fragment-length=DOUBLE  Estimated average fragment length
    standard_deviation  = FRAG_STD_DEV, #-s, --sd=DOUBLE               Estimated standard deviation of fragment length
    output_dir          = KALLISTO_DIR + "/{sample}",
  output:
    abundance_h5  = KALLISTO_DIR + "/{sample}/abundance.h5",
    abundance_tsv = KALLISTO_DIR + "/{sample}/abundance.tsv",
    run_info      = KALLISTO_DIR + "/{sample}/run_info.json"
  log : LOGS + "/{sample}_kallisto.log"
  threads: 1
  run:
    start_log(log[0],"kallisto_quantif")
    options = "--single --fragment-length {params.fragment_length} --sd {params.standard_deviation}"
    shell("{KALLISTO} quant -i {input.index} -o {params.output_dir} " + options + " {input.reads} 2>> {log}")
    end_log(log[0],"kallisto_quantif")

# 1.4 Merge all transcripts counts from kallisto abundance files
rule transcript_counts:
  input:
    kallisto_outputs  = expand("{kallisto_dir}/{sample}/abundance.tsv", sample = SAMPLE_NAMES, kallisto_dir = KALLISTO_DIR)
  output:
    TRANSCRIPT_COUNTS
  run:
    extracted_counts  = expand("<(echo -e 'feature\t{sample}' && tail -n+2 {kallisto_dir}/{sample}/abundance.tsv | cut -f1,4)", sample = SAMPLE_NAMES, kallisto_dir = KALLISTO_DIR)
    shell("{MERGE_COUNTS} {extracted_counts} | gzip -c > {output}")

# 1.5 Create a conversion table from transcript id to gene ids
if 'transcript_to_gene' not in config:
    rule transcript_to_gene_mapping:
        input: REF_TRANSCRIPT_FASTA
        output: TRANSCRIPT_TO_GENE_MAPPING
        run:
            mapping = open(output[0], 'w')
            if(input[0].endswith('.gz')):
                opener = gzip.open
            else:
                opener = open
            with opener(input[0], 'rt') as f:
                for line in f:
                    if line[0] == ">":
                        fields = line[1:].split("|",2)
                        mapping.write("\t".join([fields[0],fields[1]]) + "\n")

# 1.6 Convert transcript counts to gene counts
rule gene_counts:
  input:
    transcript_counts = TRANSCRIPT_COUNTS,
    transcript_to_gene_mapping = TRANSCRIPT_TO_GENE_MAPPING
  output:
    GENE_COUNTS
  run:
    # Load the conversion hash
    conversion_hash = {}
    with open(input['transcript_to_gene_mapping'], 'r') as f:
      for line in f:
        transcript_id, gene_id = line.split()
        conversion_hash[transcript_id] = gene_id
    # Summarize transcript into gene counts
    gene_counts = {}
    header = ""
    with gzip.open(input['transcript_counts'], 'rt') as f:
      header = f.readline().rstrip()
      for line in f:
        counts = line.split()
        transcript_id = counts[0].split("|",1)[0]
        gene_id = conversion_hash[transcript_id]
        counts[1:] = [ float(i) for i in counts[1:] ]
        if gene_id in gene_counts:
          gene_counts[gene_id] = [ sum(i) for i in zip(gene_counts[gene_id], counts[1:]) ]
        else:
          gene_counts[gene_id] = counts[1:]
    # print Gene counts
    with gzip.open(output[0], 'wb') as f:
      f.write(bytes(header + "\n",'UTF-8'))
      for gene_id in gene_counts:
        f.write(bytes(gene_id + "\t" + "\t".join([str(int(x)) for x in gene_counts[gene_id]]) + "\n",'UTF-8'))

# 1.7 Differential expression with DESEQ2
rule differential_gene_expression:
  input:
    gene_counts = GENE_COUNTS,
    sample_conditions = SAMPLE_CONDITIONS
  params:
    condition_col = CONDITION_COL,
    condition_A = CONDITION_A,
    condition_B = CONDITION_B
  output:
    differentially_expressed_genes  = DEGS,
    dist_matrix			            = DIST_MATRIX,
    norm_counts		                    = NORMALIZED_COUNTS,
    pca_design			            = PCA_DESIGN
  log : LOGS + "/DESeq2_diff_gene_exp.log"
  shell: 
        """
        Rscript {DESEQ2_REF_TRANSCRIPTS} \
        {input.gene_counts} \
        {input.sample_conditions} \
        {params.condition_col} \
        {params.condition_A} \
        {params.condition_B} \
        {output.differentially_expressed_genes} \
        {output.dist_matrix} \
        {output.norm_counts} \
        {output.pca_design} \
        {log}
        """

###############################################################################
#
# STEP 2: KMER COUNTS
#         Compiple DEkupl counter and count k-mers on all the samples
#
rule jellyfish_count:
  input:
    r1 = FASTQ_DIR + "/{sample}" + R1_SUFFIX,
    r2 = FASTQ_DIR + "/{sample}" + R2_SUFFIX
  output: COUNTS_DIR + "/{sample}.jf"
  log:
    exec_time = LOGS + "/{sample}_jellyfishRawCounts_exec_time.log"
  threads: MAX_CPU_JELLYFISH
  resources: ram = MAX_MEM_SORT
  run:
    options = "-L 2 -m {KMER_LENGTH} -s 10000 -t {threads} -o {output} -F 2"

    start_log(log['exec_time'], "jellyfish_count (raw counts)")

    r1_pipe = "{ZCAT} {input.r1}" if input.r1.endswith(".gz") else "cat {input.r1}"
    r2_pipe = "{ZCAT} {input.r2}" if input.r2.endswith(".gz") else "cat {input.r2}"

    if LIB_TYPE == "rf":
      options += " <(%s | {REVCOMP}) <(%s)" % (r1_pipe, r2_pipe)
      shell("echo -e \"R1 is rev comp\n\" >>{log.exec_time}")
    elif LIB_TYPE == "fr":
      options += " <(%s) <(%s | {REVCOMP})" % (r1_pipe, r2_pipe)
      shell("echo -e \"R2 is rev comp\n\" >>{log.exec_time}")
    elif LIB_TYPE == "unstranded":
      options += " -C <(%s) <(%s)" % (r1_pipe, r2_pipe)
    else:
      sys.exit('Unknown library type')

    shell("{JELLYFISH_COUNT} " + options)
    end_log(log['exec_time'], "jellyfish_count")

rule jellyfish_count_single_end:
  input:
    reads = FASTQ_DIR + "/{sample}.fastq.gz"
  output: COUNTS_DIR + "/{sample}.jf"
  log:
    exec_time = LOGS + "/{sample}_jellyfishRawCounts_exec_time.log"
  threads: MAX_CPU_JELLYFISH
  resources: ram = MAX_MEM_SORT
  run:
    options = "-L 2 -m {KMER_LENGTH} -s 10000 -t {threads} -o {output} -F 2"

    if LIB_TYPE == "unstranded":
      options += " -C"

    start_log(log['exec_time'], "jellyfish_count (raw counts)")
    shell("{JELLYFISH_COUNT} " + options + " <({ZCAT} {input.reads})")
    end_log(log['exec_time'], "jellyfish_count")

rule jellyfish_dump:
  input: COUNTS_DIR + "/{sample}.jf"
  output: COUNTS_DIR + "/{sample}.txt.gz"
  threads: MAX_CPU_SORT
  resources: ram = MAX_MEM_SORT
  log :
    exec_time = LOGS + "/{sample}_jellyfishDumpRawCounts_exec_time.log"
  run:
    start_log(log['exec_time'], "jellyfish_dump")
    shell("{JELLYFISH_DUMP} -c {input} | {SORT} -k 1 -S {resources.ram}M --parallel {threads}| pigz -p {threads} -c > {output}")
    end_log(log['exec_time'], "jellyfish_dump")

rule join_counts:
  input:
    fastq_files = expand("{counts_dir}/{sample}.txt.gz",counts_dir=COUNTS_DIR,sample=SAMPLE_NAMES),
    binary = JOIN_COUNTS
  params:
    sample_names = "\t".join(SAMPLE_NAMES)
  output: RAW_COUNTS
  log :
    exec_time = LOGS + "/joinRawCounts_exec_time.log"
  run:
    shell("echo 'tag\t{params.sample_names}' | gzip -c > {output}")
    start_log(log['exec_time'], "join_counts")
    shell("{JOIN_COUNTS} -r {MIN_REC} -a {MIN_REC_AB} {input.fastq_files} | gzip -c >> {output}")
    end_log(log['exec_time'], "join_counts")

###############################################################################
#
# STEP 3: FILTER-OUT KNOWN K-MERS
#         Default: Download gencode transcripts set and remove the k-mer occuring this
#         set from the one found in the experimental data
#

# 3.2 Counts k-mer of all transcript (for further filtration)
rule ref_transcript_count:
  input: REF_TRANSCRIPT_FASTA
  output: temp(REF_TRANSCRIPT_FASTA + ".jf")
  threads: MAX_CPU_JELLYFISH
  resources: ram = MAX_MEM_JELLYFISH
  run:
    options = "-m {KMER_LENGTH} -s 10000 -t {threads} -o {output}"
    if LIB_TYPE == "unstranded":
      options += " -C"
    if(input[0].endswith('.gz')):
      shell("{JELLYFISH_COUNT} " + options + " <({ZCAT} {input})")
    else:
      shell("{JELLYFISH_COUNT} " + options + " {input}")

rule ref_transcript_dump:
  input: REF_TRANSCRIPT_FASTA + ".jf"
  output: REF_TRANSCRIPT_COUNTS
  log :
    exec_time = LOGS + "/jellyfishDumpRefTrancriptCounts_exec_time.log"
  threads: MAX_CPU_SORT
  resources: ram = MAX_MEM_SORT
  run:
    start_log(log['exec_time'], "ref_transcript_dump")
    shell("{JELLYFISH_DUMP} -c {input} | {SORT} -k 1 -S {resources.ram}M --parallel {threads}| pigz -p {threads} -c > {output}")
    end_log(log['exec_time'], "ref_transcript_dump")

# 3.3 Filter counter k-mer that are present in the transcriptome set
rule filter_transcript_counts:
  input:
    counts = RAW_COUNTS,
    ref_transcript_counts = REF_TRANSCRIPT_COUNTS
  output: MASKED_COUNTS
  log:
    exec_time = LOGS + "/filter_transcript_counts_exec_time.log"
  run:
    start_log(log['exec_time'], "filter_transcript_counts")
    shell("{DIFF_FILTER} {input.ref_transcript_counts} {input.counts} 2>> {log.exec_time} | gzip -c > {output}")
    end_log(log['exec_time'], "filter_transcript_counts")

###############################################################################
#
# STEP 4: SELECT DIFFERENTIALLY EXPRESSED K-MERS
#         Apply a T-test on all new k-mers to select only those that are
#         differentially expressed.
#
rule test_diff_counts:
  input:
    counts = MASKED_COUNTS if DATA_TYPE == "RNA-Seq" else RAW_COUNTS,
    sample_conditions = SAMPLE_CONDITIONS_FULL,
    binary = TTEST_FILTER # this is just here to compile T-test. This rule also includes DESeq2 and Poisson tests
  output:
    diff_counts = DIFF_COUNTS,
    pvalue_all  = PVALUE_ALL,
    #tmp_dir     = TMP_DIR + "/test_diff"
    #tmp_dir     = temp(TMP_DIR + "/test_diff")
  params:
    conditionA  = CONDITION_A,
    conditionB  = CONDITION_B,
    pvalue_threshold = PVALUE_MAX,
    log2fc_threshold = LOG2FC_MIN,
    chunk_size = CHUNK_SIZE,
    tmp_dir = TMP_DIR + "/test_diff",
  threads: MAX_CPU
  log: LOGS + "/test_diff_counts.logs"
  shell: 
        """
        Rscript {TEST_DIFF_SCRIPT} \
        {input.binary} \
        {input.counts} \
        {input.sample_conditions} \
        {params.pvalue_threshold} \
        {params.log2fc_threshold} \
        {params.conditionA} \
        {params.conditionB} \
        {threads} \
        {params.chunk_size} \
        {params.tmp_dir} \
        {output.diff_counts} \
        {output.pvalue_all} \
        {log}
        """

rule merge_tags:
  input:
    counts = DIFF_COUNTS,
    binary = MERGE_TAGS
  output:
    MERGED_DIFF_COUNTS
  log:
    exec_time = LOGS + "/merge_tags_exec_time.log"
  run:
    options = "-k {KMER_LENGTH}"

    if LIB_TYPE == "unstranded":
      options += " -n"

    start_log(log['exec_time'], "merge_tags")
    shell("{MERGE_TAGS} " + options + " {input.counts} | gzip -c > {output}")
    end_log(log['exec_time'], "merge_tags")
