import os
import gzip
from snakemake.utils import R

__author__ = "Jérôme Audoux (jerome.audoux@inserm.fr)"

configfile: "config.json"

# COMMON VARIABLES
SAMPLE_NAMES    = [i['name'] for i in config["samples"]]
CONDITION_COL   = "condition"
CONDITION_A     = config['Ttest']['condition']['A']
CONDITION_B     = config['Ttest']['condition']['B']

# DIRECTORIES
BIN_DIR         = "bin"
TMP_DIR         = config['tmp_dir']+"/dekupl_tmp"
FASTQ_DIR       = config['fastq_dir']
GENE_EXP_DIR    = "gene_expression"
KALLISTO_DIR    = GENE_EXP_DIR + "/kallisto"
COUNTS_DIR      = "kmer_counts"
KMER_DE_DIR     = CONDITION_A + "_vs_" + CONDITION_B + "_kmer_counts"
METADATA_DIR    = "metadata"
REFERENCE_DIR   = "references"

# FILES
RAW_COUNTS                  = COUNTS_DIR    + "/raw-counts.tsv.gz"
NO_GENCODE_COUNTS           = COUNTS_DIR    + "/noGENCODE-counts.tsv.gz"
DIFF_COUNTS                 = KMER_DE_DIR   + "/diff-counts.tsv.gz"
MERGED_DIFF_COUNTS          = KMER_DE_DIR   + "/merged-diff-counts.tsv.gz"
SAMPLE_CONDITIONS           = METADATA_DIR  + "/sample_conditions.tsv"
SAMPLE_CONDITIONS_FULL      = METADATA_DIR  + "/sample_conditions_full.tsv"
GENCODE_FASTA               = REFERENCE_DIR + "/gencode.v24.transcripts.fa.gz"
GENCODE_COUNTS              = REFERENCE_DIR + "/gencode.v24.transcripts.tsv.gz"
TRANSCRIPT_TO_GENE_MAPPING  = REFERENCE_DIR + "/transcript_to_gene_mapping.tsv"
KALLISTO_INDEX              = REFERENCE_DIR + "/gencode.v24.transcripts-kallisto.idx"
TRANSCRIPT_COUNTS           = KALLISTO_DIR  + "/transcript_counts.tsv.gz"
GENE_COUNTS                 = KALLISTO_DIR  + "/gene_counts.tsv.gz"
DEGS                        = GENE_EXP_DIR  + "/" + CONDITION_A + "vs" + CONDITION_B + "-DEGs.tsv"
NORMALIZATION_FACTORS       = GENE_EXP_DIR  + "/normalization_factors.tsv"

# binaries
REVCOMP         = BIN_DIR + "/revCompFastq.pl"
DEKUPL_COUNTER  = BIN_DIR + "/dekupl-counter"
DIFF_FILTER     = BIN_DIR + "/diffFilter.pl"
TTEST_FILTER    = BIN_DIR + "/TtestFilter"
KALLISTO        = BIN_DIR + "/kallisto"
JOIN_COUNTS     = BIN_DIR + "/joinCounts"
MERGE_COUNTS    = BIN_DIR + "/mergeCounts.pl"
MERGE_TAGS      = BIN_DIR + "/mergeTags"
JELLYFISH       = "jellyfish"
JELLYFISH_COUNT = JELLYFISH + " count"
JELLYFISH_DUMP  = JELLYFISH + " dump"

rule all:
  input: MERGED_DIFF_COUNTS

###############################################################################
#
# COMPILATION
#
rule compile_joinCounts:
  output: JOIN_COUNTS
  run:
    shell("cd share/joinCounts && make")
    shell("ln -s -f ../share/joinCounts/joinCounts bin/")

rule compile_mergeCounts:
  output: MERGE_TAGS
  run:
    shell("cd share/mergeTags && make")
    shell("ln -s -f ../share/mergeTags/mergeTags bin/")

rule compile_TtestFilter:
  input: "share/TtestFilter/TtestFilter.c"
  output: TTEST_FILTER
  run:
    shell("cd share/TtestFilter && make")
    shell("ln -s -f ../share/TtestFilter/TtestFilter bin/")

###############################################################################
#
# UTILS
#
# Create a tabulated file with the sample name and conditions
rule sample_conditions:
  output: SAMPLE_CONDITIONS
  run:
    with open(output[0], "w") as f:
      f.write("\t".join(["sample",CONDITION_COL]) + "\n")
      for sample in config["samples"]:
        f.write("\t".join([sample["name"],sample[CONDITION_COL]]) + "\n")

rule sample_conditions_full:
  output: 
    SAMPLE_CONDITIONS_FULL
  input:
    sample_conditions     = SAMPLE_CONDITIONS,
    normalization_factors  = NORMALIZATION_FACTORS
  shell: "join --header {input.sample_conditions} {input.normalization_factors} > {output}"


###############################################################################
#
# STEP 1: DIFFERENTIAL GENE EXPRESSION
#         Download kallisto, and quantify gene expression for all
#         the samples
#
# 1.1 Download a pre-compiled version of Kallisto from github
rule download_kallisto:
  output: 
    kallisto_symlink = KALLISTO,
    kallisto_tarball = temp("share/kallisto.tar.gz")
  run:
    shell("wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz -O {output.kallisto_tarball}")
    shell("tar -xzf {output.kallisto_tarball} -C share")
    shell("ln -s ../share/kallisto_linux-v0.43.0/kallisto bin/kallisto")

# 1.2 Create a Kallisto index of the reference transrciptome
rule kallisto_index:
  input: 
    transcripts   = GENCODE_FASTA,
    kallisto_bin  = KALLISTO,
  output: 
    KALLISTO_INDEX
  shell: "{KALLISTO} index -i {output} {input.transcripts}"

# 1.3 Generic rule to quantify a sample with kallisto
rule kallito_quantif:
  input:
    r1 = FASTQ_DIR + "/{sample}_1.fastq.gz",
    r2 = FASTQ_DIR + "/{sample}_2.fastq.gz",
    index = KALLISTO_INDEX
  output:
    KALLISTO_DIR + "/{sample}"
  threads: 1
  shell: "{KALLISTO} quant -i {input.index} -o {output} {input.r1} {input.r2}"

# 1.4 Merge all transcripts counts from kallisto abundance files
rule transcript_counts:
  input: 
    kallisto_outputs  = expand("{kallisto_dir}/{sample}", sample = SAMPLE_NAMES, kallisto_dir = KALLISTO_DIR)
  output:
    TRANSCRIPT_COUNTS
  run: 
    extracted_counts  = expand("<(echo -e 'feature\t{sample}' && tail -n+2 {kallisto_dir}/{sample}/abundance.tsv | cut -f1,4)", sample = SAMPLE_NAMES, kallisto_dir = KALLISTO_DIR)
    shell("{MERGE_COUNTS} {extracted_counts} | gzip -c > {output}")

# 1.5 Create a conversion table from transcript id to gene ids
rule transcript_to_gene_mapping:
  input: GENCODE_FASTA
  output: TRANSCRIPT_TO_GENE_MAPPING
  run:
    mapping = open(output[0], 'w') 
    with gzip.open(input[0], 'rt') as f:
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
        transcript_id, trail = counts[0].split("|",1)
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
  output:
    differentially_expressed_genes = DEGS,
    normalization_factors          = NORMALIZATION_FACTORS
  run:
    R("""
    library(DESeq2) 

    # Load counts data
    countsData = read.table("{input.gene_counts}",
                            header=T,
                            row.names=1)

    # Load col data with sample specifications
    colData = read.table("{input.sample_conditions}",
                         header=T,
                         row.names=1)
    
    write(colnames(countsData),stderr())
    write(rownames(colData),stderr())

    colData = colData[colnames(countsData),,drop=FALSE]

    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData=countsData,
                                  colData=colData,
                                  design = ~ {CONDITION_COL})
    dds <- DESeq(dds)

    # Write normalization factors
    size_factors = data.frame(sample = names(sizeFactors(dds)), 
                              normalization_factor = sizeFactors(dds),
                              row.names=NULL)
    write.table(size_factors,
                file="{output.normalization_factors}",
                sep="\t",quote=FALSE, row.names = FALSE)

    write(resultsNames(dds),stderr())

    # Write DEGs
    res <- results(dds, contrast = c("{CONDITION_COL}","{CONDITION_A}","{CONDITION_B}"))
    write.table(res,file="{output.differentially_expressed_genes}",sep="\t",quote=FALSE)
    """)

###############################################################################
#
# STEP 2: KMER COUNTS
#         Compiple DEkupl counter and count k-mers on all the samples
#
rule jellyfish_count:
  input: 
    r1 = FASTQ_DIR  + "/{sample}_1.fastq.gz", 
    r2 = FASTQ_DIR  + "/{sample}_2.fastq.gz"
  output: COUNTS_DIR + "/{sample}.jf"
  threads: 10
  shell: "{JELLYFISH_COUNT} -L 2 -m {config[kmer_length]} -s 10000 -t {threads} -o {output} <(zcat {input.r1} | {REVCOMP}) <(zcat {input.r2})"

rule jellyfish_dump:
  input: COUNTS_DIR + "/{sample}.jf"
  output: COUNTS_DIR + "/{sample}.txt.gz"
  threads: 10
  resources: ram=10
  shell: "{JELLYFISH_DUMP} -c {input} | sort -k 1 -S {resources.ram}G --parallel {threads}| pigz -p {threads} -c > {output}"

rule join_counts:
  input: 
    fastq_files = expand("{counts_dir}/{sample}.txt.gz",counts_dir=COUNTS_DIR,sample=SAMPLE_NAMES),
    binary = JOIN_COUNTS
  params:
    sample_names = "\t".join(SAMPLE_NAMES)
  output: RAW_COUNTS
  run:
    shell("echo 'tag\t{params.sample_names}' | gzip -c > {output}")
    shell("""{JOIN_COUNTS} -r {config[dekupl_counter][min_recurrence]} \
          -a {config[dekupl_counter][min_recurrence_abundance]} \
          {input.fastq_files} | gzip -c >> {output}""")

###############################################################################
#
# STEP 3: FILTER-OUT KNOWN K-MERS
#         Download gencode transcripts set and remove the k-mer occuring this
#         set from the one found in the experimental data
#
# 3.1 Download the gencode transcripts in fasta format
rule gencode_download:
  output: GENCODE_FASTA
  shell: "wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.transcripts.fa.gz -O {output}"

# 3.2 Counts k-mer of all gencode transcript (for further filtration)
rule gencode_count:
  input: GENCODE_FASTA
  output: temp(GENCODE_FASTA + ".jf")
  threads: 10
  shell: """{JELLYFISH_COUNT} -m {config[kmer_length]} \
            -s 10000 -t {threads} -o {output} <(zcat {input})"""

rule gencode_dump:
  input: GENCODE_FASTA + ".jf"
  output: GENCODE_COUNTS
  threads: 10
  resources: ram=4
  shell: "{JELLYFISH_DUMP} -c {input} | sort -k 1 -S {resources.ram}G --parallel {threads}| pigz -p {threads} -c > {output}"

# 3.3 Filter counter k-mer that are present in the gencode set
rule filter_gencode_counts:
  input:
    counts = RAW_COUNTS,
    gencode_counts = GENCODE_COUNTS
  output: NO_GENCODE_COUNTS
  shell: "{DIFF_FILTER} {input.gencode_counts} {input.counts} | gzip -c > {output}"

###############################################################################
#
# STEP 4: SELECT DIFFERENTIALLY EXPRESSED K-MERS
#         Apply a T-test on all new k-mers to select only those that are
#         differentially expressed.
#
rule test_diff_counts:
  input: 
    counts = NO_GENCODE_COUNTS,
    sample_conditions = SAMPLE_CONDITIONS_FULL,
    binary = TTEST_FILTER
  output: DIFF_COUNTS
  shell: """{TTEST_FILTER} \
            -p {config[Ttest][pvalue_threshold]} \
            -f {config[Ttest][log2fc_threshold]} \
            {input.counts} {input.sample_conditions} \
            {CONDITION_A} {CONDITION_B} | gzip -c > {output}"""

rule merge_tags:
  input:
    counts = DIFF_COUNTS,
    binary = MERGE_TAGS
  output:
    MERGED_DIFF_COUNTS
  shell: "{MERGE_TAGS} -k {config[kmer_length]} {input} | gzip -c > {output}"

###############################################################################
#
# STEP 4: ANNOTATED DIFF K-MERS
#         Apply a T-test on all new k-mers to select only those that are
#         differentially expressed.
