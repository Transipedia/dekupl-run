import os
import gzip
from snakemake.utils import R

__author__ = "Jérôme Audoux (jerome.audoux@inserm.fr)"

configfile: "config.json"

# COMMON VARIABLES
TMP_FOLDER      = config['tmp_dir']+"/dekupl_tmp"
FASTQ_DIR       = config['fastq_dir']
KALLISTO_DIR    = "kallisto"
SAMPLE_NAMES    = [i['name'] for i in config["samples"]]
CONDITION_COL   = "condition"
CONDITION_A     = config['Ttest']['condition']['A']
CONDITION_B     = config['Ttest']['condition']['B']

# FILES
RAW_COUNTS                  = "raw-counts.tsv.gz"
NO_GENCODE_COUNTS           = "noGENCODE-counts.tsv.gz"
DIFF_COUNTS                 = CONDITION_A + "vs" + CONDITION_B + "-counts.tsv.gz"
MERGED_DIFF_COUNTS          = "merged-" + DIFF_COUNTS
SAMPLE_CONDITIONS           = "sample_conditions.tsv"
SAMPLE_CONDITIONS_FULL      = "sample_conditions_full.tsv"
GENCODE_FASTA               = "annotations/gencode.v24.transcripts.fa.gz"
GENCODE_COUNTS              = "annotations/gencode.v24.transcripts.tsv.gz"
TRANSCRIPT_TO_GENE_MAPPING  = "annotations/transcript_to_gene_mapping.tsv"
KALLISTO_INDEX              = "annotations/gencode.v24.transcripts-kallisto.idx"
TRANSCRIPT_COUNTS           = KALLISTO_DIR + "/transcript_counts.tsv.gz"
GENE_COUNTS                 = KALLISTO_DIR + "/gene_counts.tsv.gz"
DEGS                        = CONDITION_A + "vs" + CONDITION_B + "-DEGs.tsv"
NORMALIZATION_FACTORS       = "normalization_factors.tsv"

# Debug
#config['dekupl_counter']['min_recurrence'] = 2
#config['dekupl_counter']['min_recurrence_abundance'] = 1

# binaries
REVCOMP         = config['bin_dir'] + "/revCompFastq.pl"
DEKUPL_COUNTER  = config['bin_dir'] + "/dekupl-counter"
DIFF_FILTER     = config['bin_dir'] + "/diffFilter.pl"
TTEST_FILTER    = config['bin_dir'] + "/TtestFilter.R"
KALLISTO        = config['bin_dir'] + "/kallisto"
MERGE_COUNTS    = config['bin_dir'] + "/mergeCounts.pl"
MERGE_TAGS      = config['bin_dir'] + "/mergeTags.pl"

rule all:
  input: MERGED_DIFF_COUNTS

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
    kallisto_outputs  = expand("{kallisto_dir}/{sample}", sample = SAMPLE_NAMES, kallisto_dir = KALLISTO_DIR),
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
    size_factors = data.frame(sample = names(sizeFactors(dds)), normalization_factor = sizeFactors(dds), row.names=NULL)
    write.table(size_factors,file="{output.normalization_factors}",sep="\t",quote=FALSE, row.names = FALSE)

    # Write DEGs
    res <- results(dds, contrast=c("{CONDITION_COL}","{CONDITION_A}","{CONDITION_B}"))
    write.table(res,file="{output.differentially_expressed_genes}",sep="\t",quote=FALSE)

    #res_selection <- subset(res,padj < padj_threshold & abs(log2FoldChange) > log2foldchange_threshold)
    #write.table(res_selection,file=params$DEGS_calling_output,sep="\t",quote=FALSE)
    
    """)

###############################################################################
#
# STEP 2: KMER COUNTS
#         Compiple DEkupl counter and count k-mers on all the samples
#
# 2.1 Compile dekupl counter and create a symlink in bin directory
rule build_dekupl_counter:
  output: DEKUPL_COUNTER
  threads: 10
  run: 
    os.chdir("share/dekupl-counter")
    if not os.path.exists("build"):
      os.mkdir("build")
    os.chdir("build")
    shell("cmake -DNONCANONICAL=1 ..")
    shell("make -j {threads}")
    os.chdir("../../../bin")
    shell("ln -s ../share/dekupl-counter/build/tools/dekupl-counter .")

# 2.2 Reverse complement left mate for stranded dekupl counts
rule revcomp_pairs:
  input:  FASTQ_DIR + "/{sample}_1.fastq.gz"
  output: temp(TMP_FOLDER+"/{sample}_1.fastq.gz")
  shell:  "{REVCOMP} <(zcat {input}) | gzip -c > {output}"

# 2.3 Counts the k-mers on all the samples together
rule dekupl_counter:
  input: 
    fastq_files = expand("{tmp_folder}/{sample}_1.fastq.gz {fastq_dir}/{sample}_2.fastq.gz".split(), sample = SAMPLE_NAMES, fastq_dir = FASTQ_DIR, tmp_folder = TMP_FOLDER),
    dekupl_binary = DEKUPL_COUNTER
  output: 
    counts = RAW_COUNTS, 
    tmp_folder = temp(TMP_FOLDER),
    tmp_output = temp(TMP_FOLDER + "/counts.h5")
  threads: 10
  run:
    # Create a list of fastq files separated with comma
    fastq_list = ','.join(input.fastq_files)

    # Print nice headers
    shell("echo tag\t{SAMPLE_NAMES} | gzip -c > {output.counts}")
    shell("""{DEKUPL_COUNTER} \
          -kmer-size {config[kmer_length]} \
          -min-recurrence {config[dekupl_counter][min_recurrence]} \
          -min-recurrence-abundance {config[dekupl_counter][min_recurrence_abundance]} \
          -paired-end -out-tmp {TMP_FOLDER} -nb-cores {threads} \
          -max-memory {config[max_memory]} -max-disk {config[max_disk]} \
          -out {output.tmp_output} \
          -in {fastq_list} | gzip -c >> {output.counts}""")


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
rule gencode_counts:
  input: GENCODE_FASTA
  output: 
    counts = GENCODE_COUNTS,
    tmp_folder = temp(TMP_FOLDER),
    tmp_output = temp(TMP_FOLDER + "/counts.h5")
  threads: 10
  shell: """{DEKUPL_COUNTER} \
            -kmer-size {config[kmer_length]} \
            -min-recurrence 1 \
            -min-recurrence-abundance 1 \
            -abundance-min 1 \
            -max-memory {config[max_memory]} -max-disk {config[max_disk]} \
            -out-tmp {TMP_FOLDER} -nb-cores {threads} \
            -out {output.tmp_output} \
            -in {input} | gzip -c > {output.counts}"""

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
    sample_conditions = SAMPLE_CONDITIONS_FULL
  output: DIFF_COUNTS
  shell: "{TTEST_FILTER} {input.counts} {input.sample_conditions} {CONDITION_COL} {CONDITION_A} {CONDITION_B} {config[Ttest][pvalue_threshold]} config[Ttest][log2fc_threshold] | gzip -c > {output}"

rule merge_tags:
  input:
    DIFF_COUNTS
  output:
    MERGED_DIFF_COUNTS
  shell: "{MERGE_TAGS} {input} | gzip -c > {output}"

###############################################################################
#
# STEP 4: ANNOTATED DIFF K-MERS
#         Apply a T-test on all new k-mers to select only those that are
#         differentially expressed.
