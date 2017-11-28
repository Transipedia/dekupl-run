#!/bin/python3

#######################################################################
# This file is part of Dekupl. 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Jérôme Audoux (jerome.audoux@inserm.fr), Copyright (C) 2017  
#######################################################################
	
import os
import gzip
from snakemake.utils import R

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
OUTPUT_DIR      = config['output_dir']
FASTQ_DIR       = config['fastq_dir']

# DIRECTORIES
BIN_DIR         = "bin"
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
default_file                = "".join(REFERENCE_DIR + "/gencode.v24.transcripts.fa.gz")
REF_TRANSCRIPT_FASTA        = config['transcript_fasta'] if 'transcript_fasta' in config else default_file
REF_TRANSCRIPT_COUNTS       = REFERENCE_DIR + "/" + getbasename(REF_TRANSCRIPT_FASTA) + ".tsv.gz"
TRANSCRIPT_TO_GENE_MAPPING  = REFERENCE_DIR + "/transcript_to_gene_mapping.tsv"
KALLISTO_INDEX              = REFERENCE_DIR + "/" + getbasename(REF_TRANSCRIPT_FASTA) + "-kallisto.idx"
TRANSCRIPT_COUNTS           = KALLISTO_DIR  + "/transcript_counts.tsv.gz"
GENE_COUNTS                 = KALLISTO_DIR  + "/gene_counts.tsv.gz"
DEGS                        = GENE_EXP_DIR  + "/" + CONDITION_A + "vs" + CONDITION_B + "-DEGs.tsv"
CHECKING_PLOTS              = KMER_DE_DIR   + "/checking_plots.pdf"
DIST_MATRIX                 = GENE_EXP_DIR  + "/clustering_of_samples.pdf"
NORMALIZED_COUNTS           = GENE_EXP_DIR  + "/normalized_counts.tsv"
PCA_DESIGN                  = GENE_EXP_DIR  + "/pca_design.tsv"

# binaries
REVCOMP         = BIN_DIR + "/revCompFastq.pl"
DEKUPL_COUNTER  = BIN_DIR + "/dekupl-counter"
DIFF_FILTER     = BIN_DIR + "/diffFilter.pl"
TTEST_FILTER    = BIN_DIR + "/TtestFilter"
KALLISTO        = BIN_DIR + "/kallisto"
JOIN_COUNTS     = BIN_DIR + "/joinCounts"
MERGE_COUNTS    = BIN_DIR + "/mergeCounts.pl"
MERGE_TAGS      = BIN_DIR + "/mergeTags"
COMPUTE_NF      = BIN_DIR + "/compute_norm_factors.R"
JELLYFISH       = "jellyfish"
JELLYFISH_COUNT = JELLYFISH + " count"
JELLYFISH_DUMP  = JELLYFISH + " dump"
PIGZ            = "pigz"
ZCAT            = "gunzip -c"
SORT            = "sort"

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
if LIB_TYPE not in ['rf', 'fr', 'unstranded']:
    sys.exit("Invalid value for 'lib_type', possible choices are: 'rf', 'rf' and 'unstranded'")

rule all:
  input: MERGED_DIFF_COUNTS, DEGS


###############################################################################
#
# SOFTWARE INSTALLATION
#
rule compile_joinCounts:
  output: JOIN_COUNTS
  run:
    shell("cd share/joinCounts && make")
    shell("ln -s -f ../share/joinCounts/joinCounts bin/")

rule compile_mergeTags:
  output: MERGE_TAGS
  input: "share/mergeTags/mergeTags.c"
  run:
    shell("cd share/mergeTags && make")
    shell("ln -s -f ../share/mergeTags/mergeTags bin/")

rule compile_TtestFilter:
  input: "share/TtestFilter/TtestFilter.c"
  output: TTEST_FILTER
  run:
    shell("cd share/TtestFilter && make")
    shell("ln -s -f ../share/TtestFilter/TtestFilter bin/")

rule download_kallisto:
  output:
    kallisto_symlink = KALLISTO,
    kallisto_tarball = temp("share/kallisto.tar.gz")
  run:
    shell("wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz -O {output.kallisto_tarball}")
    shell("tar -xzf {output.kallisto_tarball} -C share")
    shell("ln -s ../share/kallisto_linux-v0.43.0/kallisto bin/kallisto")


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
    raw_counts = RAW_COUNTS
  output: 
    nf      = NORMALIZATION_FACTORS,
    tmp_dir = temp(TMP_DIR + "/NF")
  log: LOGS + "/compute_norm_factors.log"
  script: COMPUTE_NF

rule sample_conditions_full:
  output:
    SAMPLE_CONDITIONS_FULL
  input:
    sample_conditions     = SAMPLE_CONDITIONS,
    normalization_factors  = NORMALIZATION_FACTORS
  shell: "join --header {input.sample_conditions} {input.normalization_factors} > {output}"


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
  output:
    dir           = KALLISTO_DIR + "/{sample}",
    abundance_h5  = KALLISTO_DIR + "/{sample}/abundance.h5",
    abundance_tsv = KALLISTO_DIR + "/{sample}/abundance.tsv",
    run_info      = KALLISTO_DIR + "/{sample}/run_info.json"
  log : LOGS + "/{sample}_kallisto.log"
  threads: 1
  shell: """

         echo -e \"******\" >{log}
         echo -e \"start of rule kallisto_quantif : $(date)\n\" >>{log}

         {KALLISTO} quant -i {input.index} -o {output.dir} {input.r1} {input.r2} 2>> {log}

         echo -e \"\nend of rule kallisto_quantif : $(date)\n\" >>{log}
         echo -e \"******\" >>{log}

         """

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
  input: REF_TRANSCRIPT_FASTA
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
    differentially_expressed_genes  = DEGS,
    dist_matrix			            = DIST_MATRIX,
    norm_counts		                    = NORMALIZED_COUNTS,
    pca_design			            = PCA_DESIGN
  log : LOGS + "/DESeq2_diff_gene_exp.log"
  run:
    R("""
    library(DESeq2)
    library(RColorBrewer)
    library(pheatmap)
    library(ggplot2)

    write(date(),file="{log}")

    # Load counts data
    countsData = read.table("{input.gene_counts}",
                            header=T,
                            row.names=1,
			    check.names=FALSE)

    # Load col data with sample specifications
    colData = read.table("{input.sample_conditions}",
                         header=T,
                         row.names=1,
			 check.names=FALSE)

    write(colnames(countsData),stderr())
    write(rownames(colData),stderr())

    colData = colData[colnames(countsData),,drop=FALSE]

    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData=countsData,
                                  colData=colData,
                                  design = ~ {CONDITION_COL})
    dds <- DESeq(dds)

    #normalized counts
    NormCount<- as.data.frame(counts(dds, normalized=TRUE ))

    #writing in a file normalized counts
    normalized_counts<-data.frame(id=row.names(NormCount),NormCount,row.names=NULL)
    write.table(normalized_counts,file="{output.norm_counts}", sep="\t",row.names=F, col.names=T, quote=F)

    write(resultsNames(dds),stderr())

    # Write DEGs
    res <- results(dds, contrast = c("{CONDITION_COL}","{CONDITION_A}","{CONDITION_B}"))
    write.table(res,file="{output.differentially_expressed_genes}",sep="\t",quote=FALSE)

    rld<-rlog(dds)
    sampleDists<-dist(t(assay(rld) ) )
    sampleDistMatrix<-as.matrix( sampleDists )
    rownames(sampleDistMatrix)<-colnames(rld)
    colnames(sampleDistMatrix)<-colnames(rld)
    colours=colorRampPalette(rev(brewer.pal(9,"Blues")) )(255)

    pdf("{output.dist_matrix}",width=15,height=10)

    pheatmap(sampleDistMatrix,
	main="clustering of samples",
	clustering_distance_rows=sampleDists,
	clustering_distance_cols=sampleDists,
	col=colours,
	fontsize = 14)

    data <- plotPCA(rld,ntop=nrow(rld),returnData=TRUE)
    write.table(data,"{output.pca_design}",row.names=F, col.names=T, quote=F,sep="\t")

    print(ggplot(data,aes(PC1,PC2,color=condition))+geom_point()+geom_text(aes(label=name),hjust=0,vjust=0))

    dev.off()

    write(date(),file="{log}",append=T)
    """)

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

    shell("echo -e \"******\" >{log.exec_time}")
    shell("echo -e \"start of rule jellyfish_count (raw counts) : $(date)\n\" >>{log.exec_time}")

    if LIB_TYPE == "rf":
      options += " <({ZCAT} {input.r1} | {REVCOMP}) <({ZCAT} {input.r2})"
      shell("echo -e \"R1 is rev comp\n\" >>{log.exec_time}")
    elif LIB_TYPE == "fr":
      options += " <({ZCAT} {input.r1}) <({ZCAT} {input.r2} | {REVCOMP})"
      shell("echo -e \"R2 is rev comp\n\" >>{log.exec_time}")
    elif LIB_TYPE == "unstranded":
      options += " -C <({ZCAT} {input.r1}) <({ZCAT} {input.r2})"
    else:
      sys.exit('Unknown library type')

    shell("{JELLYFISH_COUNT} " + options)
    shell("echo -e \"\nend of rule jellyfish_count : $(date)\n\" >>{log.exec_time}")
    shell("echo -e \"******\" >>{log.exec_time}")

rule jellyfish_dump:
  input: COUNTS_DIR + "/{sample}.jf"
  output: COUNTS_DIR + "/{sample}.txt.gz"
  threads: MAX_CPU_SORT
  resources: ram = MAX_MEM_SORT
  log :
    exec_time = LOGS + "/{sample}_jellyfishDumpRawCounts_exec_time.log"
  shell: """

         echo -e \"******\" >{log.exec_time}
         echo -e \"start of rule jellyfish_dump : $(date)\n\" >>{log.exec_time}

         {JELLYFISH_DUMP} -c {input} | {SORT} -k 1 -S {resources.ram}G --parallel {threads}| pigz -p {threads} -c > {output}

         echo -e \"\nend of rule jellyfish_dump : $(date)\n\" >>{log.exec_time}
         echo -e \"******\" >>{log.exec_time}

         """

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
    shell("""

           echo -e \"******\" >{log.exec_time}
           echo -e \"start of rule join_counts : $(date)\n\" >>{log.exec_time}

           {JOIN_COUNTS} -r {MIN_REC} -a {MIN_REC_AB} \
          {input.fastq_files} | gzip -c >> {output}

          echo -e \"\nend of rule dekupl_counter : $(date)\n\" >>{log.exec_time}
          echo -e \"******\" >>{log.exec_time}

          """)

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
    shell("{JELLYFISH_COUNT} " + options + " <({ZCAT} {input})")

rule ref_transcript_dump:
  input: REF_TRANSCRIPT_FASTA + ".jf"
  output: REF_TRANSCRIPT_COUNTS
  log :
    exec_time = LOGS + "/jellyfishDumpRefTrancriptCounts_exec_time.log"
  threads: MAX_CPU_SORT
  resources: ram = MAX_MEM_SORT
  shell: """

         echo -e \"******\" >{log.exec_time}
         echo -e \"start of ref_transcript_dump : $(date)\n\" >>{log.exec_time}


        {JELLYFISH_DUMP} -c {input} | {SORT} -k 1 -S {resources.ram}G --parallel {threads}| pigz -p {threads} -c > {output}

        echo -e \"\nend of rule ref_transcript_dump : $(date)\n\" >>{log.exec_time}
        echo -e \"******\" >>{log.exec_time}

        """

# 3.3 Filter counter k-mer that are present in the transcriptome set
rule filter_transcript_counts:
  input:
    counts = RAW_COUNTS,
    ref_transcript_counts = REF_TRANSCRIPT_COUNTS
  output: MASKED_COUNTS
  log:
    exec_time = LOGS + "/filter_transcript_counts_exec_time.log"
  shell: """
         echo -e \"******\" >{log.exec_time}
         echo -e \"start of filter_transcript_counts : $(date)\n\" >>{log.exec_time}

         {DIFF_FILTER} {input.ref_transcript_counts} {input.counts} | gzip -c > {output}


         echo -e \"\nend of filter_transcript_counts : $(date)\n\" >>{log.exec_time}
         echo -e \"******\" >>{log.exec_time}

         """

###############################################################################
#
# STEP 4: SELECT DIFFERENTIALLY EXPRESSED K-MERS
#         Apply a T-test on all new k-mers to select only those that are
#         differentially expressed.
#
rule test_diff_counts:
  input:
    counts = MASKED_COUNTS,
    sample_conditions = SAMPLE_CONDITIONS_FULL,
    binary = TTEST_FILTER
  output:
    diff_counts = DIFF_COUNTS,
    pvalue_all  = PVALUE_ALL,
    tmp_dir     = TMP_DIR + "/test_diff"
    #tmp_dir     = temp(TMP_DIR + "/test_diff")
  params:
    conditionA  = CONDITION_A,
    conditionB  = CONDITION_B,
    pvalue_threshold = PVALUE_MAX,
    log2fc_threshold = LOG2FC_MIN,
    chunk_size = CHUNK_SIZE,
  threads: MAX_CPU
  log: LOGS + "/test_diff_counts.logs"
  script: TEST_DIFF_SCRIPT

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

    shell("echo -e \"******\" >{log.exec_time}")
    shell("echo -e \"start of merge_tags : $(date)\n\" >>{log.exec_time}")
    shell("{MERGE_TAGS} " + options + " {input.counts} | gzip -c > {output}")
    shell("echo -e \"\nend of merge_tags : $(date)\n\" >>{log.exec_time}")
    shell("echo -e \"******\" >>{log.exec_time}")
