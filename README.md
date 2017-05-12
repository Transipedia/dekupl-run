# dekupl-run

DE-kupl is a pipeline that finds differentially expressed k-mers between RNA-Seq datasets.

Dekupl-run handles the first part of the [DE-kupl pipeline]('https://github.com/Transipedia/dekupl') from raw FASTQ to
the production of contigs from differentially expressed k-mers.

## Dependencies

Before using Dekupl-run, install these dependencies:

- Snakemake
- jellyfish
- pigz
- [gsnap](http://research-pub.gene.com/gmap/)
- CMake
- R: 
  * DESEq2 : open R and execute :
    `> source("https://bioconductor.org/biocLite.R")`
    `> biocLite("DESeq2")`
- Python: 
  * rpy2 : `pip3 install rpy2`
- Perl: 
  * CracTools::Utils : `cpanm install CracTools::Utils`

## Installation and usage

1. Clone this repository including submodules : `git clone --recursive git@github.com:Transipedia/dekupl-run.git`
2. Edit the config.json file to add the list of your samples, their conditions and the location their FASTQ file. See next section for parameters description.
3. Run the pipeline with then `snakemake -jNB_THREADS -p` command. Replace `NB_THREADS` with the number of threads.

## Configuration (config.json)

- **fastq_dir**:  Location of FASTQ files
- **nb_threads**: Default number of thread to use (unless specified in the
  snakemake command-line
- **kmer_length**: Length of k-mers (default: 31). This value shoud not exceed
  32.
- **lib_type**: Paired-end library type (default: `rf`). You can specify either `rf` for reverse-forward strand-specific libraries, `fr` for strand-specific forward-reverse, or `unstranded` for unstranded libraries.
- **output_dir**: Location of DE-kupl results (default: `DEkupl_result`).
- **dekupl_counter**:
  * *min_reccurence*: Minimum number of samples to support a k-mer
  * *min_recurrence_abundance*: Min abundance threshold to consider a k-mer in
    the reccurency filter.
- **Ttest**:
  * *condition*: Specify A and B conditions.
  * *pvalue_threshold*: Min p-value (adjusted) to consider a k-mer as DE. Only
    DE k-mers are selected for assembly.
  * *log2fc_threshold*: Min Log2 Fold Change to consider a k-mer as DE.
- **Samples**: An array of samples. Each sample is described by a `name` and a
  `condition`. The FASTQ files for a sample will be located using the following
  command
    `fastq_dir/sample_name_{1,2}.fastq.gz`

## FAQ
- if new samples are added to the config.json, make sure to remove the `sample_conditions.tsv` file in order to force SnakeMake to re-make all targets that depends on this file

## TODO

- Create a dekupl binary with two commands :
  - `dekupl build_index {genome}`:
    This command will download reference files and create all indexes
  - `dekupl run {dekupl_index} {config.yml} {output_dir}`:
    This command will run the dekupl pipeline
