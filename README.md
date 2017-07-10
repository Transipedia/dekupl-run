# dekupl-run

DE-kupl is a pipeline that finds differentially expressed k-mers between RNA-Seq datasets.

Dekupl-run handles the first part of the [DE-kupl pipeline](https://github.com/Transipedia/dekupl) from raw FASTQ to
the production of contigs from differentially expressed k-mers.

## Dependencies

Before using Dekupl-run, install these dependencies:

- Snakemake
- jellyfish
- pigz
- CMake
- R: 
  * DESEq2 : open R and execute :
    `> source("https://bioconductor.org/biocLite.R")`
    `> biocLite("DESeq2")`
  * RColorBrewer
  * pheatmap
- Python: 
  * rpy2 : `pip3 install rpy2`

## Installation and usage

1. Clone this repository including submodules : `git clone --recursive https://github.com/Transipedia/dekupl-run.git`
2. Edit the config.json file to add the list of your samples, their conditions and the location their FASTQ file. See next section for parameters description.
3. Run the pipeline with then `snakemake -jNB_THREADS --resources ram=MAX_MEMORY -p` command. Replace `NB_THREADS` with the number of threads and `MAX_MEMORY` with the maximum memory (in Megabyte) you want DEkupl to allocate.
4. Once Dekupl-run has been fully executed, DE contigs produced by Dekupl-run
   (under `DEkupl_results/A_vs_B_kmer_counts/merged-diff-counts.tsv.gz`)
   can be annotate using [Dekupl-annotation](https://github.com/Transipedia/dekupl-annotation)

## Configuration (config.json)

- **fastq_dir**:  Location of FASTQ files
- **nb_threads**: Default number of thread to use (unless specified in the
  snakemake command-line
- **kmer_length**: Length of k-mers (default: 31). This value shoud not exceed
  32.
- **diff_method**: Method used for differential testing (default: DESeq2). Possible choices are 'Ttest' which is fast and 'DESeq2' which is more sensitive but longer to run.
- **lib_type**: Paired-end library type (default: `rf`). You can specify either `rf` for reverse-forward strand-specific libraries, `fr` for strand-specific forward-reverse, or `unstranded` for unstranded libraries.
- **output_dir**: Location of DE-kupl results (default: `DEkupl_result`).
- **tmp_dir**: Temporary directory to use (default: `./` aka current directory)
- **r1_suffix**: Suffix to use for the FASTQ with left mate. Set `r2_suffix` for the second FASTQ.
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

## Output files

The output directory of a DE-kupl will have the following content :

```
├── {A}_vs_{B}_kmer_counts
│   ├── diff-counts.tsv.gz
│   ├── merged-diff-counts.tsv.gz
├── gene_expression
│   ├── {A}vs{B}-DEGs.tsv
├── kmer_counts
│   ├── normalization_factors.tsv
│   ├── raw-counts.tsv.gz
│   ├── noGENCODE-counts.tsv.gz
│   ├── {sample}.jf
│   ├── {sample}.txt.gz
│   ├── ...
├── metadata
│   ├── sample_conditions.tsv
│   ├── sample_conditions_full.tsv
```

The following table describes the output files produced by DE-kupl :

FileName | Description
---------|------------
`diff-counts.tsv.gz` | Contains k-mers counts from `noGENCODE-counts.tsv.gz` that have passed the differential testing. Output format is a tsv with the following columns: `kmer pvalue meanA meanB log2FC [SAMPLES]`.
`merged-diff-counts.tsv.gz` | Contains assembled k-mers from `diff-counts.tsv.gz`. Output format is a tsv with the following columns: `nb_merged_kmers contig kmer pvalue meanA meanB log2FC [SAMPLES]`.
`raw-counts.tsv.gz` | Containins raw k-mer counts of all libraries that have been filtered with the reccurency filters.
`noGENCODE-counts.tsv.gz` | Containtains k-mer counts filtered from `raw-counts.tsv` with the k-mers from the reference transcription (ex: GENCODE by default).
`sample_conditions_full.tsv` | Tabulated file with samples names, conditions and normalization factors. `sample_conditions.tsv` is the sample

## FAQ
- if new samples are added to the config.json, make sure to remove the `metadata` folder in order to force SnakeMake to re-make all targets that depends on this file

## TODO

- Create a dekupl binary with two commands :
  - `dekupl build_index {genome}`:
    This command will download reference files and create all indexes
  - `dekupl run {dekupl_index} {config.yml} {output_dir}`:
    This command will run the dekupl pipeline
