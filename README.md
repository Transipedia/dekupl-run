![dekupl-annot-logo](dekupl-run-logo.png)

# dekupl-run [![pipeline status](https://gitlab.com/transipedia/dekupl-run/badges/master/pipeline.svg)](https://gitlab.com/transipedia/dekupl-run/commits/master)


DE-kupl is a pipeline that finds differentially expressed k-mers between RNA-Seq datasets under The MIT License.

Dekupl-run handles the first part of the [DE-kupl pipeline](https://github.com/Transipedia/dekupl) from raw FASTQ to
the production of contigs from differentially expressed k-mers.

## Usage

Dekupl-run is a pipeline build with Snakemake. It works with a [configuration file](#configuration) that you will use to set the list of samples and their conditions as well as parameters for the test.

1. **Create a config.json** with the list of your samples, their conditions and the location their FASTQ file. See next section for parameters description.
2. **Run the pipeline**. Replace `CONFIG_JSON` with the config file you have created, `NB_THREADS` with the number of threads and `MAX_MEMORY` with the maximum memory (in Megabyte) you want DEkupl to allocate. This command line can varry depending of the installation (docker, singularity, manual, etc).
   `dekupl-run --configfile CONFIG_JSON -jNB_THREADS --resources ram=MAX_MEMORY -p`

3. **Explore results**. Once Dekupl-run has been successfully executed, DE contigs produced by Dekupl-run
   are located under `DEkupl_results/A_vs_B_kmer_counts/merged-diff-counts.tsv.gz`. They can be annoted
   can be annotate using [Dekupl-annotation](https://github.com/Transipedia/dekupl-annotation) and the vizualized with [Dekupl-viewer](https://github.com/Transipedia/dekupl-viewer).

## Installation

We recommand tu use [conda](https://anaconda.org/) to install dekupl-run, but you can also use Docker, Singularity and manual installatioN.

### Option1 : Use dekupl-run with conda

- **Step 1: Install conda.** If you do not have a conda distribution installed, we recommend to install miniconda as follow. See [Miniconda website](https://conda.io/miniconda.html) for other installation instructions (ex. for OSX).
    ```
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ``` 
- **Step 2: Install dekupl-run**. This will create a conda environment dekupl (if missing) and install dekupl run inside, the order of the parameters is important.
    ```
    conda install -n dekupl -y -m --override-channels -c transipedia \
     -c bioconda -c conda-forge -c https://repo.anaconda.com/pkgs/main \
     -c https://repo.anaconda.com/pkgs/free \
     -c https://repo.anaconda.com/pkgs/pro dekupl-run
    ```
- **Step 3: Run dekupl-run**. We first activate the conda environement where dekupl-run was installed, then we run the software.
    ```
    source activate dekupl
    dekupl-run --configfile my-config.json  -jNB_THREADS --resources ram=MAX_MEMORY -p
    ```

### Option 2: Use dekupl-run with Docker

- **Step 1: Retrieve the docker image.**
    ```
    docker pull transipedia/dekupl-run
    ```
- **Step 2: Run dekupl-run**.
    You may need to mount some volumes :
    - Your `my-config.json` to `/dekupl/my-config.json`
    - Your fastq_dir (the one defined in your `config.json`) to `/dekupl/FASTQ_DIR`
    - Your output_dir (the one defined in your `config.json`) to `/dekupl/OUTPUT_DIR`
    - Any other necessary folder depending on your `config.json`
    ```
    docker run --rm -v ${PWD}/my-config.json:/dekupl/my-config.json \
    -v ${PWD}/data:/dekupl/data  -v ${PWD}/results:/dekupl/results \
    transipedia/dekupl-run --configfile my-config.json  \
    -jNB_THREADS --resources ram=MAX_MEMORY -p
    ```

### Option 3: Use dekupl-run with singularity

We can create a singularity container from the docker image. Two method are available, they should both work. 

A difference with docker image is that with Singularity, you don't need to mount any volumes with singularity, but you must have your config.json and your inputs file in the directory where you are running dekupl-run.

- **Method 1**
  ```
    singularity pull docker://transipedia/dekupl-run
    ./dekupl-run.simg --configfile my-config.json -jNB_THREADS --resources ram=MAX_MEMORY -p
    ```
- **Method 2**
    ```
    singularity build dekupl-run.img docker://transipedia/dekupl-run
    singularity run ./dekupl-run.img --configfile my-config.json -jNB_THREADS --resources ram=MAX_MEMORY -p
    ```

### Option 4: Build and run yourself (not recommended)

- **Step 1: Install dependancies**. Before using Dekupl-run, install these dependencies:
    - Snakemake, jellyfish, pigz, CMake, boost
    - R packages (DESEq2, RColorBrewer, pheatmap, foreach, doParallel)
    `Rscript install_r_packages.R`
- **Step 2: . Clone this repository including submodules.**
  `git clone --recursive https://github.com/Transipedia/dekupl-run.git`
- **Step 3: Edit config file & run dekupl-run with Snakemake.**
  `snakemake -jNB_THREADS --resources ram=MAX_MEMORY -p`

## Configuration

Here is an example of a minimal config file with only mandatory informations. You can copy this base and adapt it to your needs (see following paragraph).

```
{
  "fastq_dir": "data",

  "dekupl_counter": {
    "min_recurrence": 2,
    "min_recurrence_abundance": 5
  },

  "diff_analysis": {
    "condition" : {
      "A": "A",
      "B": "B"
    },
    "pvalue_threshold": 0.05,
    "log2fc_threshold": 2
  },

  "samples": [{
      "name": "sample1",
      "condition": "A"
    }, {
      "name" : "sample2",
      "condition" : "A"
    }, {
      "name" : "sample3",
      "condition" : "B"
    }, {
      "name" : "sample4",
      "condition" : "B"
    }
  ]
}
```

### Parameters FAQ

**How can I use DEkupl-run with non-human data ?**
You need to specify your own FASTA using the `transcript_fasta` option as well as file with mapping of transcript_id to gene_id with the `transcript_to_gene` option.

**How can I use DEkupl-run with single-end reads?**
You need to specify the value *"single"* for the parameter `lib_type`. You can also specify frangments length (see section [Configuration for single-end libraries](#configuration-for-single-endlibraries))

### General configuration parameters

- **fastq_dir**:  Location of FASTQ files
- **nb_threads**: Default number of thread to use (unless specified in the
  snakemake command-line
- **kmer_length**: Length of k-mers (default: 31). This value shoud not exceed
  32.
- **diff_method**: Method used for k-mer differential testing (default: DESeq2). Possible choices are 'Ttest' which is fast and 'DESeq2' which is more sensitive but longer to run.
- **gene_diff_method**: Method used for gene differential testing (default: 'DESeq2' or 'limma-voom' if number of samples > 100). Possible choices are 'DESeq2' and 'limma-voom'. 'limma-voom' is a faster alternative for large cohorts.
- **lib_type**: Paired-end library type (default: `rf`). You can specify either `rf` for reverse-forward strand-specific libraries, `fr` for strand-specific forward-reverse, or `unstranded` for unstranded libraries.
- **output_dir**: Location of DE-kupl results (default: `DEkupl_result`).
- **tmp_dir**: Temporary directory to use (default: `./` aka current directory)
- **r1_suffix**: Suffix to use for the FASTQ with left mate. Set `r2_suffix` for the second FASTQ.
- **dekupl_counter**:
  * *min_recurrence*: Minimum number of samples to support a k-mer
  * *min_recurrence_abundance*: Min abundance threshold to consider a k-mer in
    the reccurency filter.
- **diff_analysis**:
  * *condition*: Specify A and B conditions.
  * *pvalue_threshold*: Min p-value (adjusted) to consider a k-mer as DE. Only
    DE k-mers are selected for assembly.
  * *log2fc_threshold*: Min Log2 Fold Change to consider a k-mer as DE.
- **Samples**: An array of samples. Each sample is described by a `name` and a
  `condition`. The FASTQ files for a sample will be located using the following
  command
    `fastq_dir/sample_name_{1,2}.fastq.gz`
- **transcript_fasta**: The reference transcriptome to be used for masking. By default DEKupl-run uses the human Gencode transcriptome for masking. To change this, add to the config.json file:
`"transcript_fasta":my_transciptome.fa`
- **transcript_to_gene**: This is a two column tabulated file, with the transcript ID in the first column and the gene ID in the second column. The file is not mandatory if the FASTA transcriptome is from Gencode, were the gene ID can be extracted from the sequence names in the FASTA. An example of this file can be found here : [tests/gencode.v24.transcripts.head1000.mapping.tsv](tests/gencode.v24.transcripts.head1000.mapping.tsv).

### Configuration for single-end libraries

For single-end libraries please specify the following parameters :

- **lib_type**: You can either set the lib_type to `single` in the case of single-end strand-specific library or `unstranded` for single-end unstranded libraries.
- **fragment_length** : The estimated fragment length (necessary for kallisto quantification). Default value is `200`.
- **fragment_standard_deviation** : The estimated standard deviation of fragment length (necessary for kallisto quantification). Default value is `30`.

*Notes* :
The fastq files for single-end samples will be located using the following path : `{fastq_dir}/{sample_name}.fastq.gz`
If present, parameters **r1_suffix** and **r2_suffix** will be ignored.

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

## Whole-genome data

It is now possible to run DE-kupl-style analysis on whole-genome data, i.e. without using a reference transcriptome.
To do so, please change `data_type` to `WGS` in `config.json`.

## FAQ
- if new samples are added to the config.json, make sure to remove the `metadata` folder in order to force SnakeMake to re-make all targets that depends on this file
- Snakemake uses Rscript, not R. If a R module is not installed, type `which Rscript` and `which R` and make sure they point to the same installation of R.
- For OSX support you need to install the coreutils package with HomeBrew `brew install coreutils`. This package provide Linux versions of famus Unix command like "sort", "join", etc.