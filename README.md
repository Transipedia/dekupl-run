# DE-kupl

DE-kupl is a pipline, that finds differentially expressed k-mers between RNA-Seq dataset.

# INSTALL

Before using Dekupl, install these dependencies.

- Snakemake
- jellyfish
- pigz
- gsnap
- CMake
- R: 
  * DESEq2 : open R and execute :
    `> source("https://bioconductor.org/biocLite.R")`
    `> biocLite("DESeq2")`
- Python: 
  * rpy2 : `pip3 install rpy2`
- Perl: 
  * CracTools::Utils : `cpanm install CracTools::Utils`

# USAGE

1. Clone this repository to your cluster
2. Edit the config.json file to add the list of your samples, their conditions and the location their FASTQ files
3. Run the pipeline with then `snakemake -jNB_THREADS -p` command 

# FAQ

- if new samples are added to the config.json, make sure to remove the sample_conditions.tsv file in order to force SnakeMake to re-make all targets that depends on this file

# TODO

- Create a dekupl binary with two commands :
  - `dekupl build_index {genome}`:
    This command will download reference files and create all indexes
  - `dekupl run {dekupl_index} {config.yml} {output_dir}`:
    This command will run the dekupl pipeline
