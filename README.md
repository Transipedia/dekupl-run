# Dependencies

- Snakemake
- CMake
- R: DESEq2 :   
- Python: rpy2 pip3 install rpy2 --user
- Perl: CracTools::Utils : cpanm install CracTools::Utils
- Kallisto (downloaded automatically)

# FAQ

- if new samples are added to the config.json, make sure to remove the sample_conditions.tsv file in order to force SnakeMake to re-make all targets that depends on this file
