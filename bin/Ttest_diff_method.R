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

Ttest                     = snakemake@input$binary
no_GENCODE                = snakemake@input$counts
normalization_factor_path = snakemake@input$sample_conditions
pvalue_threshold          = snakemake@params$pvalue_threshold
log2fc_threshold          = snakemake@params$log2fc_threshold
conditionA                = snakemake@params$conditionA
conditionB                = snakemake@params$conditionB

output_diff_counts  = snakemake@output$diff_counts
output_pvalue_all   = snakemake@output$pvalue_all
output_tmp          = snakemake@output$tmp_dir
output_log          = snakemake@log[[1]]

# Function for logging to the output
logging <- function(str) {
  sink(file=paste(output_log), append=TRUE, split=TRUE)
  print(paste(Sys.time(),str))
  sink()
}

dir.create(output_tmp, showWarnings = FALSE)

logging(paste("Start Ttest_diff_methods"))
logging(paste("pvalue_threshold", pvalue_threshold))
logging(paste("log2fc_threshold", log2fc_threshold))

system(paste(Ttest ,"-p", pvalue_threshold, "-f", log2fc_threshold, "-r", output_pvalue_all, no_GENCODE, normalization_factor_path, conditionA, conditionB, "| gzip -c > ",output_diff_counts))

logging(paste("End Ttest_diff_methods"))
