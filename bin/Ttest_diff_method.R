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

system(paste(Ttest ,"-p", pvalue_threshold, "-f", log2fc_threshold, "-r", output_pvalue_all, no_GENCODE, normalization_factor_path, conditionA, conditionB, "| gzip -c > ",output_diff_counts))

logging(paste("End Ttest_diff_methods"))
