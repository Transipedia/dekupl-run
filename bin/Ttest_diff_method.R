Ttest=snakemake@input$Ttest_filter
no_GENCODE=snakemake@input$counts
normalization_factor_path = snakemake@input$sample_conditions
pvalue_threshold=snakemake@config$Ttest$pvalue_threshold
log2fc_threshold=snakemake@config$Ttest$log2fc_threshold

output_diff_counts=snakemake@output$diff_counts
output_pvalue_all=snakemake@output$pvalue_all

output_log=snakemake@log[[1]]

## GET THE PATH OF KMER_DE_DIR IN ORDER TO OUTPUT THE raw_pvals.txt FILE DIRECTLY IN IT
working_dir=strsplit(output_pvalue_all,"raw_pvals.txt")[[1]]

setwd(working_dir)

sink(file=paste(output_log), append=TRUE, split=TRUE)
print(paste(Sys.time(),"Start Ttest_diff_methods"))
sink()

system(paste(Ttest," -p ",pvalue_threshold, " -f ",log2fc_threshold," ",no_GENCODE," ",normalization_factor_path," A B | gzip -c > ",output_diff_counts,sep=""))

sink(file=paste(output_log), append=TRUE, split=TRUE)
print(paste(Sys.time(),"End Ttest_diff_methods"))
sink()
