library("data.table")

  ## SET UP VARS
data=snakemake@input$gene_counts
output_tmp=snakemake@config$tmp_dir
output_norm_factors=snakemake@output

output_log=snakemake@log[[1]]

dir.create(output_tmp, showWarnings = FALSE)

setwd(output_tmp)

sink(paste(output_log), append=TRUE, split=TRUE)
print(paste(Sys.time(),"Start normalization factors computation"))
sink()

  ## SAMPLING
  ## SELECT 33% OF THE TOTAL NUMBER OF K-MERS FOR THE SAMPLING
system(paste("zcat ",data," | awk '{if(NR % 3 ==0 || NR ==1){print $0}}' > selected_kmers",sep=""))

selected_kmers_counts <- data.frame(fread(paste("selected_kmers")))

sink(file=paste(output_log), append=TRUE, split=TRUE)
print(paste("Number of kmers :",nrow(selected_kmers_counts)-1))
sink()

loggeomeans <- rowMeans(log(selected_kmers_counts[,2:ncol(selected_kmers_counts)]))

normFactors <- apply(selected_kmers_counts[,2:ncol(selected_kmers_counts)], 2, function(cnts) { exp(median((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))})

 ## WRITING THE NORMALIZATION FACTORS
 ## SAMPLE NAMES \t NORMALIZATION FACTOR
write.table(data.frame(sample=names(normFactors),normalization_factor=as.vector(normFactors)),
             file = paste(output_norm_factors),
             sep="\t",
             quote=FALSE,
             row.names = FALSE)

sink(file=paste(output_log), append=TRUE, split=TRUE)
print(paste(Sys.time(),"Normalization factors computation done"))
sink()
