library("data.table")

## SET UP VARS
data                = snakemake@input$raw_counts
output_tmp          = snakemake@output$tmp_dir
output_norm_factors = snakemake@output$nf
output_log          = snakemake@log[[1]]

selected_kmers      = paste(output_tmp,"/selected_kmers.txt.gz",sep="")

# Function for logging to the output
logging <- function(str) {
  sink(file=paste(output_log), append=TRUE, split=TRUE)
  print(paste(Sys.time(),str))
  sink()
}

dir.create(output_tmp, showWarnings = FALSE)
#setwd(output_tmp)

logging("Start normalization factors computation")

## SAMPLING
## SELECT 33% OF THE TOTAL NUMBER OF K-MERS FOR THE SAMPLING
system(paste("zcat", data ," | awk '{if(NR % 3 ==0 || NR ==1){print $0}}' |gzip -c >", selected_kmers))

selected_kmers_counts <- read.table(selected_kmers, header=T, stringsAsFactor=F, check.names=F)

logging(paste("Number of kmers :",nrow(selected_kmers_counts)-1))

loggeomeans <- rowMeans(log(selected_kmers_counts[,2:ncol(selected_kmers_counts)]))

normFactors <- apply(selected_kmers_counts[,2:ncol(selected_kmers_counts)], 2, function(cnts) { exp(median((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))})

## WRITING THE NORMALIZATION FACTORS
## SAMPLE NAMES \t NORMALIZATION FACTOR
write.table(data.frame(sample=names(normFactors),normalization_factor=as.vector(normFactors)),
             file       = output_norm_factors,
             sep        = "\t",
             quote      = FALSE,
             row.names  = FALSE)

logging("Normalization factors computation done")
