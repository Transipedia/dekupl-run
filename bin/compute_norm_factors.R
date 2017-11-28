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

selected_kmers_counts <- read.table(selected_kmers, header=T, stringsAsFactor=F)

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
