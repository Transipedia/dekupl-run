#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

if (length(args)==0) {
  stop("Usage; TtestFilter.R counts.tsv conditions.tsv > counts-filtered.tsv")
}

# Retrieve arguments
counts_file <- args[1]
conditions_file <- args[2]
condition_col <- args[3]
conditionA <- args[4]
conditionB <- args[5]
pvalue_threshold <- args[6]
log2fc_threshold <- args[7]

write(paste("condition A: ",conditionA), stderr())
write(paste("condition A: ",conditionB), stderr())

# Parse the sample conditions files
sample_conditions <- read.table(conditions_file,header=T,row.names=1)

#sample_conditions$condition = as.factor(sample_conditions$condition)
sample_conditions[[condition_col]] = as.factor(sample_conditions[[condition_col]])

# Parse the counts file
con  <- file(counts_file, open = "r")

header_line <- readLines(con, n = 1)

header <- strsplit(header_line,"\\s+")
header <- header[[1]]

nb_samples <- length(header) - 1
samples <- header[2:(nb_samples +1)]

normalization_factor <- as.numeric(sample_conditions[header[2:(nb_samples+1)],]$normalization_factor)

conditionA_indicies <- grep(paste(rownames(sample_conditions[ sample_conditions[[condition_col]] == conditionA, ]),
                                  collapse="|"),
                            header[2:(nb_samples+1)])

conditionB_indicies <- grep(paste(rownames(sample_conditions[ sample_conditions[[condition_col]] == conditionB, ]),
                                  collapse="|"),
                            header[2:(nb_samples+1)])

write(paste("Group A: ", paste(samples[conditionA_indicies],collapse=", ")), stderr())
write(paste("Group B: ", paste(samples[conditionB_indicies],collapse=", ")), stderr())

cat("tag\tpvalue\tmeanA\tmeanB\tlog2FC\t",paste(samples,collapse="\t"),"\n")

while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  values <- strsplit(line,"\\s+")
  values <- values[[1]]

  tag    <- values[1]
  values <- as.numeric(values[2:(nb_samples + 1)])
  values <- values / normalization_factor

  test <- t.test(values[conditionA_indicies],values[conditionB_indicies])
  #test <- wilcox.test(values[conditionA_indicies],values[conditionB_indicies])
  meanA = test$estimate[1]
  meanB = test$estimate[2]
  log2FC = log2(meanA/meanB)
  if(!is.nan(test$p.value) && test$p.value < pvalue_threshold && abs(log2FC) >= log2fc_threshold) {
    #write(test$p.value,stderr())
    #write(paste(values[conditionA_indicies],collapse=", "),stderr())
    #write(paste(values[conditionB_indicies],collapse=", "),stderr())
    cat(tag,"\t",test$p.value,"\t",meanA,"\t",meanB,"\t",log2FC,"\t",paste(round(values,2),collapse="\t"),"\n")
    #Sys.sleep(1)
  }
}

