repository <- "https://cloud.r-project.org"
install.packages("RColorBrewer", repos=repository)
install.packages("pheatmap", repos=repository)
install.packages("foreach", repos=repository)
install.packages("doParallel", repos=repository)
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
