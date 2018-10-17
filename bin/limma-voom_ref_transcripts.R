#######################################################################
# The MIT License
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

library(limma)
library(edgeR)
# library(RColorBrewer)
# library(pheatmap)
# library(ggplot2)

args <- commandArgs(TRUE)

# Get parameters for the test
gene_counts       = args[1]#snakemake@input$gene_counts
sample_conditions = args[2]#snakemake@input$sample_conditions
condition_col     = args[3]#snakemake@params$condition_col
condition_A       = args[4]#snakemake@params$condition_A
condition_B       = args[5]#snakemake@params$condition_B

# Get output files
differentially_expressed_genes  = args[6]#snakemake@output$differentially_expressed_genes
norm_counts		                  = args[7]#snakemake@output$norm_counts
output_log                      = args[8]#snakemake@log[[1]]
#dist_matrix			              = args[10]#snakemake@output$dist_matrix
#pca_design			                = args[11]#snakemake@output$pca_design

write(date(),file=output_log)

# Debug files
# gene_counts <- "DEkupl_result_transcript_to_gene_mapping/gene_expression/kallisto/gene_counts.tsv.gz"
# sample_conditions <- "DEkupl_result_transcript_to_gene_mapping/metadata/sample_conditions_full.tsv"

# Load counts data
countsData = read.table(gene_counts,
                        header=T,
                        row.names=1,
                        check.names=FALSE)

# Load col data with sample specifications
colData = read.table(sample_conditions,
                     header=T,
                     row.names=1,
                     check.names=FALSE)

## remove genes with 0 count
nullGenes   <- rownames(countsData[rowSums(countsData)==0,])
countsData  <- countsData[rowSums(countsData)!=0,]

## perform limma-voom
dge <- DGEList(countsData, group=colData$condition)
v   <- voom(dge, plot=TRUE)
fit <- lmFit(v)
fit <- eBayes(fit, robust=FALSE)

# writing in a file normalized counts
# FIXME add back genes with 0 counts
normalized_counts <- data.frame(id=row.names(v$E),v$E,row.names=NULL)
write.table(normalized_counts, file=norm_counts, sep="\t", row.names=F, col.names=T, quote=F)

# Write DEGs to differentially_expressed_genes file
results <- topTable(fit, number = Inf)
results <- results[,c("AveExpr", "logFC", "P.Value","adj.P.Val")]
colnames(results) <- c("baseMean", "log2FoldChange", "pvalue", "padj")
# Add missing columns and re-order them to match DESeq2 output
# #baseMean        log2FoldChange  lfcSE   stat    pvalue  padj
results$lfcSE <- NA
results$stat <- NA
results <- results[,c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
write.table(results,file=differentially_expressed_genes,sep="\t",quote=FALSE)

# sampleDists<-dist(t(v$E) )
# sampleDistMatrix<-as.matrix( sampleDists )
# rownames(sampleDistMatrix)<-colnames(v$E)
# colnames(sampleDistMatrix)<-colnames(v$E)
# colours=colorRampPalette(rev(brewer.pal(9,"Blues")) )(255)
# 
# pdf(dist_matrix,width=15,height=10)
# 
# pheatmap(sampleDistMatrix,
#          main="clustering of samples",
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colours,
#          fontsize = 14)
# 
# data <- plotPCA(v$E,ntop=nrow(v$E),returnData=TRUE)
# write.table(data,pca_design,row.names=F, col.names=T, quote=F,sep="\t")
# 
# print(ggplot(data,aes(PC1,PC2,color=condition))+geom_point()+geom_text(aes(label=name),hjust=0,vjust=0))
# 
# dev.off()
