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

library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)

# Get parameters for the test
gene_counts       = snakemake@input$gene_counts
sample_conditions = snakemake@input$sample_conditions
condition_col     = snakemake@params$condition_col
condition_A       = snakemake@params$condition_A
condition_B       = snakemake@params$condition_B

# Get output files
differentially_expressed_genes  = snakemake@output$differentially_expressed_genes
dist_matrix			                = snakemake@output$dist_matrix
norm_counts		                  = snakemake@output$norm_counts
pca_design			                = snakemake@output$pca_design
output_log                      = snakemake@log[[1]]

write(date(),file=output_log)

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

write(colnames(countsData),stderr())
write(rownames(colData),stderr())

colData = colData[colnames(countsData),,drop=FALSE]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=countsData,
                              colData=colData,
                              design = ~ condition)
dds <- DESeq(dds)

#normalized counts
NormCount<- as.data.frame(counts(dds, normalized=TRUE ))

#writing in a file normalized counts
normalized_counts<-data.frame(id=row.names(NormCount),NormCount,row.names=NULL)
write.table(normalized_counts,file=norm_counts, sep="\t",row.names=F, col.names=T, quote=F)

write(resultsNames(dds),stderr())

# Write DEGs
res <- results(dds, contrast = c(condition_col,condition_A,condition_B))
write.table(res,file=differentially_expressed_genes,sep="\t",quote=FALSE)

rld<-rlog(dds)
sampleDists<-dist(t(assay(rld) ) )
sampleDistMatrix<-as.matrix( sampleDists )
rownames(sampleDistMatrix)<-colnames(rld)
colnames(sampleDistMatrix)<-colnames(rld)
colours=colorRampPalette(rev(brewer.pal(9,"Blues")) )(255)

pdf(dist_matrix,width=15,height=10)

pheatmap(sampleDistMatrix,
main="clustering of samples",
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colours,
fontsize = 14)

data <- plotPCA(rld,ntop=nrow(rld),returnData=TRUE)
write.table(data,pca_design,row.names=F, col.names=T, quote=F,sep="\t")

print(ggplot(data,aes(PC1,PC2,color=condition))+geom_point()+geom_text(aes(label=name),hjust=0,vjust=0))

dev.off()

write(date(),file=output_log,append=T)
