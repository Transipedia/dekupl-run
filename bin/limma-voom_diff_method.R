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
library("foreach")
library("doParallel")
library("limma")
library("edgeR")
library("DESeq2")

args <- commandArgs(TRUE)

# Get parameters for the test
binary                    = args[1]#snakemake@input$binary
kmer_counts               = args[2]#snakemake@input$counts
sample_conditions         = args[3]#snakemake@input$sample_conditions
pvalue_threshold          = args[4]#snakemake@params$pvalue_threshold
log2fc_threshold          = args[5]#snakemake@params$log2fc_threshold
nb_core                   = args[6]#snakemake@threads
chunk_size                = as.numeric(args[7])#snakemake@params$chunk_size

# Get output files
output_tmp                = args[8]#snakemake@output$tmp_dir
pre_output_diff_counts    = args[9]#snakemake@output$diff_counts
pre_output_pvalue_all     = args[10]#snakemake@output$pvalue_all
output_log                = args[11]#snakemake@log[[1]]
seed                      = args[12]#snakemake@params$seed

# Get conditions and contrast
nb_condition              = as.numeric(args[13])#snakemake@nb_condition
conditions                = args[14:(14+nb_condition-1)]#snakemake@conditions
contrast                  = args[(14+nb_condition):length(args)]#snakemake@contrast

# Temporary files
output_tmp_chunks         = paste(output_tmp,"/tmp_chunks/",sep="")
output_tmp_LimmaVoom      = paste(output_tmp,"/tmp_LimmaVoom/",sep="")
header_kmer_counts        = paste(output_tmp,"/header_kmer_counts.txt",sep="")
tmp_concat                = paste(output_tmp,"/tmp_concat.txt",sep="")
pre_adj_pvalue            = paste(output_tmp,"/adj_pvalue",sep="")
pre_dataLimmaVoomAll      = paste(output_tmp,"/dataLimmaVoomAll",sep="")
pre_dataLimmaVoomFiltered = paste(output_tmp,"/dataLimmaVoomFiltered",sep="")

# Create directories
dir.create(output_tmp, showWarnings = FALSE, recursive = TRUE)
dir.create(output_tmp_chunks, showWarnings = FALSE, recursive = TRUE)
dir.create(output_tmp_LimmaVoom, showWarnings = FALSE, recursive = TRUE)

# Function for logging to the output
logging <- function(str) {
  sink(file=paste(output_log), append=TRUE, split=TRUE)
  print(paste(Sys.time(),str))
  sink()
}
# Return the number of line in the last files of the directory
nbLineLastFile <- function(dir) {
  return(as.numeric(system(paste("zcat", paste(dir, "/$(ls ", dir, " |sort -n|grep subfile|tail -1)",sep=""), "| wc -l"), intern=TRUE)))
}
# Return the number of files in the directory
nbFiles <- function(dir) {
  return(as.numeric(system(paste("ls ", dir, "|grep subfile | wc -l"), intern=TRUE)))
}

logging("Start LimmaVoom_diff_methods")

# Check the chunk size
if(chunk_size > 1000000){
  logging(paste("Chunks too large for LimmaVoom computations, reduce from",chunk_size,"to 1 000 000"))
  chunk_size = 1000000
}

# Set the number of cores to use
registerDoParallel(cores=nb_core)

# CLEAN THE TMP FOLDER
system(paste("rm -f ", output_tmp_chunks, "/*", sep=""))

# SAVE THE HEADER INTO A FILE
system(paste("zcat", kmer_counts, "| head -1 | cut -f2- >", header_kmer_counts))

# SHUFFLE AND SPLIT THE MAIN FILE INTO CHUNKS WITH AUTOINCREMENTED NAMES
if(seed == 'fixed'){
    system(paste("zcat", kmer_counts, " >tmp_shuff; cat tmp_shuff| tail -n +2 | shuf --random-source=tmp_shuff | awk -v", paste("chunk_size=", chunk_size,sep=""), "-v", paste("output_tmp_chunks=",output_tmp_chunks,sep=""),
             "'NR%chunk_size==1{OFS=\"\\t\";x=++i\"_subfile.txt.gz\"}{OFS=\"\";print | \"gzip >\" output_tmp_chunks x}'"))
    system("rm tmp_shuff")
}else{
    system(paste("zcat", kmer_counts, "| tail -n +2 | shuf | awk -v", paste("chunk_size=", chunk_size,sep=""), "-v", paste("output_tmp_chunks=",output_tmp_chunks,sep=""),
                 "'NR%chunk_size==1{OFS=\"\\t\";x=++i\"_subfile.txt.gz\"}{OFS=\"\";print | \"gzip >\" output_tmp_chunks x}'"))
}
logging("Shuffle and split done")

nb_line_last_file = nbLineLastFile(output_tmp_chunks)
nb_files = nbFiles(output_tmp_chunks)

# IF THE LAST FILE HAS LESS THAN HALF OF THE CHUNK SIZE
# CONCATENATE THE LAST 2 FILES AND THEN SPLIT
if(nb_files > 1 && nb_line_last_file < (chunk_size/2)) {

  ## CONCATENATE THE 2 FILES
  logging(paste("The last file has",nb_line_last_file,"line(s) it will be concatenated to the second last one"))

  before_last_file = paste(output_tmp_chunks, (nb_files - 1), "_subfile.txt.gz", sep="")
  last_file        = paste(output_tmp_chunks, (nb_files), "_subfile.txt.gz", sep="")

  #CONCATENATE THE LAST 2 FILES INTO A TMP FILE
  system(paste("cat", before_last_file, last_file, ">", tmp_concat, sep=" "))

  #NUMBER OF LINE OF THE TMP FILE
  nb_line_last_file = chunk_size + nb_line_last_file

  logging(paste("The last file has", nb_line_last_file, "line(s) it will be splitted in two"))

  ### DIVIDE IN TWO PARTS
  system(paste("zcat", tmp_concat, "| head -n", floor(as.integer(nb_line_last_file/2)), "| gzip >", before_last_file))
  system(paste("zcat", tmp_concat, "| tail -n", paste("+", floor(as.integer(nb_line_last_file/2 + 1)), sep=""), "| gzip >", last_file))
  system(paste("rm", tmp_concat))
}

## LOAD THE FILENAMES OF THE DIFFERENT CHUNKS
lst_files = system(paste("find",output_tmp_chunks,"-iname \"*_subfile.txt.gz\" | sort -n"), intern = TRUE)

logging("Split done")

## LOAD THE HEADER
header = as.character(unlist(read.table(file = header_kmer_counts, sep = "\t", header = FALSE)))

## LOADING PRIOR KNOWN NORMALISATION FACTORS
colData = read.table(sample_conditions,header=T,row.names=1)

##PREPATATION OF LIMMA-VOOM
#Group
group=colData$condition
#Design
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
#Contrast
if (contrast[1] == 'NA'){
   x.contrast=paste(conditions[2],conditions[1],sep="-")
}else{
   x.contrast=contrast
}
print(x.contrast)
print(group)
print(colnames(design))
contr.matrix <-makeContrasts(contrasts=x.contrast,levels=colnames(design))
print(contr.matrix)

logging(paste("Foreach of the", length(lst_files),"files"))

## LimmaVoom ANALYSIS ON EACH CHUNKS
invisible(foreach(i=1:length(lst_files)) %dopar% {
            ##READ AND FORMAT DATA
            countData = read.table(lst_files[i],header=F,stringsAsFactors=F)
            #SET TAGS AS ROWNAMES
            rownames(countData)=countData[,1]
            #REMOVE THE TAG AS A COLUMN
            countData=countData[,2:ncol(countData)]
            names(countData)=header
            #FORMATING DATA
            dge <- DGEList(count=countData,group=group)

            #REPLACE SIZE FACTORS by SIZE FACTORS COMPUTED ON THE ALL DATASET
            normFactors <- c(t(matrix(colData$normalization_factor)))
            dge$samples$norm.factors <- normFactors
            
            ##RUN LIMMA-VOOM AND COLLECT Limma-voom results
            
            #RUN Limma-voom
            v <- voom(dge,design)
            fitLimmaVoom <-lmFit(v)
            fitLimmaVoom <- contrasts.fit(fitLimmaVoom, contrasts=contr.matrix)
            fitLimmaVoom <-eBayes(fitLimmaVoom, robust=FALSE)
            
            for (j in 1:length(x.contrast)) {
               resLimmaVoom <-topTable(fitLimmaVoom, number = nrow(countData), coef=j,sort.by=NULL,resort.by=NULL)
               #COLLECT COUNTS  (in cpm)
               NormCount<- as.data.frame(cpm(dge, log=F,normalized.lib.sizes=TRUE))
               names(NormCount) <- colnames(countData)
               NormCount=cbind(NormCount,ID=rownames(NormCount))
               
               # WRITE A TSV WITH THIS FORMAT FOR THE CURRENT CHUNK
               # Kmer_ID, mean,..., mean, log2FC, NormCount
               #recuperation des moyennes selon les conditions presentes dans le contraste
               #conditions presentes dans le contraste
               print(x.contrast[j])
               cond_contrast=intersect(unique(strsplit(x.contrast[j], "-|\\+|/|\\(|\\)")[[1]]), conditions)
               #calculs des moyennes (en cpm)
               means=matrix(NA,length(rownames(resLimmaVoom)),length(cond_contrast))
               for (k in 1:length(cond_contrast)){
                   means[,k]=rowMeans(NormCount[,rownames(subset(colData,condition == cond_contrast[k]))])
               }
               means=cbind(as.data.frame(means, row.names=rownames(NormCount)),rownames(NormCount))
               names(means)=c(cond_contrast,"ID")
               #ecriture
               dataframe_tmp=data.frame(ID=rownames(resLimmaVoom),log2FC=resLimmaVoom$logFC)
               dataframe_tmp=merge(means, dataframe_tmp, by="ID")

               write.table(merge(dataframe_tmp, NormCount, by="ID"),
                           file=gzfile(paste(output_tmp_LimmaVoom,i,"_dataLimmaVoom_part_tmp",j,".txt.gz", sep="")),
                           sep="\t",quote=FALSE,
                           row.names = FALSE,
                           col.names = FALSE)

               # WRITE PVALUES FOR THE CURRENT CHUNK
               write.table(data.frame(ID=rownames(resLimmaVoom),pvalue=resLimmaVoom$P.Value),
                           file=gzfile(paste(output_tmp_LimmaVoom,i,"_pvalue_part_tmp",j,".txt.gz",sep="")),
                           sep="\t",quote=FALSE,
                           row.names = FALSE,
                           col.names = FALSE)
            }
            # Remove processed chunk
            system(paste("rm",lst_files[i]))

}) #END FOREACH

system(paste("rm -rf", output_tmp_chunks))

logging("Foreach done")
print(x.contrast)
for (j in 1:length(x.contrast)){
	print(x.contrast[j])
	output_pvalue_all     = paste(pre_output_pvalue_all,j,".txt.gz",sep="")
	dataLimmaVoomAll      = paste(pre_dataLimmaVoomAll,j,".txt.gz",sep="")
	adj_pvalue            = paste(pre_adj_pvalue,j,".txt.gz",sep="")
	dataLimmaVoomFiltered = paste(pre_dataLimmaVoomFiltered,j,".txt.gz",sep="")
	output_diff_counts    = paste(pre_output_diff_counts,j,".tsv.gz",sep="")

	#MERGE ALL CHUNKS PVALUE INTO A FILE
	pvalueToFind=paste("'*_pvalue_part_tmp",j,".txt.gz'",sep="")
	system(paste("find", output_tmp_LimmaVoom, "-name",pvalueToFind,"| xargs cat >",output_pvalue_all))
	logging(paste("Pvalues merged into",output_pvalue_all))

	#MERGE ALL CHUNKS LimmaVoom INTO A FILE
	resultsToFind=paste("'*_dataLimmaVoom_part_tmp",j,".txt.gz'",sep="")
	system(paste("find", output_tmp_LimmaVoom, "-name",resultsToFind,"| xargs cat >", dataLimmaVoomAll))
	logging(paste("LimmaVoom results merged into",dataLimmaVoomAll))

	#CREATE AND WRITE THE ADJUSTED PVALUE UNDER THRESHOLD WITH THEIR ID
	pvalueAll         = read.table(output_pvalue_all, header=F, stringsAsFactors=F)
	names(pvalueAll)  = c("tag","pvalue")
	adjPvalue         = p.adjust(as.numeric(as.character(pvalueAll[,"pvalue"])),"BH")

	adjPvalue_dataframe = data.frame(tag=pvalueAll$tag,
									 pvalue=adjPvalue)

	write.table(adjPvalue_dataframe,
				file=gzfile(adj_pvalue),
				sep="\t",
				quote=FALSE,
				col.names = FALSE,
				row.names = FALSE)

	logging(paste("Pvalues are adjusted for",x.contrast[j]))

	#LEFT JOIN INTO dataLimmaVoomAll
	#GET ALL THE INFORMATION (ID,MEAN,...,MEAN,LOG2FC,COUNTS) FOR DE KMERS
	cond_contrast=intersect(unique(strsplit(x.contrast[j], "-|\\+|/|\\(|\\)")[[1]]), conditions)
	num_column_logFC=paste("(abs(\\$", 3+length(cond_contrast), ")", sep="") #3=tag+pvalue+log2FC
	system(paste("echo \"join <(zcat ", adj_pvalue,") <(zcat ", dataLimmaVoomAll," ) | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {if", num_column_logFC, ">=", log2fc_threshold, " && \\$2 <= ", pvalue_threshold, ") print \\$0}' | tr ' ' '\t' | gzip > ", dataLimmaVoomFiltered, "\" | bash", sep=""))
	system(paste("rm", adj_pvalue, dataLimmaVoomAll))

	logging("Get counts for pvalues that passed the filter")

	#CREATE THE FINAL HEADER USING ADJ_PVALUE AND DATALimmaVoomALL ONES AND COMPRESS THE FILE
	head=paste("'tag\tpvalue",paste(paste0("mean_",cond_contrast),collapse="\t"),"log2FC'",sep="\t")
	system(paste("echo",head,"| paste - ", header_kmer_counts," | gzip > ", output_diff_counts))
	system(paste("cat", dataLimmaVoomFiltered, ">>", output_diff_counts))
	system(paste("rm", dataLimmaVoomFiltered))
}

# REMOVE LimmaVoom CHUNKS RESULTS
system(paste("rm -rf", output_tmp_LimmaVoom))

logging("Analysis done")
