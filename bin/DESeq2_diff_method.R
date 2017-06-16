library("data.table")
library("foreach")
library("doParallel")
library("DESeq2")

# Get parameters for the test
no_GENCODE                = snakemake@input$counts
sample_conditions         = snakemake@input$sample_conditions
pvalue_threshold          = snakemake@params$pvalue_threshold
log2fc_threshold          = snakemake@params$log2fc_threshold
conditionA                = snakemake@params$conditionA
conditionB                = snakemake@params$conditionB
nb_core                   = snakemake@threads
chunk_size                = snakemake@params$chunk_size

# Get output files
output_tmp          = snakemake@output$tmp_dir
output_diff_counts  = snakemake@output$diff_counts
output_pvalue_all   = snakemake@output$pvalue_all
output_log          = snakemake@log[[1]]


# Temporary files
shuffle_path              = paste(output_tmp,"/SHUFFLE_tmp.gz",sep="")
output_tmp_chunks         = paste(output_tmp,"/tmp_chunks/",sep="")
output_tmp_DESeq2         = paste(output_tmp,"/tmp_DESeq2/",sep="")
header_no_GENCODE         = paste(output_tmp,"/header_no_GENCODE.txt",sep="")
tmp_concat                = paste(output_tmp,"/tmp_concat.txt",sep="")
header_adj_pvalue         = paste(output_tmp,"/header_adj_pvalue.txt",sep="")
adj_pvalue                = paste(output_tmp,"/adj_pvalue.txt.gz",sep="")
sorted_adj_pvalue_tmp     = paste(output_tmp,"/sorted_adj_pvalue_tmp.txt",sep="")
dataDESeq2All             = paste(output_tmp,"/dataDESeq2All.txt.gz",sep="")
header_dataDESeq2All      = paste(output_tmp,"/header_dataDESeq2All.txt",sep="")
sorted_dataDESeq2All_tmp  = paste(output_tmp,"/sorted_dataDESeq2All_tmp.txt",sep="")
dataDESeq2Filtered        = paste(output_tmp,"/dataDESeq2Filtered.txt.gz",sep="")
final_header_tmp          = paste(output_tmp,"/final_header_tmp.txt.gz",sep="")

# Create directories
dir.create(output_tmp, showWarnings = FALSE)
dir.create(output_tmp_chunks, showWarnings = FALSE)
dir.create(output_tmp_DESeq2, showWarnings = FALSE)

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

logging("Start DESeq2_diff_methods")

# Check the chunk size
if(chunk_size > 1000000){
  logging(paste("Chunks too large for DESeq2 computations, reduce from",chunk_size,"to 1 000 000"))
  chunk_size = 1000000
}

# Set the number of cores to use
registerDoParallel(cores=nb_core)

# SHUFFLE THE KMERS (AVOIDING SHUFFLE THE HEADER)
system(paste("zcat", no_GENCODE, "| tail -n +2 | shuf | gzip -c >" ,shuffle_path))
logging("Shuffle done")

# SAVE THE HEADER INTO A FILE
system(paste("zcat", no_GENCODE, "| head -1 >", header_no_GENCODE))

# SPLIT THE MAIN FILE INTO CHUNKS WITH AUTOINCREMENTED NAMES
# FIXME make sure the directory is empty
# TODO these files should be gzipped to improve space efficency

system(paste("zcat", shuffle_path, "| awk -v", paste("chunk_size=", chunk_size,sep=""), "-v", paste("output_tmp_chunks=",output_tmp_chunks,sep=""),
             "'NR%chunk_size==1{OFS=\"\\t\";x=++i\"_subfile.txt.gz\"}{OFS=\"\";print | \"gzip >\" output_tmp_chunks x}'"))

# Remove shuffle path
system(paste("rm",shuffle_path))

nb_line_last_file = nbLineLastFile(output_tmp_chunks)
nb_files = nbFiles(output_tmp_chunks)

# CONCATENATE THE LAST FILE CREATED WITH THE SECOND LAST ONE AND THEN REMOVE IT
# IN ORDER TO AVOID TOO SHORT FILES
# IF THE MERGED FILE IS TOO LARGE, IT WILL BE DIVIDED IN TWO (jerome: it does not look like that we are doing it...)

if(nb_files > 1 && nb_line_last_file < (chunk_size/2)) {

  ## CONCATENATE THE 2 FILES
  logging(paste("The last file has",nb_line_last_file,"line(s) it will be concatenated to the second last one"))

  before_last_file = paste(output_tmp_chunks, (nb_files - 2), "_subfile.txt.gz", sep="")
  last_file        = paste(output_tmp_chunks, (nb_files - 1), "_subfile.txt.gz", sep="")

  system(paste("zcat", before_last_file, ">",  tmp_concat))
  system(paste("zcat", last_file,        ">>", tmp_concat))

  nb_line_last_file = chunk_size + nb_line_last_file

  logging(paste("The last file has",nb_line_last_file,"line(s) it will be splitted in two"))

  ### DIVIDE IN TWO PARTS
  system(paste("head -n", floor(nb_line_last_file/2),                          tmp_concat, "| gzip >", before_last_file))
  system(paste("tail -n", paste("+", floor(nb_line_last_file/2 + 1), sep=""),  tmp_concat, "| gzip >", last_file))
  system(paste("rm", tmp_concat))
}

## LOAD THE FILENAMES OF THE DIFFERENT CHUNKS
lst_files = system(paste("find",output_tmp_chunks,"-iname \"*_subfile.txt.gz\" | sort -n"), intern = TRUE)

logging("Split done")

## LOAD THE HEADER
header = as.character(unlist(read.table(file = header_no_GENCODE, sep = "\t", header = FALSE)))

logging(paste("Foreach of the", length(lst_files),"files"))

## LOADING PRIOR KNOWN NORMALISATION FACTORS
colData = read.table(sample_conditions,header=T,row.names=1)

## DESeq2 ANALYSIS ON EACH CHUNKS
invisible(foreach(i=1:length(lst_files)) %dopar% {
            bigTab = read.table(lst_files[i],header=F,stringsAsFactors=F)
            #SET TAGS AS ROWNAMES
            rownames(bigTab)=bigTab[,1]
            #REMOVE THE TAG AS A COLUMN
            bigTab=bigTab[,2:ncol(bigTab)]
            names(bigTab)=header[2:length(header)]

            countData = as.matrix(bigTab)

            dds <- DESeqDataSetFromMatrix(countData,
                                          colData = colData,
                                          design = ~ condition)

            NormCount_names = colnames(bigTab)
            rm(bigTab);gc()

            #RUN DESEQ2 AND COLLECT DESeq2 results

            #REPLACE SIZE FACTORS by SIZE FACTORS COMPUTED ON
            #THE ALL DATASET
            normFactors <- matrix(colData$normalization_factor,
                                  ncol=ncol(dds),nrow=nrow(dds),
                                  dimnames=list(1:nrow(dds),1:ncol(dds)),
                                  byrow = TRUE)

            normalizationFactors(dds) <- normFactors

            #RUN DESeq2
            dds <- estimateDispersionsGeneEst(dds)
            dds <- estimateDispersionsFit(dds)
            dds <- estimateDispersionsMAP(dds)
            dds <- nbinomWaldTest(dds)
            resDESeq2 <- results(dds, pAdjustMethod = "none")

            #COLLECT COUNTS
            NormCount<- as.data.frame(counts(dds, normalized=TRUE))
            names(NormCount) <- NormCount_names

            # WRITE A TSV WITH THIS FORMAT FOR THE CURRENT CHUNK
            # Kmer_ID
            # meanA
            # meanB
            # log2FC
            # NormCount
            write.table(data.frame(ID=rownames(resDESeq2),
                                   meanA=rowMeans(NormCount[,rownames(subset(colData, condition == conditionA))]),
                                   meanB=rowMeans(NormCount[,rownames(subset(colData, condition == conditionB))]),
                                   log2FC=resDESeq2$log2FoldChange,
                                   NormCount),
                        file=gzfile(paste(output_tmp_DESeq2,i,"_dataDESeq2_part_tmp.gz", sep="")),
                        sep="\t",quote=FALSE,
                        row.names = FALSE,
                        col.names = TRUE)

            # WRITE PVALUES FOR THE CURRENT CHUNK
            write.table(data.frame(ID=rownames(resDESeq2),pvalue=resDESeq2$pvalue),
                        file=gzfile(paste(output_tmp_DESeq2,i,"_pvalue_part_tmp.gz",sep="")),
                        sep="\t",quote=FALSE,
                        row.names = FALSE,
                        col.names = FALSE)

            # Remove processed chunk
            system(paste("rm",lst_files[i]))

}) #END FOREACH

system(paste("rm -rf", output_tmp_chunks))

logging("Foreach done")

#MERGE ALL CHUNKS PVALUE INTO A FILE
system(paste("find", output_tmp_DESeq2, "-name '*_pvalue_part_tmp.gz' | xargs cat >", output_pvalue_all))

logging(paste("Pvalues merged into",output_pvalue_all))

#MERGE ALL CHUNKS DESeq2 INTO A FILE
system(paste("find", output_tmp_DESeq2, "-name '*_dataDESeq2_part_tmp.gz' | xargs zcat | awk '{OFS=\"\\t\"}{if(NR==1){print}else{if($1 !~ \"ID\"){print}}}' | gzip >", dataDESeq2All))

logging(paste("DESeq2 results merged into",dataDESeq2All))

# Remove DESeq2 chunked results
system(paste("rm -rf", output_tmp_DESeq2))

# SAVE THE COLNAMES INTO HEADER.TXT AND REMOVE THE COLNAME ID
system(paste("zcat", dataDESeq2All, "| head -1 | cut -f2- >", header_dataDESeq2All))

#CREATE AND WRITE THE ADJUSTED PVALUE UNDER THRESHOLD WITH THEIR ID
pvalueAll         = read.table(output_pvalue_all, header=F, stringsAsFactors=F)
names(pvalueAll)  = c("ID","pvalue")
adjPvalue         = p.adjust(as.numeric(as.character(pvalueAll[,"pvalue"])),"BH")

adjPvalue_dataframe = data.frame(ID=pvalueAll$ID,
                                 pvalue=adjPvalue)

logging(paste("Removing k-mer with pvalue >",pvalue_threshold,"and abs(log2FC) <",log2fc_threshold))
logging(str(adjPvalue_dataframe))

adjPvalue_dataframe = subset(adjPvalue_dataframe, pvalue <= pvalue_threshold)

adjPvalue_dataframe = na.omit(adjPvalue_dataframe)

write.table(adjPvalue_dataframe,
            file=gzfile(adj_pvalue),
            sep="\t",
            quote=FALSE,
            col.names = TRUE,
            row.names = FALSE)

logging("Pvalues are adjusted")

#SAVE THE HEADER
system(paste("zcat", adj_pvalue, "|head -1 >", header_adj_pvalue))

#LEFT JOIN INTO dataDESeq2All
#GET ALL THE INFORMATION (ID,MEAN_A,MEAN_B,LOG2FC,COUNTS) FOR DE KMERS
system(paste("zcat", adj_pvalue, "| tail -n +2 | sort -k1,1 >", sorted_adj_pvalue_tmp))
system(paste("rm",adj_pvalue))
system(paste("zcat", dataDESeq2All,"| tail -n +2 | sort -k1,1 >", sorted_dataDESeq2All_tmp))
system(paste("rm",dataDESeq2All))
system(paste("join", sorted_adj_pvalue_tmp, sorted_dataDESeq2All_tmp, "| awk 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs($5) >=", log2fc_threshold, ") print $0}' | tr ' ' '\t' | gzip >", dataDESeq2Filtered))
system(paste("rm",sorted_adj_pvalue_tmp, sorted_dataDESeq2All_tmp))

logging("Get counts for pvalues that passed the filter")

#CREATE THE FINAL HEADER USING ADJ_PVALUE AND DATADESeq2ALL ONES AND COMPRESS THE FILE
system(paste("paste", header_adj_pvalue, header_dataDESeq2All, "| gzip >", final_header_tmp))
system(paste("cat", final_header_tmp, dataDESeq2Filtered, ">", output_diff_counts))
system(paste("rm", dataDESeq2Filtered))

logging("Analysis done")
