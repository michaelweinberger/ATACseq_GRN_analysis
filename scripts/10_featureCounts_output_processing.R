### Run this script after mapping RNA-seq data to process featureCounts output




### user defined variables

# unpack variables passed from parent shell script
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)
args

for (e in args) {
  argname <- e[1]
  argval <- e[2]
  # regular expression to delete initial \" and trailing \"
  argval <- gsub("(^\\\"|\\\"$)", "", argval)
  assign(argname, argval)
}





### packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("rtracklayer")) BiocManager::install("rtracklayer")




### functions

## function to read in and process featureCounts file
# feature_counts_file_raw: Raw featureCounts output file
# gtf_file: Genome .gtf file
# out_dir: Directory containing processed featureCounts file
fc_process <- function(feature_counts_file_raw, gtf_file, out_dir) {
  
  # read in featureCounts output file 
  data <- read.table(feature_counts_file_raw, skip=1, sep="\t", header=TRUE)
  data <- data[,c(1,6:ncol(data))]
  data$Geneid <- gsub("\\.[0-9][0-9]", "", data$Geneid)
  data$Geneid <- gsub("\\.[0-9]", "", data$Geneid)
  
  print("Dimensions of raw featureCounts file: ")
  print(dim(data))
  print("This is what the file looks like: ")
  print(data, max=40)
  
  # create a mapping of ensembl gene id to gene name
  gtf_df <- as.data.frame(rtracklayer::import(gtf_file))
  genes <- unique(gtf_df[,c("gene_id","gene_name")])
  genes <- genes[!is.na(genes$gene_name),]

  # add the gene names
  data_1 <- merge(genes, data, by.x="gene_id", by.y="Geneid")
  
  # extract sample names from original column names 
  # (the full file paths of bam files analysed via featureCounts)
  # also cut off the ".rmdup.bam" file name ending
  # featureCounts replaced "/" with "."
  partitions <- lengths(regmatches(colnames(data)[3:ncol(data)], 
                                   gregexpr("\\.", colnames(data)[3:ncol(data)])))[1]
  sample_names <- matrix(unlist(strsplit(colnames(data)[3:ncol(data)], "\\.")), 
                         ncol = (partitions + 1), byrow = TRUE)[,(partitions - 1)]
  
  # set column names
  colnames(data_1) <- c("gene_id", "gene_name", "gene_length", sample_names)
  
  print("Dimensions of featureCounts file with gene names: ")
  print(dim(data_1))
  print("This is what the file looks like: ")
  print(data_1, max=40)
  
  # save file
  write.table(data_1, paste(out_dir,"/featureCounts_final.txt",sep=""), sep="\t", 
              quote=FALSE, row.names=FALSE)
  
  return(data_1)
}




### Analysis

feature_counts <- fc_process(feature_counts_file_raw=feature_counts_file_raw, 
                             gtf_file=gtf_file,
                             out_dir=out_dir)







