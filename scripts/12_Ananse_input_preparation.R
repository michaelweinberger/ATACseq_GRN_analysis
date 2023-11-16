### Run this script to prepare input data for Ananse analysis



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

peak_width <- as.numeric(peak_width)
chrom_number <- as.numeric(chrom_number)

conditions_ATAC <- unlist(read.table(paste(out_dir, "/conditions_ananse_ATAC.txt", sep=""), sep=" "))
conditions_RNA <- unlist(read.table(paste(out_dir, "/conditions_ananse_RNA.txt", sep=""), sep=" "))




### packages

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("scater")) BiocManager::install("scater")
if (!require("edgeR")) install.packages("edgeR")
if (!require("DESeq2")) BiocManager::install("DESeq2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("TFBSTools")) BiocManager::install("TFBSTools")
if (!require("GenomicFeatures")) BiocManager::install("GenomicFeatures")
if (!require("rtracklayer")) BiocManager::install("rtracklayer")




### Functions

## subfunction to generate TPM values
to_tpm <- function(read_count, gene_length) {
  gene_length_kb <- gene_length / 1000
  rpk <- read_count / gene_length_kb
  total_counts <- sum(rpk) / 1e6
  tpm <- t(t(rpk) / total_counts)
  return(tpm)
}



## function to generate TPM values of bulk RNAseq data 
tpm_bulk <- function(feature_counts_file, condition, out_dir) {

  # read in bulk RNA-seq data
  all_data <- read.delim(feature_counts_file, header=TRUE, sep="\t", 
                         stringsAsFactors=FALSE, check.names=FALSE)
  
  all_data <- all_data[grep("ENS",all_data$gene_id),]
  rownames(all_data) <- make.names(all_data$gene_name, unique=TRUE)

  # subset read count columns to those belonging to condition of interest
  cols2keep <- colnames(all_data)[grep(condition,colnames(all_data))]
  all_data <- all_data[,c("gene_id", "gene_name", "gene_length", cols2keep)]
  #print(paste("RNA-seq raw counts of ",condition,":",sep=""))
  #print(all_data, max=30)

  # identify columns with count data (all columns except gene_id, gene_name, gene_length columns)
  count_cols <- colnames(all_data)[!colnames(all_data) %in% c("gene_id", "gene_name", "gene_length")]
  print(paste("Columns with count data for ", condition, ": ", count_cols, sep=""))

  # convert raw count values to TPM and save TPM values for each column separately
  for (col in count_cols) {
    all_data[,col] <- to_tpm(read_count=all_data[,col], gene_length=all_data$gene_length)
    write.table(all_data[,c("gene_name",col)], file=paste(out_dir,"/",col,"_ANANSE_TPM.txt",sep=""), 
                row.names=FALSE, col.names=c("target_id", "tpm"), quote = FALSE, sep="\t")
  }

  #print(paste("RNA-seq TPM counts of ",condition,":",sep=""))
  #print(all_data, max=40)

  # generate mean TPM of all samples for condition of interest
  #all_data_1 <- cbind(all_data[,"gene_name"], as.data.frame(rowMeans(all_data[,count_cols])))
  #write.table(all_data_1, file=paste(out_dir,"/",condition,"_ANANSE_mean_TPM.txt",sep=""), 
  #            row.names=FALSE, col.names=c("target_id", "tpm"), quote = FALSE, sep="\t")

  return(all_data)
}



## function to generate gene position bed file for Ananse Network
ananse_gene_pos <- function(gtf_file, genome, out_dir) {

  # create a mapping of ensembl gene id to gene name
  gtf_df <- as.data.frame(rtracklayer::import(gtf_file))
  genes_names <- unique(gtf_df[,c("gene_id","gene_name")])

  # extract named genes
  gtf_df_genes <- gtf_df[gtf_df$type=="gene",]
  gtf_df_genes <- gtf_df_genes[!is.na(gtf_df_genes$gene_name),]
  gtf_df_genes$gene_name <- make.names(gtf_df_genes$gene_name, unique=TRUE)

  # exclude rows with gene name containing "unm_" or "NC_"
  gtf_df_genes <- gtf_df_genes[!grepl("unm_*",gtf_df_genes$gene_name),]
  gtf_df_genes <- gtf_df_genes[!grepl("NC_*",gtf_df_genes$gene_name),]

  # replace underscore in gene names with period
  gtf_df_genes$gene_name <- gsub("_", "\\.", gtf_df_genes$gene_name)

  # keep only gene name, not gene id
  gtf_df_genes <- gtf_df_genes[,c("seqnames","start","end","gene_name","width","strand")]


  # add dummy columns to fit 12-column format required by Ananse
  x <- ncol(gtf_df_genes)
  for (i in seq(1,6)) {
    colname <- paste("col",(i+x),sep="_")
    gtf_df_genes[,colname] <- rep(0,nrow(gtf_df_genes))
  } 

  write.table(gtf_df_genes, paste(out_dir,"/ANANSE_",genome,"_gene_positions.bed",sep=""), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  return(gtf_df_genes)
}






### Analysis

print("Performing Ananse input preparation")

# generate single RNAseq condition TPM value files
for (condition in conditions_RNA) {
  tpm_values <- tpm_bulk(feature_counts_file=feature_counts_file, 
                         condition=condition, out_dir=out_dir)
}

# generate gene position bed file for Ananse Network
gene_pos <- ananse_gene_pos(gtf_file=gtf_file, genome=genome, out_dir=out_dir)







