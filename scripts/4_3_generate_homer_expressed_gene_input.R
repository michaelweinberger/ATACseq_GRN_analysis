### Run this script to generate .gtf files that have been subset to expressed genes,
# as defined by input featureCounts file
# was run in R 4.0.3




### user defined variables

# unpack variables passed from parent shell script
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

for (e in args) {
  argname <- e[1]
  argval <- e[2]
  # regular expression to delete initial \" and trailing \"
  argval <- gsub("(^\\\"|\\\"$)", "", argval)
  assign(argname, argval)
}




### packages

if (!require("biomaRt")) BiocManager::install("biomaRt")
if (!require("edgeR")) BiocManager::install("edgeR")




### Functions

## function to extract genes expressed above cutoff in samples of specified condition
# feature_counts_file: Location of processed FeatureCounts output file, file should contain 
#	as first three columns: "gene_id", "gene_name" and "gene_length"
# condition: Part of column name(s) identifying samples belonging to the same condition
# rpkm_cutoff: Minimum RPKM value for a gene to be considered expressed
# out_dir: Directory for output text file listing gene IDs of expressed genes
extract_expressed_genes <- function(feature_counts_file, condition, rpkm_cutoff=2.5, out_dir) {

  all_data <- read.table(feature_counts_file, header=TRUE, sep="\t")

  # generate RPKM values
  rpkm_meta <- cbind(as.data.frame(all_data$gene_name), as.data.frame(all_data$gene_length))
  colnames(rpkm_meta) <- c("gene_name", "gene_length")
  y_norm <- DGEList(counts=all_data[,4:ncol(all_data)], genes=rpkm_meta) 
  rpkm_data <- rpkm(y_norm, y_norm$genes$gene_length) 
  rownames(rpkm_data) <- all_data$gene_id

  # subset data to samples of condition
  rpkm_tmp <- as.data.frame(rpkm_data[,grep(condition, colnames(rpkm_data))])
  rownames(rpkm_tmp) <- rownames(rpkm_data)
  print(head(rpkm_tmp))

  # extract expressed genes
  if (ncol(rpkm_tmp) == 1) {
    expressed_genes <- rownames(rpkm_tmp[rpkm_tmp > rpkm_cutoff,])
  } else if (ncol(rpkm_tmp) > 1) {
    expressed_genes <- rownames(rpkm_tmp[rowMeans(rpkm_tmp) > rpkm_cutoff,])
  } else {
    print(paste("Error: Condition", condition, "not found in RPKM data"))
  }

  write.table(expressed_genes, 
              file=paste(out_dir, "/Expressed_genes_", condition, "_rpkm_cutoff_", rpkm_cutoff, ".txt", sep=""), 
              row.names=FALSE, col.names=FALSE, quote = FALSE)

  return(expressed_genes)
}




### Analysis

conditions <- unlist(read.table(paste(conditions_dir, "/conditions.txt", sep=""), sep=" "))

# extract expressed genes
for (i in seq(1,length(conditions))) {
  genes_expressed <- extract_expressed_genes(feature_counts_file=feature_counts_file, 
                                             condition=conditions[i], rpkm_cutoff=2.5, 
                                             out_dir=out_dir)
}











