
### Run this script after adding gene annotation to peaks using HOMER
### to perform DESeq2 differential peak accessibility analyses




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

if (!require("DiffBind")) install.packages("DiffBind")
if (!require("DESeq2")) install.packages("DESeq2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("rlist")) install.packages("rlist")




### functions

## subfunction to extract significantly enriched peaks and annotate to closest expressed gene
annotate_peaks <- function(peak_df, annotation_df, annotation_col) {
  peak_df$PeakID <- rownames(peak_df)
  peak_df <- merge(peak_df, annotation_df[,c("PeakID","Feature","CpG","GC",annotation_col)], by="PeakID")
  peak_df[,annotation_col] <- make.names(peak_df[,annotation_col], unique=TRUE)
  peak_df$Peak_gene_ID <- paste(peak_df$PeakID, peak_df[,annotation_col], sep="__")
  rownames(peak_df) <- peak_df$Peak_gene_ID
  peak_df <- peak_df[order(peak_df$padj, decreasing=FALSE),]
  return(peak_df)
}



## function to run DESeq2 comparison of samples in condition_1 versus samples in condition_2,
# as identified by colnames
# batch_list (optional): list of integers defining batch numbers of samples
deseq2_run <- function(count_df, rpkm_df, annotation_df, annotation_col_1, annotation_col_2,
                       condition_1, condition_2, batch_list, volcano_label="top_diff", 
                       n_top_diff=10, n_peaks_heatmap=10, out_dir) {
  for (i in seq(1,length(condition_1))) {
    tmp <- grep(condition_1[i], colnames(count_df))
    if (i == 1) {
      colpos_1 <- tmp
    } else {
      colpos_1 <- list.append(colpos_1, tmp)
    }   
  }
  
  for (i in seq(1,length(condition_2))) {
    tmp <- grep(condition_2[i], colnames(count_df))
    if (i == 1) {
      colpos_2 <- tmp
    } else {
      colpos_2 <- list.append(colpos_2, tmp)
    }
  }
  
  # generate sample metadata  
  condition_list_1 <- paste(condition_1, collapse="_")
  condition_list_2 <- paste(condition_2, collapse="_")
  sample_info <- as.data.frame(c(rep(condition_list_1, length(colpos_1)), 
                                 rep(condition_list_2, length(colpos_2))))
  
  if (missing(batch_list)) {
    colnames(sample_info) <- "condition"
  } else {
    sample_info$batch <- as.factor(batch_list)
    colnames(sample_info) <- c("condition", "batch")
  }
  rownames(sample_info) <- c(colnames(count_df)[colpos_1], colnames(count_df)[colpos_2])
  print(sample_info)
  
  # construct DESeq object
  tmp <- cbind(as.data.frame(count_df[,colpos_1]), as.data.frame(count_df[,colpos_2]))
  colnames(tmp) <- c(colnames(count_df)[colpos_1], colnames(count_df)[colpos_2])
  
  if (missing(batch_list)) {
    dds <- DESeqDataSetFromMatrix(countData=tmp, colData=sample_info, design= ~ condition)
  } else {
    dds <- DESeqDataSetFromMatrix(countData=tmp, colData=sample_info, design= ~ batch + condition)  
  }
  rownames(dds) <- rownames(count_df)
  
  # filter out rows with only 0 or 1 counts
  dds <- dds[rowSums(counts(dds)) > 1,]
  
  # ensure correct comparison
  dds$condition <- relevel(dds$condition, ref=condition_list_2)
  #print(dds)
  
  # run analysis
  dds <- DESeq(dds)
  res <- results(dds, alpha=0.05)
  res$padj[is.na(res$padj)] <- 1
  print(summary(res))
  saveRDS(res, paste(out_dir, "/DESeq2_ATAC_", condition_list_1, "_vs_", condition_list_2, sep=""))
  
  # extract significantly enriched DAPs and annotate to closest expressed gene
  resSig_up <- as.data.frame(subset(res, log2FoldChange >= 0))
  resSig_up <- annotate_peaks(resSig_up, annotation_df, annotation_col_1)
  write.csv(resSig_up, file = paste(out_dir, "/DESeq2_ATAC_", 
                                    condition_list_1, "_vs_", condition_list_2, "__", 
                                    condition_list_1, "_enriched.csv", sep=""), row.names=TRUE)
  resSig_down <- as.data.frame(subset(res, log2FoldChange < 0))
  resSig_down <- annotate_peaks(resSig_down, annotation_df, annotation_col_2)
  write.csv(resSig_down, file = paste(out_dir, "/DESeq2_ATAC_", 
                                      condition_list_1, "_vs_", condition_list_2, "__", 
                                      condition_list_2, "_enriched.csv", sep=""), row.names=TRUE)
  
  colnames(resSig_up)[ncol(resSig_up)-1] <- "gene_name"
  colnames(resSig_down)[ncol(resSig_down)-1] <- "gene_name"
  resSig_up_down <- rbind(resSig_up, resSig_down)
  
  
  # volcano plot
  resName <- data.frame(cbind(rownames(res), res))
  colnames(resName)[1] <- "PeakID"
  resName$up_down <- resName$log2FoldChange / abs(resName$log2FoldChange)
  resName <- resName[order(resName$up_down,resName$padj,decreasing=FALSE),]
  resName$Significant <- ifelse(resName$padj < 0.05, "Adjusted p-value < 0.05", "Not Sig")
  
  # generate annotation labels
  if (volcano_label=="top_diff") {
    up_start <- length(which(resName$up_down == -1)) + 1
    label_data <- rbind(resName[1:n_top_diff,],resName[up_start:(up_start + (n_top_diff-1)),])
    label_data$gene_name <- rep(1,nrow(label_data))
    for (i in 1:n_top_diff) {
      tmp <- resSig_down[grep(paste(label_data[i,"PeakID"], "_", sep=""), 
                              paste(resSig_down[,"PeakID"], "_", sep="")), "gene_name"]
      if (length(tmp) > 0) {
        label_data[i,"gene_name"] <- tmp
      }
    }
    for (i in (n_top_diff+1):nrow(label_data)) {
      tmp <- resSig_up[grep(paste(label_data[i,"PeakID"], "_", sep=""), 
                            paste(resSig_up[,"PeakID"], "_", sep="")), "gene_name"]
    if (length(tmp) > 0) {
        label_data[i,"gene_name"] <- tmp
      }
    }
    label_data <- label_data[label_data$gene_name!=1,]
  } else {
    for (i in 1:length(volcano_label)) {
      tmp <- resSig_up_down[grep(paste("_", volcano_label, "_", sep=""), 
                                 paste("_", resSig_up_down$gene_name, "_", sep="")),]
      if (i==1) {
        label_data <- tmp
        
      } else {
        label_data <- rbind(label_data, tmp)
      }
    }
    label_data$Peak_gene_ID <- rownames(label_data)
  }
  #print(label_data)
  
  limit_max <- as.integer(max(-log10(resName$padj))) + 10
  limit_max_x <- as.integer(max(abs(resName$log2FoldChange))) + 5
  limit_max_x <- round(limit_max_x / 5) * 5
  title <- paste(condition_list_1, "vs", condition_list_2, sep="_")
  
  # plot dots only + save as png
  ggplot(resName, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(size = 20, colour = "black", face = "plain", angle=0, hjust=0.5),
    axis.text.y = element_text(size = 20, colour = "black", face = "plain"),
    axis.title.x = element_text(size = 20, vjust = -1, hjust=0.5),
    axis.title.y = element_text(size = 20, vjust = 1),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "none",
    plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")
  ) +
  scale_y_continuous(name="Adjusted p-value [-log10]", 
                     limits=c(0,limit_max), 
                     breaks=seq(0, limit_max, by=10)) +
  scale_x_continuous(name="Fold change [log2]", limits=c(-limit_max_x, limit_max_x), 
                     breaks=seq(-limit_max_x, limit_max_x, by=5)) +
  ggtitle(title)
  
  ggsave(filename=paste(out_dir, "/DESeq2_ATAC_", condition_list_1, "_vs_", 
                        condition_list_2, "_volcano_plot_dots_only.png", sep=""), 
         device="png", dpi=300)
  
  
  # plot with gene annotation + save as pdf
  ggplot(resName, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(size = 20, colour = "black", face = "plain", angle=0, hjust=0.5),
    axis.text.y = element_text(size = 20, colour = "black", face = "plain"),
    axis.title.x = element_text(size = 20, vjust = -1, hjust=0.5),
    axis.title.y = element_text(size = 20, vjust = 1),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "none",
    plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")
  ) +
  scale_y_continuous(name="Adjusted p-value [-log10]", 
                     limits=c(0,limit_max), 
                     breaks=seq(0, limit_max, by=10)) +
  scale_x_continuous(name="Fold change [log2]", limits=c(-limit_max_x, limit_max_x), 
                     breaks=seq(-limit_max_x, limit_max_x, by=5)) +
  ggtitle(title) +
  geom_text_repel(
    data = label_data,
    aes(label = gene_name),
    size = 7,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = Inf
  )
  
  ggsave(filename=paste(out_dir, "/DESeq2_ATAC_", condition_list_1, "_vs_", 
                        condition_list_2, "_volcano_plot.pdf", sep=""), 
         device="pdf", dpi=300, useDingbats=FALSE)
  
  
  # heatmap of top peaks enriched in either condition
  label_data <- rbind(resSig_up_down[1:n_peaks_heatmap,], 
                      resSig_up_down[(nrow(resSig_up)+1):(nrow(resSig_up)+n_peaks_heatmap),])
  
  log_rpkm_count_df <- as.data.frame(log2(rpkm_df + 1))
  # combine peak ids and annotated gene names + add as row names
  log_rpkm_count_df$PeakID <- rownames(rpkm_df)
  log_rpkm_count_df <- merge(log_rpkm_count_df, resSig_up_down[,c("PeakID", "Peak_gene_ID")], by="PeakID")
  rownames(log_rpkm_count_df) <- make.names(log_rpkm_count_df$Peak_gene_ID, unique=TRUE)
  log_rpkm_count_df$Peak_gene_ID <- NULL
  log_rpkm_count_df$PeakID <- NULL
  
  pdf(file = paste(out_dir, "/DESeq2_ATAC_", condition_list_1, "_vs_", 
                   condition_list_2, "_heatmap_top.pdf", sep=""))
    pheatmap(log_rpkm_count_df[label_data$Peak_gene_ID,], 
    cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, show_colnames=TRUE,
    fontsize = 10, color=colorRampPalette(rev(brewer.pal(n = 7, name ="PuOr")))(100))
  dev.off()
  
  return(res)
}




### Analysis

conditions <- unlist(read.table(paste(conditions_dir, "/conditions.txt", sep=""), sep=" "))

## Load pre-constructed consensus peak set + gene annotation
all_data_ATAC_all <- readRDS(paste(out_dir,"/Diffbind_consensus_peak_set",sep=""))
all_data_ATAC_all$PeakID <- rownames(all_data_ATAC_all)
all_data_ATAC_all <- all_data_ATAC_all[,c(ncol(all_data_ATAC_all),1:(ncol(all_data_ATAC_all)-1))]


# extract raw counts to compare with DESeq2
n_samples <- (ncol(all_data_ATAC_all)-4) / 2
count_data_ATAC_all <- as.data.frame(all_data_ATAC_all[,5:(4+n_samples)])
rownames(count_data_ATAC_all) <- rownames(all_data_ATAC_all)

# extract FPKM counts
rpkm_data_ATAC_all <- as.data.frame(all_data_ATAC_all[,(5+n_samples):ncol(all_data_ATAC_all)])
rownames(rpkm_data_ATAC_all) <- rownames(all_data_ATAC_all)


ann_names <- readRDS(paste(out_dir,"/Diffbind_annotated_peaks_gene_names",sep=""))



## Loop over conditions and perform pairwise DESeq2 comparisons
# only run if there are at least 2 conditions
if (length(conditions) > 1) {

  # Initialise list with conditions that remain to be compared
  conditions_2 <- conditions

  # do not loop over last condition, because it would not have anything to compare it to
  for (i in seq(1,(length(conditions)-1))) {

    # identify column containing gene annotation for condition 1 
    # check if expressed-genes peak annotation exists for condition 1
    if(paste("gene_name_", conditions[i], sep="") %in% colnames(ann_names)) {
      annotation_col_1 <- paste("gene_name_", conditions[i], sep="")
    } else {
      annotation_col_1 <- paste("gene_name_", "all", sep="")
    }

    # gradually shrink conditions_2 list to avoid duplicate comparisons
    conditions_2 <- conditions_2[conditions_2 != conditions[i]]

    for (k in conditions_2) {

      # identify column containing gene annotation for condition 2 
      # check if expressed-genes peak annotation exists for condition 2
      if(paste("gene_name_", k, sep="") %in% colnames(ann_names)) {
        annotation_col_2 <- paste("gene_name_", k, sep="")
      } else {
        annotation_col_2 <- paste("gene_name_", "all", sep="")
      }

      # run differential gene expression analysis
      res <- deseq2_run(count_df=count_data_ATAC_all, rpkm_df=rpkm_data_ATAC_all, 
                        annotation_df=ann_names, 
                        annotation_col_1=annotation_col_1,
                        annotation_col_2=annotation_col_2, 
                        condition_1=conditions[i], condition_2=k, 
                        volcano_label="top_diff", n_top_diff=10, n_peaks_heatmap=10,
                        out_dir=out_dir)
    }
  }
}




