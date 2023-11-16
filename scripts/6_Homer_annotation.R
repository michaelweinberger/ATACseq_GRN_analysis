
### Run this script after adding gene annotation to peaks using HOMER
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

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("DiffBind")) BiocManager::install("DiffBind")
if (!require("DescTools")) BiocManager::install("DescTools")
if (!require("GenomicRanges")) BiocManager::install("GenomicRanges")
if (!require("ChIPpeakAnno")) BiocManager::install("ChIPpeakAnno")
if (!require("DESeq2")) install.packages("DESeq2")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("stringr")) install.packages("stringr")







### functions

## function to process homer annotated peak file
# homer_output_file: Output of Homer's "annotatePeaks.pl" command
# condition: string that will be used to name gene annotation columns in output file
# out_dir: Directory for output file
homer_process <- function(homer_output_file, condition, out_dir) {
  print(paste("Processing ",homer_output_file,sep=""))

  # read in file and adjust column names
  ann_all <- read.table(homer_output_file, header=TRUE, sep="\t")
  colnames(ann_all)[1] <- "PeakID"
  colnames(ann_all)[colnames(ann_all)=="Annotation"] <- "Feature"
  colnames(ann_all)[colnames(ann_all)=="Distance.to.TSS"] <- "Distance_nearest_TSS"
  colnames(ann_all)[colnames(ann_all)=="Nearest.PromoterID"] <- paste("TranscriptID_", condition, sep="")
  colnames(ann_all)[colnames(ann_all)=="Entrez.ID"] <- paste("gene_id_", condition, sep="")
  colnames(ann_all)[colnames(ann_all)=="Gene.Name"] <- paste("gene_name_", condition, sep="")
  colnames(ann_all)[colnames(ann_all)=="CpG."] <- "CpG"
  colnames(ann_all)[colnames(ann_all)=="GC."] <- "GC"

  # subset data 
  ann_all <- ann_all[,c("PeakID", "Chr", "Start", "End", "Strand",
                        "Feature", "Distance_nearest_TSS",
                        "CpG", "GC",
                        paste("TranscriptID_", condition, sep=""),
                        paste("gene_id_", condition, sep=""),
                        paste("gene_name_", condition, sep=""))]

  # fill empty gene names 
  ann_all[ann_all[,paste("gene_name_", condition, sep="")]=="",
          paste("gene_name_", condition, sep="")] <- "NA"

  # sort data
  ann_all <- ann_all[order(ann_all$PeakID),]
  
  # save data
  print(paste("Saving processed data as ",
              paste(out_dir,"/Diffbind_annotated_peaks_gene_names",sep=""), sep=""))
  saveRDS(ann_all, paste(out_dir,"/Diffbind_annotated_peaks_gene_names_", condition,sep=""))
  
  return(ann_all)
}



## function to plot genomic feature frequencies and GC / CpG content
# homer_annotation_file: File path of output of "homer_process" function
# out_dir: Directory for output plots
plot_feature_freq_gc_cpg <- function(homer_annotation_file, out_dir) {
  print(paste("Plotting genomic feature frequencies and GC,CpG content for file ",
              homer_annotation_file,sep=""))
  
  ann_names <- readRDS(homer_annotation_file)
  
  # plot frequency of genomic feature annotation in entire peak set
  feat_all <- ann_names[,c("PeakID","Feature","Distance_nearest_TSS","CpG","GC","gene_name_all")]
  feat_all$Feature <- as.character(feat_all$Feature)
  feat_all$Feature[grep("^Intergenic",feat_all$Feature)] <- "a_intergenic"
  feat_all$Feature[grep("^intron",feat_all$Feature)] <- "e_intron"
  feat_all$Feature[grep("^exon",feat_all$Feature)] <- "d_exon"
  feat_all$Feature[grep("^promoter-TSS",feat_all$Feature)] <- "b_promoter"
  feat_all$Feature[grep("^TTS",feat_all$Feature)] <- "g_tts"
  feat_all$Feature[grep("^3'",feat_all$Feature)] <- "f_utr3"
  feat_all$Feature[grep("^5'",feat_all$Feature)] <- "c_utr5"
  
  # plot bar charts of genomic feature frequencies and GC / CpG contents
  feat_all_count <- feat_all %>% group_by(Feature) %>% summarise(count=n())
  feat_all_count$gc <- as.data.frame(feat_all %>% group_by(Feature) %>% summarise(mean=mean(GC)))[,2]
  feat_all_count$cpg <- as.data.frame(feat_all %>% group_by(Feature) %>% summarise(mean=mean(CpG)))[,2]
  feat_all_count$freq <- feat_all_count$count / sum(feat_all_count$count) * 100
  
  gg <- ggplot(feat_all_count, aes(x=Feature,y=freq)) +
  geom_bar(aes(fill=Feature),stat="identity", width=0.8) +
  scale_fill_manual(breaks=c("a_intergenic","b_promoter","c_utr5","d_exon","e_intron","f_utr3","g_tts"),
    name="",
    labels=c("Intergenic","Promoter","5'UTR","Exon","Intron","3'UTR","TTS"),
    values=c("darkblue", "darkgreen", "dodgerblue", "firebrick2", "maroon4", "darkorange", "saddlebrown")) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x = element_text(size = 22, colour = "black",face = "plain", angle=30, hjust=1),
    axis.text.y = element_text(size = 22, colour = "black",face = "plain"),
    axis.title.y = element_text(size = 22, face = "plain", vjust = 0),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "right",
    legend.text = element_text(size = 18, colour = "black"),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    strip.text.x = element_text(size=16, face="plain"),
    strip.background = element_blank()) +
  scale_y_continuous(name="Relative number of peaks [%]", limits=c(0, 60), breaks=seq(0,60,by=10)) +
  scale_x_discrete(name="", labels=c("Intergenic","Promoter","Exon","Intron","TTS")) +
  labs(title="")
  
  ggsave(gg, filename=paste(out_dir,"/Diffbind_annotated_peaks_gene_names_genomic_feat_freq.pdf",sep=""),
         device="pdf", dpi=300, useDingbats=FALSE)

  # plot gc and cpg content as violin plots
  gg <- ggplot(feat_all, aes(x=Feature,y=GC*100)) +
  geom_violin(aes(fill=Feature), trim = FALSE) +
  scale_fill_manual(breaks=c("a_intergenic","b_promoter","c_utr5","d_exon","e_intron","f_utr3","g_tts"),
    name="",
    labels=c("Intergenic","Promoter","5'UTR","Exon","Intron","3'UTR","TTS"),
    values=c("darkblue", "darkgreen", "dodgerblue", "firebrick2", "maroon4", "darkorange", "saddlebrown")) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x = element_text(size = 22, colour = "black",face = "plain", angle=30, hjust=1),
    axis.text.y = element_text(size = 22, colour = "black",face = "plain"),
    axis.title.y = element_text(size = 22, face = "plain", vjust = 0),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "right",
    legend.text = element_text(size = 18, colour = "black"),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    strip.text.x = element_text(size=16, face="plain"),
    strip.background = element_blank()) +
  scale_y_continuous(name="GC content [%]", limits=c(0, 100), breaks=seq(0,100,by=20)) +
  scale_x_discrete(name="", labels=c("Intergenic","Promoter","Exon","Intron","TTS")) +
  labs(title="") +
  geom_boxplot(width = 0.2)
  
  ggsave(gg, filename=paste(out_dir,"/Diffbind_annotated_peaks_gene_names_genomic_feat_gc.pdf",sep=""),
         device="pdf", dpi=300, useDingbats=FALSE)
  
  gg <- ggplot(feat_all, aes(x=Feature,y=CpG*100)) +
  geom_violin(aes(fill=Feature), trim = FALSE) +
  scale_fill_manual(breaks=c("a_intergenic","b_promoter","c_utr5","d_exon","e_intron","f_utr3","g_tts"),
    name="",
    labels=c("Intergenic","Promoter","5'UTR","Exon","Intron","3'UTR","TTS"),
    values=c("darkblue", "darkgreen", "dodgerblue", "firebrick2", "maroon4", "darkorange", "saddlebrown")) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x = element_text(size = 22, colour = "black",face = "plain", angle=30, hjust=1),
    axis.text.y = element_text(size = 22, colour = "black",face = "plain"),
    axis.title.y = element_text(size = 22, face = "plain", vjust = 0),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "right",
    legend.text = element_text(size = 18, colour = "black"),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    strip.text.x = element_text(size=16, face="plain"),
    strip.background = element_blank()) +
  scale_y_continuous(name="CpG content [%]", limits=c(0, 20), breaks=seq(0,20,by=4)) +
  scale_x_discrete(name="", labels=c("Intergenic","Promoter","Exon","Intron","TTS")) +
  labs(title="") +
  geom_boxplot(width = 0.2)
  
  ggsave(gg, filename=paste(out_dir,"/Diffbind_annotated_peaks_gene_names_genomic_feat_cpg.pdf",sep=""),
         device="pdf", dpi=300, useDingbats=FALSE)
  
  return(feat_all)
}



## function to perform TSS enrichment analysis
# homer_tss_file: File path to output of Homer's "annotatePeaks.pl tss" command
# conditions: Names of tag directories used in Homer analysis
# out_dir: Directory for output plot
plot_tss_enrichment <- function(homer_tss_file, conditions, out_dir) {
  print(paste("Plotting TSS accessibility enrichment using file ",
              homer_tss_file,sep=""))
  
  # read in Homer output
  tss <- read.table(homer_tss_file, header=TRUE, sep="\t")
  cols2keep <- c(1, seq(2,ncol(tss) - 1,3))
  tss <- tss[,cols2keep]
  colnames(tss)[1] <- "Distance"

  # split colnames over dots
  num_split <- str_count(colnames(tss)[2], "\\.")
  colnames(tss)[2:ncol(tss)] <- matrix(unlist(strsplit(colnames(tss)[2:ncol(tss)], "\\.")), 
                                       byrow=TRUE, ncol=(num_split+1))[,num_split]
  
  # check which columns in tss belong to which condition +
  # compute mean coverage values for each condition
  for (i in seq(1,length(conditions))) {
    tmp <- grep(conditions[i], colnames(tss))
    name <- paste(conditions[i],"_cov",sep="")
    if (length(tmp) == 1) {
      tss[,name] <- tss[,tmp]
    } else if (length(tmp) > 1) {
      tss[,name] <- rowMeans(tss[,tmp])
    }
  }
  
  # generate dataframe for ggplot
  start <- ncol(tss)-(length(conditions)-1)
  for (i in seq(start,ncol(tss))) {
    tmp <- cbind(tss[,c(1,i)], rep(colnames(tss)[i],nrow(tss)))
    colnames(tmp) <- c("distance","coverage", "condition")
    if (i == start) {
      tss_plot <- tmp
    } else {
      tss_plot <- rbind(tss_plot, tmp)
    }
  }
  tss_plot$distance <- as.factor(tss_plot$distance)
  
  # generate plot
  max_value <- ceiling(max(tss_plot$coverage))
  colours <- c("darkorange", "saddlebrown", "maroon2", "darkmagenta", "green3", "darkgreen", "black")
  
  gg <- ggplot(tss_plot, aes(x=distance,y=coverage, group=condition)) +
  geom_line(aes(color=condition), size=1.2) +
  scale_color_manual(values=colours[1:length(conditions)]) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x = element_text(size = 22, colour = "black",face = "plain", angle=30, hjust=1),
    axis.text.y = element_text(size = 22, colour = "black",face = "plain"),
    axis.title.y = element_text(size = 22, face = "plain", vjust = 0),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "right",
    legend.text = element_text(size = 18, colour = "black"),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    strip.text.x = element_text(size=16, face="plain"),
    strip.background = element_blank()) +
  scale_y_continuous(name="Fragment depth per bp per peak", limits=c(0, max_value), breaks=seq(0,max_value,by=1)) +
  scale_x_discrete(name="Distance from TSS", breaks=c("-2000","-1000","0","1000","2000")) +
  labs(title="")
  
  ggsave(gg, filename=paste(out_dir,"/Diffbind_annotated_peaks_gene_names_tss_enrichment.pdf",sep=""),
         device="pdf", dpi=300, useDingbats=FALSE)

  return(tss_plot)
}






### Analysis

conditions <- unlist(read.table(paste(conditions_dir, "/conditions.txt", sep=""), sep=" "))



## after homer annotation analysis, process annotated peak file
# annotated using all genes
homer_output_file <- paste(homer_dir, "/Diffbind_consensus_peak_set_annotated_", 
                           genome, "_subset.txt", sep="")
homer_annotation <- homer_process(homer_output_file=homer_output_file, condition="all",
                                  out_dir=out_dir)



## OPTIONAL: add gene annotation using expressed genes in individual conditions
# collect expressed gene files
homer_output_file_list <- list.files(path=homer_dir, 
                                     pattern="Diffbind_consensus_peak_set_annotated_Expressed_genes")

if (length(homer_output_file_list) > 0) {
  for (file in homer_output_file_list) {

    # build file path
    homer_output_file <- paste(homer_dir, file, sep="/")

    # extract condition from file name
    condition <- sub("Diffbind_consensus_peak_set_annotated_Expressed_genes_", "", file)
    condition <- sub("_rpkm_cutoff_.*", "", condition)

    # process file
    homer_annotation_tmp <- homer_process(homer_output_file=homer_output_file, condition=condition,
                                          out_dir=out_dir)

    # combine output of all-gene annotation with relevant columns of expressed-genes annotation
    homer_annotation_tmp_1 <- merge(homer_annotation, 
                                    homer_annotation_tmp[,c("PeakID",
                                                            paste("TranscriptID_", condition, sep=""),
                                                            paste("gene_id_", condition, sep=""),
                                                            paste("gene_name_", condition, sep=""))],
                                    by="PeakID", all.x=TRUE)
    homer_annotation <- homer_annotation_tmp_1
  }
}



## save final file
saveRDS(homer_annotation, 
        paste(out_dir, "/Diffbind_annotated_peaks_gene_names", sep=""))



## plot genomic feature frequencies and GC / CpG content
homer_annotation_file <- paste(out_dir,"/Diffbind_annotated_peaks_gene_names",sep="")
feature_freq_gc_cpg <- plot_feature_freq_gc_cpg(homer_annotation_file=homer_annotation_file, 
                                                out_dir=out_dir)


## Plot output of Homer TSS enrichment analysis
homer_tss_file <- paste(homer_dir,"/Tag_density_tss.txt",sep="")
tss_enrichment <- plot_tss_enrichment(homer_tss_file=homer_tss_file, conditions=conditions, 
                                      out_dir=out_dir)














