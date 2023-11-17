
### Run this script after mapping ATAC-seq data and calling peaks with MACS




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






### packages

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("DiffBind")) BiocManager::install("DiffBind")
if (!require("DescTools")) BiocManager::install("DescTools")
if (!require("dplyr")) BiocManager::install("dplyr")
if (!require("GenomicRanges")) BiocManager::install("GenomicRanges")
if (!require("ChIPpeakAnno")) BiocManager::install("ChIPpeakAnno")






### functions

## subfunction to exclude peaks overlapping repetitive elements in the genome
# peak_file: MACS .narrowPeak file
# repeats_gr: Genomic ranges object specifying repetetive regions in the genome
# chrom_number: Number of main chromosomes
# cutoff: Percentage of sequence overlap with repetitive elements above which peak will be excluded
F_peak_repeat_cleanup <- function (peak_file, repeats_gr, chrom_number, cutoff=50) {

  print(paste("Removing repetitive regions from ",peak_file,sep=""))

  ## read in peaks file
  peaks <- read.table(peak_file)
  colnames(peaks) <- c("chr","start","end","name","score_1","strand",
                       "score_2","score_3","score_4","score_5")
  peaks$name <- make.names(peaks$name)

  # subset to main chromosomes
  peaks <- peaks[peaks$chr %in% seq(1:chrom_number),]
  
  peaks_gr <- toGRanges(peaks, format="BED")

  ## intersect peak positions and repetitive element positions
  tmp <- GenomicRanges::intersect(peaks_gr, repeats_gr)
  tmp$width <- width(ranges(tmp))

  # find overlaps between peaks and peak intersections with repetitive elements
  tmp_2 <- findOverlaps(peaks_gr, tmp)
  peaks_ov <- peaks_gr[from(tmp_2),]
  peak_names <- as.data.frame(names(ranges(peaks_ov)))
  names(peaks_ov) <- seq(1,length(peaks_ov))
  peaks_ov <- as.data.frame(ranges(peaks_ov))
  peaks_ov$peak <- queryHits(tmp_2)
  width_ov <- as.data.frame(tmp[to(tmp_2),c("width")])

  # combine peaks and overlaps
  peaks_width_ov <- cbind(peak_names, peaks_ov, width_ov[,"width"])
  colnames(peaks_width_ov)[ncol(peaks_width_ov)] <- "width_ov"

  # if peak overlaps multiple repetitive elements, sum those overlaps up
  peaks_width_ov <- peaks_width_ov %>% group_by(peak) %>% 
                    mutate(width_ov_total=sum(width_ov)) %>% ungroup()
  peaks_width_ov$perc <- (peaks_width_ov$width_ov_total / peaks_width_ov$width) * 100
  peaks2drop <- as.data.frame(peaks_width_ov[peaks_width_ov$perc > cutoff,], stringsAsFactors=FALSE)

  # subset peak file
  peaks2keep <- peaks[!peaks$name %in% peaks2drop[,1],]

  print(paste("Number of retained peaks: ",nrow(peaks2keep),sep=""))
  return(peaks2keep)
}



## function to generate DiffBind peak set file
# bam_dir: Directory containing input ATAC-seq bam files
# peak_dir: Directory containing MACS .narrowPeak files
# out_dir: Directory containing output rds file
# repeats_file: File specifying repetetive regions in the genome
# peak_width: Length of peaks in consensus peak set
# chr_sizes: File specifying chromosome lengths
# chrom_number: Number of main chromosomes
diffbind <- function(bam_dir, peak_dir, out_dir, repeats_file,
                     peak_width, chr_sizes, chrom_number) {
  
  # define sample input information
  bam_list <- list.files(path=bam_dir, pattern=".rmdup.bam", full.names=TRUE)
  bam_list <- bam_list[grep("*rmdup.bam_", paste(bam_list,"_",sep=""))]
  print("Using following BAM files:")
  print(bam_list)
  
  peak_list <- list.files(path=peak_dir, pattern=".narrowPeak", full.names=FALSE)
  peak_list_full <- list.files(path=peak_dir, pattern=".narrowPeak", full.names=TRUE)
  print("Using following peak BED files:")
  print(peak_list_full)
  
  samples <- rep(1,length(peak_list))
  for (i in seq(1,length(peak_list))) {
    samples[i] <- matrix(unlist(strsplit(peak_list[i], "\\.")), ncol = 6, byrow = TRUE)[1]
  }
  samples <- gsub("_peaks", "", samples)
  print("Sample names used:")
  print(samples)
  
  # process repeat file
  repeats <- read.table(repeats_file)
  repeats <- repeats[,c(6:8,11:13)]
  colnames(repeats) <- c("chr","start","end","name","type","subtype")
  
  # subset to main chromosomes
  repeats <- repeats[repeats$chr %in% paste("chr",seq(1:chrom_number),sep=""),]
  
  # change chromosome names to Ensembl format
  repeats$chr <- gsub("chr", "", repeats$chr)
  
  # transform into granges and merge overlapping elements
  repeats_gr <- toGRanges(repeats, format="BED")
  repeats_gr <- reduce(repeats_gr)
  print(head(repeats_gr))
  
  # exclude repetitive peaks from narrowPeak files
  no_repeats_dir <- paste(peak_dir, "/No_repetitive_elements", sep="")
  if (!dir.exists(no_repeats_dir)) {
    dir.create(no_repeats_dir)
  }
  
  for (i in seq(1,length(peak_list))) {
    peak_file <- peak_list[i]
    peak_file_full <- peak_list_full[i]
    cleaned_peaks <- F_peak_repeat_cleanup(peak_file_full, repeats_gr, chrom_number=chrom_number)
    write.table(cleaned_peaks, file=paste(no_repeats_dir, "/No_repetitive_elements_", peak_file, sep=""),
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }
  
  no_rep_peak_list <- list.files(path=paste(peak_dir, "/No_repetitive_elements", sep=""), 
                               pattern=".narrowPeak", full.names=TRUE)
  
  # generate .csv samplesheet containing sample data for DiffBind
  samplesheet <- data.frame(samples)
  colnames(samplesheet) <- "SampleID"
  samplesheet$bamReads <- bam_list
  samplesheet$Peaks <- no_rep_peak_list
  samplesheet$PeakCaller <- rep("bed", length(no_rep_peak_list))
  write.csv(samplesheet, file=paste(out_dir,"/Diffbind_sample_sheet.csv",sep=""), row.names = FALSE)
  
  # Construct DiffBind consensus peak set
  print("Constructing consensus peak set")
  
  # generate dba object for analysis
  ATAC_all <- dba(sampleSheet=paste(out_dir,"/Diffbind_sample_sheet.csv",sep=""))
  
  pdf(file=paste(out_dir,"/Diffbind_peak_set_correlation_before_counting.pdf",sep=""))
    plot(ATAC_all)
  dev.off()
  
  # read in BAM files and count relevant tags, recenter peaks 
  summits <- peak_width / 2
  ATAC_all <- dba.count(ATAC_all, summits=summits, score="DBA_SCORE_READS")
  
  pdf(file=paste(out_dir,"/Diffbind_peak_set_correlation_after_counting.pdf",sep=""))
    plot(ATAC_all)
  dev.off()
  
  saveRDS(ATAC_all, paste(out_dir,"/Diffbind_consensus_peak_set_raw",sep=""))
  
  
  # Modify consensus peak set before saving it
  print("Processing consensus peak set")
  all_data_ATAC_all <- cbind(as.data.frame(ATAC_all$peaks[[1]]$Chr), 
                           as.data.frame(ATAC_all$peaks[[1]]$Start), 
                           as.data.frame(ATAC_all$peaks[[1]]$End))
  colnames(all_data_ATAC_all) <- c("Chr", "Start", "End")
  
  # extract read counts for each peak
  for (i in seq(1, length(ATAC_all$peaks))) {
    all_data_ATAC_all[, paste(ATAC_all$samples$SampleID[i],"raw",sep="_")] <- as.data.frame(ATAC_all$peaks[[i]]$Reads)
  }
  
  # name peaks
  rownames(all_data_ATAC_all) <-  make.names(paste("P.", all_data_ATAC_all$Chr, ".", all_data_ATAC_all$Start, sep=""), unique = TRUE)
  all_data_ATAC_all$Chr <- as.character(all_data_ATAC_all$Chr)
  
  # clean up the data set by excluding non-standard chromosomes and peaks at start/end of chromosomes
  # create a list of main chromosomes and their lengths
  chrom_sizes <- read.table(chr_sizes)
  colnames(chrom_sizes) <- c("Chr", "Length")
  chrom_sizes <- chrom_sizes[chrom_sizes$Chr %in% seq(1:chrom_number),]
  
  # only keep peaks located in main chromosomes
  print("Only keeping peaks in chromosomes:")
  print(chrom_sizes$Chr)
  all_data_ATAC_all <- all_data_ATAC_all[all_data_ATAC_all$Chr %in% chrom_sizes$Chr,]
  
  # exclude peaks extending beyond chromosome start
  all_data_ATAC_all <- all_data_ATAC_all[all_data_ATAC_all$Start >= 0,]
  
  # exclude peaks extending beyond chromosome end
  for (i in 1:nrow(chrom_sizes)) {
    tmp <- as.data.frame(cbind(rownames(all_data_ATAC_all)[all_data_ATAC_all$Chr == chrom_sizes$Chr[i]],
        all_data_ATAC_all[all_data_ATAC_all$Chr == chrom_sizes$Chr[i],"End"]))
    tmp$V2 <- as.integer(as.character(tmp$V2))
    if (i==1) {
        peaks_to_keep <- as.data.frame(tmp[tmp$V2 < chrom_sizes[i,"Length"],"V1"])
    }
    else {
        peaks_to_keep <- rbind(peaks_to_keep,as.data.frame(tmp[tmp$V2 < chrom_sizes[i,"Length"],"V1"]))
    }
  }    
  all_data_ATAC_all <- all_data_ATAC_all[rownames(all_data_ATAC_all) %in% peaks_to_keep[,1],]
  
  # manually generate RPKM values (normalise raw counts to peak library size)
  for (x in seq(4,ncol(all_data_ATAC_all))){
    all_data_ATAC_all[,ncol(all_data_ATAC_all)+1] <- (all_data_ATAC_all[,x]/(peak_width/1000))/(colSums(all_data_ATAC_all[x])/1e+06)
    colnames(all_data_ATAC_all)[ncol(all_data_ATAC_all)] <- paste(colnames(all_data_ATAC_all)[x],"rpkm",sep="_")
  }
  
  # save final output files
  print(paste("Saving final output file as ",
              paste(out_dir, "/Diffbind_consensus_peak_set", sep=""),sep=""))
  saveRDS(all_data_ATAC_all, paste(out_dir, "/Diffbind_consensus_peak_set", sep=""))
  
  return(all_data_ATAC_all)
}




## function to generate bed file for UCSC genome browser session and gimmescan analysis
# diffbind_file: Output of "diffbind" function
# out_dir: Directory containing output bed file
ucsc_gimmescan_bed <- function(diffbind_file, out_dir) {
  print(paste("Generating BED file for UCSC genome browser and gimmescan as ",
              paste(out_dir,"/Diffbind_consensus_peak_set_ucsc_gimmescan.bed",sep=""), sep=""))
  
  diffbind_counts <- readRDS(diffbind_file)
  peak_bed_file <- cbind(diffbind_counts[,1:3], rownames(diffbind_counts))
  colnames(peak_bed_file) <- c("Chr","Start","End","Name")
  rownames(peak_bed_file) <- NULL
  
  # transform chromosome nomenclature from eNsembl to UCSC
  xx <- peak_bed_file$Chr
  xx <- paste("chr", xx, sep="")
  xx <- gsub('chrKN', 'chrUn_KN', xx)
  xx <- gsub('\\.', 'v', xx)
  xx <- gsub('chrMT', 'chrM', xx)
  peak_bed_file$Chr <- xx
  
  colnames(peak_bed_file) <- NULL
  write.table(peak_bed_file, file=paste(out_dir,"/Diffbind_consensus_peak_set_ucsc_gimmescan.bed",sep=""), 
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  
  return(peak_bed_file)
}




## function to generate bed file for homer peak annotation
# diffbind_file: Output of "diffbind" function
# out_dir: Directory containing output bed file
homer_bed <- function(diffbind_file, out_dir) {
  print(paste("Generating BED file for Homer as ",
              paste(out_dir,"/Diffbind_consensus_peak_set_for_Homer.bed",sep=""), sep=""))
  
  diffbind_counts <- readRDS(diffbind_file)
  homer_bed_file <- cbind(as.data.frame(rownames(diffbind_counts)),diffbind_counts[,1:3],
                          as.data.frame(rep("+",dim(diffbind_counts)[1])))
  colnames(homer_bed_file) <- c("PeakID","Chr","Start","End","Strand")
  rownames(homer_bed_file) <- NULL
  write.table(homer_bed_file, file=paste(out_dir,"/Diffbind_consensus_peak_set_for_Homer.bed",sep=""), 
              row.names=FALSE, quote = FALSE, sep="\t")
  
  return(homer_bed_file)
}




## function to generate chromosome size file for UCSC genome browser
# chr_sizes: File with chromosome lengths (chromosome names in ensembl format)
# genome: Name of genome used in analysis
# out_dir: Directory containing output chromosome sizes file
ucsc_chrom_sizes <- function(chr_sizes, genome, out_dir) {
  print(paste("Generating chromosome sizes file for UCSC genome browser as ",
              paste(out_dir,"/",genome,"_ucsc.chrom.sizes",sep=""), sep=""))
  
  chrom_sizes <- read.table(chr_sizes)
  colnames(chrom_sizes) <- c("Chr", "Length")
  
  # transform chromosome nomenclature from eNsembl to UCSC
  xx <- chrom_sizes$Chr
  xx <- paste("chr", xx, sep="")
  xx <- gsub('chrKN', 'chrUn_KN', xx)
  xx <- gsub('\\.', 'v', xx)
  xx <- gsub('chrMT', 'chrM', xx)
  chrom_sizes$Chr <- xx
  
  colnames(chrom_sizes) <- NULL
  write.table(chrom_sizes, file=paste(out_dir,"/",genome,"_ucsc.chrom.sizes",sep=""), 
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  
  return(chrom_sizes)
}







### Analysis

print("Generating DiffBind consensus peak set")


## generate input for DiffBind
diffbind_peak_set <- diffbind(bam_dir=bam_dir, peak_dir=peak_dir, 
                              out_dir=out_dir, repeats_file=repeats_file,
                              peak_width=peak_width, chr_sizes=chr_sizes, chrom_number=chrom_number)


## generate bed file for UCSC genome browser session and gimmescan analysis
diffbind_file <- paste(out_dir, "/Diffbind_consensus_peak_set", sep="")
ucsc_gimmescan <- ucsc_gimmescan_bed(diffbind_file=diffbind_file, out_dir=out_dir)


## generate bed file of peak positions as input for gene annotation with homer
diffbind_file <- paste(out_dir, "/Diffbind_consensus_peak_set", sep="")
homer <- homer_bed(diffbind_file=diffbind_file, out_dir=out_dir)


## create chromosome sizes file with UCSC connotation
chr_sizes_ucsc <- ucsc_chrom_sizes(chr_sizes=chr_sizes, genome=genome, out_dir=genome_dir)





