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





### packages

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("JASPAR2020")) BiocManager::install("JASPAR2020")
if (!require("TFBSTools")) BiocManager::install("TFBSTools")
if (!require("GenomicFeatures")) BiocManager::install("GenomicFeatures")
if (!require("biomaRt")) BiocManager::install("biomaRt")




### functions

# subfunction to account for multiple gene ids annotated to the same gene name
gene_names_make_unique <- function(df, col_id, col_name) {
  tmp <- unique(df[,c(col_id, col_name)])
  tmp[,col_name] <- make.names(tmp[,col_name], unique=TRUE)
  df_1 <- merge(df, tmp, by=col_id)
  colnames(df_1)[ncol(df_1)] <- col_name 
  df_1[, paste(col_name, "x", sep=".")] <- NULL
  return(df_1)
}



# subfunction to map gene names from one species to another
map_gene_names <- function(genes, mart_1, mart_2, init_1, init_2, out_dir) {
  
  mapping <- getLDS(attributes = c("ensembl_gene_id","external_gene_name"),
                    filters = "ensembl_gene_id", values = rownames(genes), 
                    mart = mart_1,
                    attributesL = c("ensembl_gene_id","external_gene_name"), 
                    martL = mart_2)
  colnames(mapping) <- c(paste("gene_id_",init_1,sep=""), paste("gene_name_",init_1,sep=""),
                         paste("gene_id_",init_2,sep=""), paste("gene_name_",init_2,sep=""))
  
  mapping_1 <- gene_names_make_unique(mapping, 
                                      paste("gene_id_",init_1,sep=""), 
                                      paste("gene_name_",init_1,sep=""))
  mapping_1 <- mapping_1[,c(paste("gene_id_",init_1,sep=""), paste("gene_name_",init_1,sep=""),
                            paste("gene_id_",init_2,sep=""), paste("gene_name_",init_2,sep=""))]
  
  write.table(mapping_1, paste(out_dir, "/Gene_mapping_", init_1, "_to_", init_2, ".txt", sep=""), 
              sep="\t", quote=FALSE, row.names=FALSE)
  
  return(mapping_1)
} 



## function to generate files mapping human or mouse gene names to zebrafish or chicken gene names
map_gene_names_species <- function(gtf_file, organism, out_dir) {
  
  # extract gene data from gtf file
  txdb <- makeTxDbFromGFF(gtf_file, format="gtf", organism=organism)
  genes <- as.data.frame(genes(txdb))
  
  # create a mapping of ensembl gene id and gene name across species
  initials_2 <- "hs"
  initials_3 <- "mm"
  ensembl_hs <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                        host = "https://dec2021.archive.ensembl.org/")
  ensembl_mm <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl",
                        host = "https://dec2021.archive.ensembl.org/")
  
  if (organism == "Danio rerio") {
    initials_1 <- "dr"
    ensembl_dr <- useMart(biomart = "ensembl", dataset = "drerio_gene_ensembl",
                          host = "https://dec2021.archive.ensembl.org/")
    dr_to_hs <- map_gene_names(genes=genes, mart_1=ensembl_dr, mart_2=ensembl_hs, 
                               init_1=initials_1, init_2=initials_2, out_dir=out_dir)
    dr_to_mm <- map_gene_names(genes=genes, mart_1=ensembl_dr, mart_2=ensembl_mm, 
                               init_1=initials_1, init_2=initials_3, out_dir=out_dir)
  
  } else if (organism == "Gallus gallus") {
    initials_1 <- "gg"
    ensembl_gg <- useMart(biomart = "ensembl", dataset = "ggallus_gene_ensembl",
                          host = "https://dec2021.archive.ensembl.org/")
    gg_to_hs <- map_gene_names(genes=genes, mart_1=ensembl_gg, mart_2=ensembl_hs, 
                               init_1=initials_1, init_2=initials_2, out_dir=out_dir)
    gg_to_mm <- map_gene_names(genes=genes, mart_1=ensembl_gg, mart_2=ensembl_mm, 
                               init_1=initials_1, init_2=initials_3, out_dir=out_dir)
  
  } else if (organism == "Homo sapiens") {
    hs_to_mm <- map_gene_names(genes=genes, mart_1=ensembl_hs, mart_2=ensembl_mm, 
                               init_1=initials_2, init_2=initials_3, out_dir=out_dir)
  
  } else if (organism == "Mus musculus") {
    mm_to_hs <- map_gene_names(genes=genes, mart_1=ensembl_mm, mart_2=ensembl_hs, 
                               init_1=initials_3, init_2=initials_2, out_dir=out_dir)
  }
  
}



# subfunction to convert human/mouse TF names in motif2factors file to chicken or zebrafish TF names
motif2factors_to_gg_dr <- function(motif2factors, target_organism_initials,
                                   target_to_hs_gene_mapping, target_to_mm_gene_mapping) {
  
  init <- target_organism_initials
  
  # loop over rows of motif2factors file
  for (i in seq(1,nrow(motif2factors))) {
    tmp <- list()
  
    # grep chicken/zebrafish factor name based on human factor name
    if (motif2factors[i,"Factor"] %in% target_to_hs_gene_mapping$gene_name_hs) {
      tmp <- target_to_hs_gene_mapping[grep(paste("_",motif2factors[i,"Factor"],"_",sep=""),
                                            paste("_",target_to_hs_gene_mapping$gene_name_hs,"_",sep="")), 
                                       paste("gene_name_",init,sep="")]
  
    # or grep chicken/zebrafish factor name based on mouse factor name
    } else if (!motif2factors[i,"Factor"] %in% target_to_hs_gene_mapping$gene_name_hs &
               motif2factors[i,"Factor"] %in% target_to_mm_gene_mapping$gene_name_mm) {
      tmp <- target_to_mm_gene_mapping[grep(paste("_",motif2factors[i,"Factor"],"_",sep=""),
                                            paste("_",target_to_mm_gene_mapping$gene_name_mm,"_",sep="")), 
                                       paste("gene_name_",init,sep="")]
    }
  
    # split chicken or zebrafish duplicated factor names into separate row entries
    if (length(tmp) > 0) {
      for (k in seq(1,length(tmp))) {
        if (i==1 & k==1) {
          motif2factors_new <- data.frame(Motif=motif2factors[i,"Motif"], Factor=tmp[k],
                                         Evidence="JASPAR", Curated="Y")
        } else {
          motif2factors_new[nrow(motif2factors_new)+1,"Motif"] <- motif2factors[i,"Motif"]
          motif2factors_new[nrow(motif2factors_new),"Factor"] <- tmp[k]
          motif2factors_new[nrow(motif2factors_new),"Evidence"] <- "JASPAR"
          motif2factors_new[nrow(motif2factors_new),"Curated"] <- "Y"
        }
      }
    }  
  }
  
  return(motif2factors_new)
}



# subfunction to convert human/mouse TF names in motif2factors file to all human or mouse TF names
motif2factors_to_hs_mm <- function(motif2factors, target_organism_initials, mm_hs_gene_mapping) {
  
  init <- target_organism_initials
  
  # loop over rows of motif2factors file
  for (i in seq(1,nrow(motif2factors))) {
    tmp <- list()
  
    # grep factor name based on human factor name
    if (motif2factors[i,"Factor"] %in% mm_hs_gene_mapping$gene_name_hs) {
      tmp <- mm_hs_gene_mapping[grep(paste("_",motif2factors[i,"Factor"],"_",sep=""),
                                     paste("_",mm_hs_gene_mapping$gene_name_hs,"_",sep="")), 
                                paste("gene_name_",init,sep="")]
  
    # or grep factor name based on mouse factor name
    } else if (!motif2factors[i,"Factor"] %in% mm_hs_gene_mapping$gene_name_hs &
               motif2factors[i,"Factor"] %in% mm_hs_gene_mapping$gene_name_mm) {
      tmp <- mm_hs_gene_mapping[grep(paste("_",motif2factors[i,"Factor"],"_",sep=""),
                                     paste("_",mm_hs_gene_mapping$gene_name_mm,"_",sep="")), 
                                paste("gene_name_",init,sep="")]
    }
  
    # split duplicated factor names into separate row entries
    if (length(tmp) > 0) {
      for (k in seq(1,length(tmp))) {
        if (i==1 & k==1) {
          motif2factors_new <- data.frame(Motif=motif2factors[i,"Motif"], Factor=tmp[k],
                                         Evidence="JASPAR", Curated="Y")
        } else {
          motif2factors_new[nrow(motif2factors_new)+1,"Motif"] <- motif2factors[i,"Motif"]
          motif2factors_new[nrow(motif2factors_new),"Factor"] <- tmp[k]
          motif2factors_new[nrow(motif2factors_new),"Evidence"] <- "JASPAR"
          motif2factors_new[nrow(motif2factors_new),"Curated"] <- "Y"
        }
      }
    }  
  }
  
  return(motif2factors_new)
}



## function to annotate binding motifs to transcription factors using JASPAR2020
map_tf_motifs_names_jaspar <- function(organism, out_dir) {
  
  # generate motif2factors file
  opts <- list()
  opts[["collection"]] <- "CORE"
  opts[["tax_group"]] <- "vertebrates"
  jaspar_motifs <- getMatrixSet(JASPAR2020, opts)
  
  tf_tfbs_names <- data.frame(Motif=NA, Factor=NA)
  k <- 1
  for (i in seq(1,length(jaspar_motifs))) {
    # split TF::TF names onto separate rows
    if (grepl("::", name(jaspar_motifs[[i]]), fixed = TRUE)) {
      tmp <- unlist(strsplit(name(jaspar_motifs[[i]]), "::"))
      for (l in seq(1, length(tmp))) {
        tf_tfbs_names[k,"Motif"] <- ID(jaspar_motifs[[i]])
        tf_tfbs_names[k,"Factor"] <- tmp[l]
        k <- k+1
      }
    } else {
      tf_tfbs_names[k,"Motif"] <- ID(jaspar_motifs[[i]])
      tf_tfbs_names[k,"Factor"] <- name(jaspar_motifs[[i]])
      k <- k+1
    }
  }
  
  # delete "(var.x)" substrings in Factor names
  tf_tfbs_names$Factor <- gsub("\\(var.*\\)", "", tf_tfbs_names$Factor)
  
  # split "EWSR1-FLI1" name and replace
  tmp <- tf_tfbs_names[grep("EWSR1-FLI1", tf_tfbs_names$Factor),]
  tmp_1 <- unlist(strsplit(tmp$Factor, "-"))
  
  for (l in seq(1, length(tmp_1))){
    tf_tfbs_names[nrow(tf_tfbs_names)+1,"Motif"] <- tmp$Motif
    tf_tfbs_names[nrow(tf_tfbs_names),"Factor"] <- tmp_1[l]
  }
  
  tf_tfbs_names <- tf_tfbs_names[!grepl("EWSR1-FLI1", tf_tfbs_names$Factor),]
  
  
  # replace human/mouse TF names in motif database with TF names of species analysed
  initials_2 <- "hs"
  initials_3 <- "mm"
  
  if (organism == "Danio rerio") {
    initials_1 <- "dr"
    dr_to_hs <- read.table(paste(out_dir, "/Gene_mapping_", initials_1, "_to_", initials_2, ".txt", sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    dr_to_mm <- read.table(paste(out_dir, "/Gene_mapping_", initials_1, "_to_", initials_3, ".txt", sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
    motif2factors_new <- motif2factors_to_gg_dr(motif2factors=tf_tfbs_names, 
                                                target_organism_initials=initials_1,
                                                target_to_hs_gene_mapping=dr_to_hs, 
                                                target_to_mm_gene_mapping=dr_to_mm)

  # subfunction to manually add motifs to motif2factors file
  add_factor <- function(df=motif2factors_new, motif, factor, evidence="JASPAR", curated="Y") {
    df[nrow(df)+1,"Motif"] <- motif
    df[nrow(df),"Factor"] <- factor
    df[nrow(df),"Evidence"] <- evidence
    df[nrow(df),"Curated"] <- curated
    return(df)
  }

  # add factors that were not automatically added
	motif2factors_new <- add_factor(motif="MA0634.1", factor="alx3")
	motif2factors_new <- add_factor(motif="MA1463.1", factor="argfx")
	motif2factors_new <- add_factor(motif="MA0602.1", factor="arid5a")
	motif2factors_new <- add_factor(motif="MA0604.1", factor="atf1")
	motif2factors_new <- add_factor(motif="MA0876.1", factor="bsx")
	motif2factors_new <- add_factor(motif="MA0878.2", factor="cdx1a")
	motif2factors_new <- add_factor(motif="MA0878.2", factor="cdx1b")
	motif2factors_new <- add_factor(motif="MA0465.2", factor="cdx2")
	motif2factors_new <- add_factor(motif="MA1473.1", factor="cdx4")
	motif2factors_new <- add_factor(motif="MA0838.1", factor="cebpg")
	motif2factors_new <- add_factor(motif="MA1636.1", factor="cebpg")
	motif2factors_new <- add_factor(motif="MA0637.1", factor="cenpb")
	motif2factors_new <- add_factor(motif="MA1474.1", factor="creb3l4")
	motif2factors_new <- add_factor(motif="MA1475.1", factor="creb3l4")
	motif2factors_new <- add_factor(motif="MA1479.1", factor="dmrtc2")
	motif2factors_new <- add_factor(motif="MA1480.1", factor="dprx")
	motif2factors_new <- add_factor(motif="MA0611.1", factor="dux")
	motif2factors_new <- add_factor(motif="MA0884.1", factor="duxa")
	motif2factors_new <- add_factor(motif="MA0864.2", factor="e2f2")
	motif2factors_new <- add_factor(motif="MA0162.4", factor="egr1")
	motif2factors_new <- add_factor(motif="MA0473.3", factor="elf1")
	motif2factors_new <- add_factor(motif="MA0641.1", factor="elf4")
	motif2factors_new <- add_factor(motif="MA0136.2", factor="elf5")
	motif2factors_new <- add_factor(motif="MA0644.1", factor="esx1")
	motif2factors_new <- add_factor(motif="MA1484.1", factor="ets2")
	motif2factors_new <- add_factor(motif="MA0762.1", factor="etv2")
	motif2factors_new <- add_factor(motif="MA0763.1", factor="etv3")
	motif2factors_new <- add_factor(motif="MA0765.2", factor="etv5a")
	motif2factors_new <- add_factor(motif="MA0765.2", factor="etv5b")
	motif2factors_new <- add_factor(motif="MA0645.1", factor="etv6")
	motif2factors_new <- add_factor(motif="MA1683.1", factor="foxa3")
	motif2factors_new <- add_factor(motif="MA0846.1", factor="foxc2")
	motif2factors_new <- add_factor(motif="MA0847.2", factor="foxd2")
	motif2factors_new <- add_factor(motif="MA0614.1", factor="foxj2")
	motif2factors_new <- add_factor(motif="MA0851.1", factor="foxj3")
	motif2factors_new <- add_factor(motif="MA0157.2", factor="foxo3a")
	motif2factors_new <- add_factor(motif="MA0157.2", factor="foxo3b")
	motif2factors_new <- add_factor(motif="MA0481.3", factor="foxp1a")
	motif2factors_new <- add_factor(motif="MA0481.3", factor="foxp1b")
	motif2factors_new <- add_factor(motif="MA0850.1", factor="foxp3a")
	motif2factors_new <- add_factor(motif="MA0850.1", factor="foxp3b")
	motif2factors_new <- add_factor(motif="MA0035.4", factor="gata1a")
	motif2factors_new <- add_factor(motif="MA0035.4", factor="gata1b")
	motif2factors_new <- add_factor(motif="MA0646.1", factor="gcm1")
	motif2factors_new <- add_factor(motif="MA0891.1", factor="gsc2")
	motif2factors_new <- add_factor(motif="MA0894.1", factor="hesx1")
	motif2factors_new <- add_factor(motif="MA0899.1", factor="hoxa10a")
	motif2factors_new <- add_factor(motif="MA0899.1", factor="hoxa10b")
	motif2factors_new <- add_factor(motif="MA0650.2", factor="hoxa13a")
	motif2factors_new <- add_factor(motif="MA0650.2", factor="hoxa13b")
	motif2factors_new <- add_factor(motif="MA0900.2", factor="hoxa2a")
	motif2factors_new <- add_factor(motif="MA0900.2", factor="hoxa2b")
	motif2factors_new <- add_factor(motif="MA1496.1", factor="hoxa4a")
	motif2factors_new <- add_factor(motif="MA1496.1", factor="hoxa4b")
	motif2factors_new <- add_factor(motif="MA1497.1", factor="hoxa6")
	motif2factors_new <- add_factor(motif="MA1498.1", factor="hoxb6b")
	motif2factors_new <- add_factor(motif="MA0909.2", factor="hoxd13a")
	motif2factors_new <- add_factor(motif="MA0909.2", factor="hoxd13b")
	motif2factors_new <- add_factor(motif="MA0910.2", factor="hoxd8")
	motif2factors_new <- add_factor(motif="MA0486.2", factor="hsf1")
	motif2factors_new <- add_factor(motif="MA0654.1", factor="isx")
	motif2factors_new <- add_factor(motif="MA0493.1", factor="klf1")
	motif2factors_new <- add_factor(motif="MA1511.1", factor="klf10")
	motif2factors_new <- add_factor(motif="MA0657.1", factor="klf13")
	motif2factors_new <- add_factor(motif="MA1515.1", factor="klf2a")
	motif2factors_new <- add_factor(motif="MA1515.1", factor="klf2b")
	motif2factors_new <- add_factor(motif="MA1518.1", factor="lhx1a")
	motif2factors_new <- add_factor(motif="MA1518.1", factor="lhx1b")
	motif2factors_new <- add_factor(motif="MA1519.1", factor="lhx5")
	motif2factors_new <- add_factor(motif="MA0058.3", factor="max")
	motif2factors_new <- add_factor(motif="MA0621.1", factor="mixa")
	motif2factors_new <- add_factor(motif="MA1523.1", factor="msantd3")
	motif2factors_new <- add_factor(motif="MA0665.1", factor="msc")
	motif2factors_new <- add_factor(motif="MA0666.1", factor="msx1a")
	motif2factors_new <- add_factor(motif="MA0666.1", factor="msx1b")
	motif2factors_new <- add_factor(motif="MA0708.1", factor="msx2a")
	motif2factors_new <- add_factor(motif="MA0708.1", factor="msx2b")
	motif2factors_new <- add_factor(motif="MA0056.2", factor="mzf1")
	motif2factors_new <- add_factor(motif="MA0057.1", factor="mzf1")
	motif2factors_new <- add_factor(motif="MA0623.2", factor="neurog1")
	motif2factors_new <- add_factor(motif="MA0669.1", factor="neurog2")
	motif2factors_new <- add_factor(motif="MA1642.1", factor="neurog2")
	motif2factors_new <- add_factor(motif="MA1643.1", factor="nfib")
	motif2factors_new <- add_factor(motif="MA0778.1", factor="nfkb2")
	motif2factors_new <- add_factor(motif="MA0048.2", factor="nhlh1")
	motif2factors_new <- add_factor(motif="MA1529.1", factor="nhlh2")
	motif2factors_new <- add_factor(motif="MA0673.1", factor="nkx2.8")
	motif2factors_new <- add_factor(motif="MA0125.1", factor="nobox")
	motif2factors_new <- add_factor(motif="MA0626.1", factor="npas2")
	motif2factors_new <- add_factor(motif="MA0160.1", factor="nr4a2a")
	motif2factors_new <- add_factor(motif="MA0160.1", factor="nr4a2b")
	motif2factors_new <- add_factor(motif="MA1544.1", factor="ovol1a")
	motif2factors_new <- add_factor(motif="MA1544.1", factor="ovol1b")
	motif2factors_new <- add_factor(motif="MA0068.2", factor="pax4")
	motif2factors_new <- add_factor(motif="MA0070.1", factor="pbx1a")
	motif2factors_new <- add_factor(motif="MA0070.1", factor="pbx1b")
	motif2factors_new <- add_factor(motif="MA1113.2", factor="pbx2")
	motif2factors_new <- add_factor(motif="MA0132.2", factor="pdx1")
	motif2factors_new <- add_factor(motif="MA1547.1", factor="pitx2")
	motif2factors_new <- add_factor(motif="MA0714.1", factor="pitx3")
	motif2factors_new <- add_factor(motif="MA0789.1", factor="pou3f4")
	motif2factors_new <- add_factor(motif="MA0066.1", factor="pparg")
	motif2factors_new <- add_factor(motif="MA0715.1", factor="prop1")
	motif2factors_new <- add_factor(motif="MA0075.3", factor="prrx2")
	motif2factors_new <- add_factor(motif="MA0729.1", factor="raraa")
	motif2factors_new <- add_factor(motif="MA0729.1", factor="rarab")
	motif2factors_new <- add_factor(motif="MA0730.1", factor="raraa")
	motif2factors_new <- add_factor(motif="MA0730.1", factor="rarab")
	motif2factors_new <- add_factor(motif="MA0857.1", factor="rarb")
	motif2factors_new <- add_factor(motif="MA0858.1", factor="rarb")
	motif2factors_new <- add_factor(motif="MA1552.1", factor="rarb")
	motif2factors_new <- add_factor(motif="MA0718.1", factor="rx1")
	motif2factors_new <- add_factor(motif="MA0717.1", factor="rx3")
	motif2factors_new <- add_factor(motif="MA0629.1", factor="rhox11")
	motif2factors_new <- add_factor(motif="MA0719.1", factor="rhoxf1")
	motif2factors_new <- add_factor(motif="MA0743.2", factor="scrt1a")
	motif2factors_new <- add_factor(motif="MA0743.2", factor="scrt1b")
	motif2factors_new <- add_factor(motif="MA0745.2", factor="snai2")
	motif2factors_new <- add_factor(motif="MA1559.1", factor="snai3")
	motif2factors_new <- add_factor(motif="MA1560.1", factor="sohlh2")
	motif2factors_new <- add_factor(motif="MA1563.1", factor="sox18")
	motif2factors_new <- add_factor(motif="MA0111.1", factor="spz1")
	motif2factors_new <- add_factor(motif="MA0805.1", factor="tbx1")
	motif2factors_new <- add_factor(motif="MA0632.2", factor="tcfl5")
	motif2factors_new <- add_factor(motif="MA1121.1", factor="tead2")
	motif2factors_new <- add_factor(motif="MA0809.2", factor="tead4")
	motif2factors_new <- add_factor(motif="MA0797.1", factor="tgif2")
	motif2factors_new <- add_factor(motif="MA0597.1", factor="thap1")
	motif2factors_new <- add_factor(motif="MA1577.1", factor="tlx2")
	motif2factors_new <- add_factor(motif="MA1123.2", factor="twist1a")
	motif2factors_new <- add_factor(motif="MA1123.2", factor="twist1b")
	motif2factors_new <- add_factor(motif="MA0748.2", factor="yy2")
	motif2factors_new <- add_factor(motif="MA0749.1", factor="zbed1")
	motif2factors_new <- add_factor(motif="MA1581.1", factor="zbtb6")
	motif2factors_new <- add_factor(motif="MA0750.2", factor="zbtb7a")
	motif2factors_new <- add_factor(motif="MA1651.1", factor="zfp42")
	motif2factors_new <- add_factor(motif="MA1583.1", factor="zfp57")
	motif2factors_new <- add_factor(motif="MA0146.2", factor="zfx")
	motif2factors_new <- add_factor(motif="MA1585.1", factor="zkscan1")
	motif2factors_new <- add_factor(motif="MA1652.1", factor="zkscan5")
	motif2factors_new <- add_factor(motif="MA1587.1", factor="znf135")
	motif2factors_new <- add_factor(motif="MA1588.1", factor="znf136")
	motif2factors_new <- add_factor(motif="MA1589.1", factor="znf140")
	motif2factors_new <- add_factor(motif="MA1654.1", factor="znf16")
	motif2factors_new <- add_factor(motif="MA1124.1", factor="znf24")
	motif2factors_new <- add_factor(motif="MA0528.2", factor="znf236")
	motif2factors_new <- add_factor(motif="MA1630.1", factor="znf281a")
	motif2factors_new <- add_factor(motif="MA1630.1", factor="znf281b")
	motif2factors_new <- add_factor(motif="MA1593.1", factor="znf317")
	motif2factors_new <- add_factor(motif="MA0130.1", factor="znf354c")
	motif2factors_new <- add_factor(motif="MA1594.1", factor="znf382")
	motif2factors_new <- add_factor(motif="MA0116.1", factor="znf423")
	motif2factors_new <- add_factor(motif="MA1656.1", factor="znf449")
	motif2factors_new <- add_factor(motif="MA1596.1", factor="znf460")
	motif2factors_new <- add_factor(motif="MA1597.1", factor="znf528")
	motif2factors_new <- add_factor(motif="MA1599.1", factor="znf682")
	motif2factors_new <- add_factor(motif="MA1600.1", factor="znf684")
	motif2factors_new <- add_factor(motif="MA1601.1", factor="znf75d")
	motif2factors_new <- add_factor(motif="MA1155.1", factor="zscan4")
	motif2factors_new <- add_factor(motif="MA0083.3", factor="srfa")
	motif2factors_new <- add_factor(motif="MA0083.3", factor="srfb")
	motif2factors_new <- add_factor(motif="MA0511.2", factor="runx2a")
	motif2factors_new <- add_factor(motif="MA0511.2", factor="runx2b")
	motif2factors_new <- add_factor(motif="MA0031.1", factor="foxd1")
  
  } else if (organism == "Gallus gallus") {
    initials_1 <- "gg"
    gg_to_hs <- read.table(paste(out_dir, "/Gene_mapping_", initials_1, "_to_", initials_2, ".txt", sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    gg_to_mm <- read.table(paste(out_dir, "/Gene_mapping_", initials_1, "_to_", initials_3, ".txt", sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    motif2factors_new <- motif2factors_to_gg_dr(motif2factors=tf_tfbs_names, 
                                                target_organism_initials=initials_1,
                                                target_to_hs_gene_mapping=gg_to_hs, 
                                                target_to_mm_gene_mapping=gg_to_mm)
  
  } else if (organism == "Homo sapiens") {
    hs_to_mm <- read.table(paste(out_dir, "/Gene_mapping_", initials_2, "_to_", initials_3, ".txt", sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    motif2factors_new <- motif2factors_to_hs_mm(motif2factors=tf_tfbs_names, 
                                                target_organism_initials=initials_2, 
                                                mm_hs_gene_mapping=hs_to_mm)
  
  } else if (organism == "Mus musculus") {
    mm_to_hs <- read.table(paste(out_dir, "/Gene_mapping_", initials_3, "_to_", initials_2, ".txt", sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    motif2factors_new <- motif2factors_to_hs_mm(motif2factors=tf_tfbs_names, 
                                                target_organism_initials=initials_3, 
                                                mm_hs_gene_mapping=mm_to_hs)
  }
  
  # filter out duplicate rows
  motif2factors_new <- motif2factors_new %>% distinct(Motif, Factor, Evidence, Curated, .keep_all = TRUE)
  
  # check if all motifs are retained in new file
  #unique(tf_tfbs_names$Motif) %in% unique(motif2factors_dr$Motif)
  
  # Ananse will fail if there are underscores in the gene names
  #grep("_", motif2factors_dr$Factor)
  
  return(motif2factors_new)
}



### Analysis

# define organism
if (species == "human") {
  organism <- "Homo sapiens"
} else if (species == "mouse") {
  organism <- "Mus musculus"
} else if (species == "chicken") {
  organism <- "Gallus gallus"
} else if (species == "zebrafish") {
  organism <- "Danio rerio"
}

# create gene name mappings between species
map_gene_names_species(gtf_file=gtf, organism=organism, out_dir=out_dir)

# adjust TF names to organism
motif2factors <- map_tf_motifs_names_jaspar(organism=organism, out_dir=out_dir)

write.table(motif2factors, paste(out_dir,"/JASPAR2020.motif2factors.txt",sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)



