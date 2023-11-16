
### Run this script to analyze the output of Ananse Network




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
if (!require("DESeq2")) BiocManager::install("DESeq2")
if (!require("edgeR")) install.packages("edgeR")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("igraph")) install.packages("igraph")




### functions

## function to extract and plot statistics from entire network
ananse_net_stats <- function(input_df, out_prefix, factors, out_dir) {

  input_df[,1] <- tolower(input_df[,1])

  # split TF target column
  input_df[,c("source","target")] <- read.table(text=input_df$tf_target, 
                                                sep="_", strip.white = TRUE)

  # network file: generate df with number of entries
  inner_conn <- length(input_df$target[input_df$target %in% input_df$source])
  net_data_1 <- as.data.frame(c(inner_conn, (nrow(input_df) - inner_conn)))
  colnames(net_data_1) <- "number_connections"
  net_data_1$type <- c("inner", "peripheral")
  net_data_1$condition <- out_prefix

  print(paste("Total number of connections:", nrow(input_df)))
  print(paste("Number of inner connections:", inner_conn))
  print(paste("Number of peripheral connections:", (nrow(input_df) - inner_conn)))

  # plot bar chart
  gg <- ggplot(net_data_1, aes(x=condition, y=number_connections)) +
  geom_bar(aes(fill=type),stat="identity", width=0.8) +
  scale_fill_manual(breaks=c("inner","peripheral"),
    name="Type of\nconnection",
    labels=c("Inner connections","Peripheral connections"),
    values=c("darkorange", "darkgreen")) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 18, colour = "black", face = "plain", angle=30, hjust=1),
    axis.text.y = element_text(size = 18, colour = "black", face = "plain"),
    axis.title.y = element_text(size = 18, face = "plain", vjust = 3),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "right",
    legend.title = element_text(size = 18, colour = "black", face = "plain"),
    legend.text = element_text(size = 18, colour = "black"),
    plot.margin = margin(0, 0, 0, 1, "cm"),
    strip.text.x = element_text(size=16, face="plain"),
    strip.background = element_blank()) +
  scale_y_continuous(name="Number of connections in network") +
  scale_x_discrete(name="") +
  labs(title="Network connections")

  ggsave(filename=paste(out_dir, "/Ananse_network_", out_prefix, "_connections.pdf", sep=""), 
         plot=gg, device="pdf", dpi=300, useDingbats=FALSE)


  # network file: generate df with number of TFs, targets
  number_tf_targets <- length(unique(input_df$target[input_df$target %in% factors]))
  net_data_2 <- as.data.frame(c(length(unique(input_df$source)),
                                number_tf_targets, 
                                (length(unique(input_df$target))-number_tf_targets)))
  colnames(net_data_2) <- "number_factors"
  net_data_2$type <- c("Source", "Target", "Target")
  net_data_2$target <- c("1_TF", "1_TF", "2_non-TF")
  net_data_2$condition <- out_prefix

  print(paste("Number of TFs:", length(unique(input_df$source))))
  print(paste("Number of targets that are TFs:", number_tf_targets))
  print(paste("Number of targets that are non-TFs:", (length(unique(input_df$target))-number_tf_targets)))

  # plot bar chart
  gg <- ggplot(net_data_2, aes(x=type, y=number_factors)) +
  geom_bar(aes(fill=target),stat="identity", width=0.8) +
  scale_fill_manual(breaks=c("1_TF","2_non-TF"),
    name="Type of\nfactor",
    labels=c("TF","non-TF"),
    values=c("darkorange", "darkgreen")) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 18, colour = "black", face = "plain", angle=30, hjust=1),
    axis.text.y = element_text(size = 18, colour = "black", face = "plain"),
    axis.title.y = element_text(size = 18, face = "plain", vjust = 3),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "right",
    legend.title = element_text(size = 18, colour = "black", face = "plain"),
    legend.text = element_text(size = 18, colour = "black"),
    plot.margin = margin(0, 0, 0, 1, "cm"),
    strip.text.x = element_text(size=16, face="plain"),
    strip.background = element_blank()) +
  scale_y_continuous(name="Number of factors in network") +
  scale_x_discrete(name="") +
  labs(title="Network factors")

  ggsave(filename=paste(out_dir, "/Ananse_network_", out_prefix, "_factors.pdf", sep=""), 
         plot=gg, device="pdf", dpi=300, useDingbats=FALSE)


  # network file: generate df with cumulative binding score per TF, number of targets per TF
  net_data_3 <- input_df %>% group_by(source) %>% summarise(sum_binding_probs=sum(prob))
  name_col_2 <- paste("sum_binding_probs_", out_prefix, sep="")
  net_data_3 <- net_data_3[order(net_data_3$sum_binding_probs, decreasing=TRUE),]
  net_data_3$source <- factor(net_data_3$source, levels = net_data_3$source)

  # plot bar chart
  gg <- ggplot(net_data_3, aes(x=source, y=sum_binding_probs)) +
  geom_bar(aes(),stat="identity", width=1) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(size = 1, colour = "black", face = "plain", angle=90, hjust=1),
    axis.text.y = element_text(size = 18, colour = "black", face = "plain"),
    axis.title.y = element_text(size = 18, face = "plain", vjust = 3),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "none",
    legend.title = element_text(size = 18, colour = "black", face = "plain"),
    legend.text = element_text(size = 18, colour = "black"),
    plot.margin = margin(0, 0.5, 0, 1, "cm"),
    strip.text.x = element_text(size=16, face="plain"),
    strip.background = element_blank()) +
  scale_y_continuous(name="Cumulative binding probability") +
  scale_x_discrete(name="") +
  labs(title="Network factor binding")

  ggsave(filename=paste(out_dir, "/Ananse_network_", out_prefix, 
                        "_factor_binding.pdf", sep=""), 
         plot=gg, device="pdf", dpi=300, useDingbats=FALSE)

  return(list(net_data_1, net_data_2, net_data_3))
}



## function to subset network
subset_network <- function(network, source_genes=NULL, target_genes=NULL, prob_cutoff=0, edge_number=0,
                           out_prefix, out_dir) {

  # split tf_target column
  network$tf_target <- tolower(network$tf_target)
  network[,c("source","target")] <- read.table(text=network$tf_target, sep="_", strip.white = TRUE)

  # create dataframe indicating which parameters have been set
  tmp_df <- data.frame(parameter=c("source","target","prob","edges"), dummy=0)

  # extract egdes to keep separately for each parameter and adjust tmp_df
  if (length(source_genes) > 0) {
    source <- network[network$source %in% source_genes,"tf_target"]
    tmp_df[tmp_df$parameter == "source", "dummy"] <- 1
  }
    
  if (length(target_genes) > 0) {
    target <- network[network$target %in% target_genes,"tf_target"]
    tmp_df[tmp_df$parameter == "target", "dummy"] <- 1
  }

  if (prob_cutoff!=0) {
    prob <- network[network$prob > prob_cutoff,"tf_target"]
    tmp_df[tmp_df$parameter == "prob", "dummy"] <- 1
  }

  if (edge_number!=0) {
    network <- network[order(network$prob, decreasing=TRUE),]
    edges <- network[1:edge_number,"tf_target"]
    tmp_df[tmp_df$parameter == "edges", "dummy"] <- 1
  }

  # find intersection of edges if multiple parameters are set
  params <- tmp_df[tmp_df$dummy == 1, "parameter"]
  if (length(params) > 1) {  
    for (i in seq(1, length(params))) {
      tmp <- eval(as.name(params[i]))
      if (i==1) {
        subnet_edges <- tmp
      } else { 
        subnet_edges <- intersect(subnet_edges, tmp)
      }
    }
  } else if (length(params) == 1) {
    subnet_edges <- eval(as.name(params))
  }
  
  # create sub network
  subnet <- network[network$tf_target %in% subnet_edges,]
  subnet <- subnet[order(subnet$prob, decreasing=TRUE),]

  # check that subnet is correct
  for (i in params) {
    if (i=="source") {
      subnet <- subnet[subnet$source %in% source_genes,]
    } else if (i=="target") {
      subnet <- subnet[subnet$target %in% target_genes,]
    } else if (i=="prob") {
      subnet <- subnet[subnet$prob > prob_cutoff,]
    } else if (i=="edges") {
      subnet <- subnet[1:edge_number,]
    }
  }

  write.table(subnet, paste(out_dir, "/Ananse_network_", out_prefix, ".txt", sep=""), sep="\t",
              row.names=FALSE, col.names=TRUE, quote=FALSE)
  return(subnet)
}



## function to plot subnetwork as graph and compute centrality
net_graph <- function(subnet, out_prefix, motif_db="jaspar2020", out_dir) {

  # create graph
  g <- graph.data.frame(subnet[,c("source", "target")], directed=TRUE, vertices=NULL)
  E(g)$weight <- subnet[,"prob"]

  #pdf(file = paste("Ananse_network_graph_", motifs, "_", out_prefix, ".pdf", sep=""))
  #  plot.igraph(g,vertex.label=V(g)$name,layout=layout.fruchterman.reingold, 
  #    vertex.label.color="black",edge.color="black",
  #    edge.width=E(net)$weight, edge.arrow.size=1)
  #dev.off()

  # centrality
  res_centr <- as.data.frame(betweenness(g))
  colnames(res_centr) <- "betweenness_centrality"
  res_centr$eigenvector_centrality <- evcent(g)$vector

  res_centr$degree_centrality_all <- degree(g, mode = "all") / (vcount(g) - 1)
  res_centr$degree_centrality_in <- degree(g, mode = "in") / (vcount(g) - 1)
  res_centr$degree_centrality_out <- degree(g, mode = "out") / (vcount(g) - 1)

  #res_centr$clustering_coefficient_weighted <- transitivity(g,vids=V(g),type="weighted",isolates="zero")

  data <- walktrap.community(g, modularity=TRUE)
  res_centr$community_random_walk <- data$membership

  res_centr <- res_centr[order(res_centr$eigenvector_centrality, decreasing=TRUE),]
  
  write.csv(res_centr, paste(out_dir, "/Ananse_network_centrality_", motif_db, "_", 
                             out_prefix, ".csv", sep=""), 
            row.names=TRUE, quote=FALSE)
  return(res_centr)
}




## function to compute the relative column sums of two dataframes with matching column names
# and plot the column sums as dot plot
plot_compare_colsums <- function(df1, df2, out_prefix1, out_prefix2, genes=10, 
                                 motifs="jaspar2020", out_dir) {
  sum1 <- as.data.frame(colSums(df1))
  colnames(sum1) <- "colsum1"
  sum1$name <- rownames(sum1)
  sum2 <- as.data.frame(colSums(df2))
  colnames(sum2) <- "colsum2"
  sum2$name <- rownames(sum2)

  sum_all <- merge(sum1, sum2, by="name", all=TRUE)
  sum_all[is.na(sum_all)] <- 0

  # normalise colsums to number of input rows
  sum_all$colsum1 <- sum_all$colsum1 / nrow(df1)
  sum_all$colsum2 <- sum_all$colsum2 / nrow(df2)

  # min max scaling
  min_max_scaling <- function(input) {
    scaled <- (input - min(input)) /(max(input)-min(input))
    return(scaled)
  }
  sum_all$colsum1 <- min_max_scaling(sum_all$colsum1)
  sum_all$colsum2 <- min_max_scaling(sum_all$colsum2)

  # prepare line to add to plot
  slope_abline <- 1
  icept_abline <- 0

  # prepare genes to be highlighted  
  if (is.character(genes)) {
    genes_df <- data.frame(name=NA, colsum1=NA, colsum2=NA)
    for (i in seq(1,length(genes))) {
      genes_df[i,] <- sum_all[grep(paste("_",genes[i],"_",sep=""), 
                                            paste("_",sum_all$name,"_",sep="")),]
    }
  } else {
    # label datapoints closest and furthest from line 
    tmp_3 <- sum_all[(sum_all$colsum1 + sum_all$colsum2) > 0.25,]
    tmp_3$y_hat <- (tmp_3$colsum1 * slope_abline) + icept_abline
    tmp_3$diff <- tmp_3$colsum2 - tmp_3$y_hat
    tmp_3 <- tmp_3[order(tmp_3$diff),]
    genes_df <- rbind(tmp_3[1:genes,], tmp_3[(nrow(tmp_3)-genes+1):nrow(tmp_3),])
    tmp_3 <- tmp_3[(tmp_3$colsum1 + tmp_3$colsum2) > 0.7,]
    tmp_3 <- tmp_3[order(abs(tmp_3$diff)),]
    genes_df <- rbind(genes_df, tmp_3[1:genes,])
 
   # label datapoints with highest x+y values
    tmp_3 <- tmp_3[order(tmp_3$colsum1 + tmp_3$colsum2, decreasing=TRUE),]
    genes_df <- rbind(genes_df, tmp_3[1:genes,])
    genes_df <- genes_df %>% distinct(name, colsum1, colsum2)
  }

  gg <- ggplot(sum_all, aes(x = colsum1, y = colsum2)) +
  geom_point(size=2) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour="white"),
    axis.text.x = element_text(size = 20, colour = "black", face = "plain", angle=0, hjust=0.5),
    axis.text.y = element_text(size = 20, colour = "black", face = "plain"),
    axis.title.x = element_text(size = 20, vjust = -1, hjust=0.5),
    axis.title.y = element_text(size = 20, vjust = 1),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.ticks = element_line(colour = "black", size = 1.2), 
    legend.position = "right",
    plot.margin = margin(0.5, 1, 0.5, 0.5, "cm")
  ) +
  scale_x_continuous(name=out_prefix1) +
  scale_y_continuous(name=out_prefix2) +
  geom_abline(slope=slope_abline, intercept=icept_abline, linetype=2, colour="grey", size=1)
  #geom_smooth(method = "lm", se = FALSE, linetype=2, colour="blue", size=1)

  if (length(genes) > 0) {
    gg <- gg +
    geom_text_repel(
    data = as.data.frame(genes_df),
    inherit.aes = FALSE,
    aes(label = name, x=colsum1, y=colsum2),
    size = 5,
    min.segment.length = 0,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = Inf
    )
  }

  ggsave(filename=paste(out_dir, "/ananse_network_dot_", motifs, "_", out_prefix1, "_vs_", 
                        out_prefix2, ".pdf", sep=""), plot=gg, device="pdf", dpi=300, 
                  useDingbats=FALSE)
  return(gg)
}




## function to plot centrality measure across multiple networks as heatmap
plot_compare_heatmap <- function(input_file_list, colnames_column, values_column, 
                                 subset_file=NULL, 
                                 out_prefix="ananse_network_heatmap_jaspar2020", out_dir) {

  condition_list <- list()

  for (i in 1:length(input_file_list)) {
    tmp <- read.csv(input_file_list[i], header=TRUE)
    tmp <- tmp[,c(colnames_column, values_column)]

    condition <- gsub("Ananse_network_", "", input_file_list[i])
    condition <- gsub("_pro_enh_sub_50k.csv", "", condition)
    condition_list <- append(condition_list, condition)
    tmp$condition <- condition

    if (i == 1) {
      heatmap_df <- tmp
    } else {
      heatmap_df <- rbind(heatmap_df, tmp)
    }
  }

  heatmap_df <- as.data.frame(pivot_wider(heatmap_df, names_from=as.name(colnames_column), 
                                          values_from=as.name(values_column), values_fill=0))
  rownames(heatmap_df) <- condition_list

  if (!is.null(subset_file)) {
    heatmap_df <- heatmap_df[,colnames(heatmap_df) %in% subset_file]
  }

  # function to compute z-score
  #zscore_rows <- function(dataframe){
  #  mean_rows <- as.vector(rowMeans(dataframe))
  #  sd_rows <- as.vector(apply(dataframe,1,sd))
  #  zscore <- (dataframe - mean_rows) / sd_rows
  #  return(zscore)
  #}
  #heatmap_df <- zscore_rows(heatmap_df)

  heatmap_df <- log2(heatmap_df-min(heatmap_df)+1)

  pdf(file = paste(out_dir, "/", out_prefix, ".pdf", sep=""),
      width = ncol(heatmap_df) / 15, height=nrow(heatmap_df) * 1)
    pheatmap(heatmap_df, 
    cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, show_colnames=TRUE,
    fontsize = 5, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
  dev.off()

  return(heatmap_df)
}







### Analysis 
conditions_ananse <- read.table(paste(ananse_input_dir, "/conditions_ananse_RNA.txt", sep=""), 
                                sep=" ")

# read in list of transcription factors
factor_file <- read.table(paste(out_dir, "/JASPAR2020.motif2factors.txt",sep=""),
                          header=TRUE)
factors <- factor_file$Factor


# define list of network files to be analysed
file_list <- list.files(path=ananse_output_dir, pattern="_jaspar2020_pro_enh.txt", full.names=TRUE)
file_list <- file_list[grep("*/Network_", file_list)]
print(file_list)


## analyse + subset networks to strongest connections
for (k in seq(1, length(file_list))) {
  print(paste("Analysing", file_list[k], "network"))
  
  n_dirs <- length(strsplit(file_list[k], split="/")[[1]])
  condition_ananse <- matrix(unlist(strsplit(file_list[k], split="/")), nrow=n_dirs)[n_dirs,]
  condition_ananse <- gsub("Network_", "", condition_ananse)
  condition_ananse <- gsub( "_jaspar2020_pro_enh.txt", "", condition_ananse)

  network <- read.table(file_list[k], header=TRUE)
  network$tf_target <- gsub("_", ".", network$tf_target)
  network$tf_target <- gsub("—", "_", network$tf_target)

  out_files <- ananse_net_stats(input_df=network, out_prefix=condition_ananse, factors=factors,
                                out_dir=ananse_output_dir)
  net_out_1 <- out_files[[1]]
  net_out_2 <- out_files[[2]]
  net_out_3 <- out_files[[3]]

  if (k==1) {
    net_out_1_all <- as.data.frame(net_out_1)
    net_out_2_all <- as.data.frame(net_out_2)
    net_out_3_all <- as.data.frame(net_out_3)
  } else {
    net_out_1_all <- as.data.frame(rbind(net_out_1_all, net_out_1))
    net_out_2_all <- as.data.frame(rbind(net_out_2_all, net_out_2))
    net_out_3_all <- as.data.frame(merge(net_out_3_all, net_out_3, by="source"))
  }

  subnet <- subset_network(network, prob_cutoff=0.75, out_prefix=paste(condition_ananse, "_pro_enh_sub", sep=""), 
                           out_dir=ananse_output_dir)
  netgraph <- net_graph(subnet, out_prefix=condition_ananse, out_dir=ananse_output_dir)
}

write.csv(net_out_1_all, paste(ananse_output_dir, "/Ananse_network_overview_all_1.csv", sep=""), 
          row.names=FALSE)
write.csv(net_out_2_all, paste(ananse_output_dir, "/Ananse_network_overview_all_2.csv", sep=""), 
          row.names=FALSE)
write.csv(net_out_3_all, paste(ananse_output_dir, "/Ananse_network_overview_all_3.csv", sep=""), 
          row.names=FALSE)




## compare networks
# only run if there are at least 2 conditions
if (nrow(conditions_ananse) > 1) {

  # Initialise list with conditions that remain to be compared
  conditions_2 <- conditions_ananse

  # do not loop over last condition, because it would not have anything to compare it to
  for (i in seq(1,(nrow(conditions_ananse)-1))) {

    # gradually shrink conditions_2 list to avoid duplicate comparisons
    conditions_2 <- as.data.frame(conditions_2[conditions_2$V1 != conditions_ananse[i,],])
    #print(conditions_2)

    for (k in seq(1,nrow(conditions_2))) {

      # read in sub-networks
      graph1 <- read.csv(paste(ananse_output_dir, "/Ananse_network_centrality_jaspar2020_",
                               conditions_ananse[i,1], ".csv", sep=""), header=TRUE)
      graph2 <- read.csv(paste(ananse_output_dir, "/Ananse_network_centrality_jaspar2020_", 
                               conditions_2[k,1], ".csv", sep=""), header=TRUE)

      # generate dot plot comparing centrality between two networks
      graph1_conn <- t(as.vector(graph1$eigenvector_centrality))
      colnames(graph1_conn) <- graph1$X
      graph2_conn <- t(as.vector(graph2$eigenvector_centrality))
      colnames(graph2_conn) <- graph2$X

      dot_plot <- plot_compare_colsums(graph1_conn, graph2_conn, 
                                       out_prefix1=paste("Eigenvector centrality ",conditions_ananse[i,1],sep=""), 
                                       out_prefix2=paste("Eigenvector centrality ",conditions_2[k,1],sep=""), 
                                       genes=10, motifs="jaspar2020", out_dir=ananse_output_dir)
    }
  }
}


# generate heatmap comparing centrality between all networks
subnet_list <- list.files(path=ananse_output_dir, pattern="Ananse_network_centrality", full.names=TRUE)

graph_heatmap <- plot_compare_heatmap(input_file_list=subnet_list,  
                                      colnames_column="X", values_column="eigenvector_centrality", 
                                      subset_file=factors, out_dir=ananse_output_dir)








