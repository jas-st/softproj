options(warn=-1)

suppressPackageStartupMessages({
  if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
  if(!require("ggtree")) install.packages("ggtree"); library(ggtree)
  if(!require("ggpubr")) install.packages("ggpubr"); library(ggpubr)
  if(!require("ggvenn")) install.packages("ggvenn"); library(ggvenn)
  if(!require("tidytree")) install.packages("tidytree"); library(tidytree)
  if(!require("treeio")) install.packages("treeio"); library(treeio)
  if(!require("phangorn")) install.packages("phangorn"); library(phangorn)
})

#read in the parameters
args <- commandArgs()
treename <- args[6]
raw_results_file<- args[7]
s<-args[8]

#directory for the results
folder<-strsplit(raw_results_file, "/")[[1]][1]

#read in the table
if (s=="true") {
  results_raw_tree <-read_csv(raw_results_file, 
                              col_types = cols(Sat_test_Cassius1 = col_number(), 
                                               Sat_test_Cassius2 = col_number(), 
                                               Chi_test = col_number(), Bowker_test = col_number(), 
                                               Stuart_test = col_number(), Internal_Symmetry = col_number(), 
                                               Proposed_test = col_number(), Alignment_Length = col_number()))
} else {
  results_raw_tree <-read_csv(raw_results_file, 
                              col_types = cols( Bowker_test = col_number(), 
                                                Stuart_test = col_number(), Internal_Symmetry = col_number(), 
                                                Proposed_test = col_number(), Alignment_Length = col_number()))
}

#clean up data and remove all pairs where computation failed/ was not possible
result_tree <- na.omit(results_raw_tree)

#stats for omitted rows
omitted <- round(100-nrow(result_tree)/nrow(results_raw_tree)*100,2)

if (omitted == 100) {
  stop("Error in evaluating the test statistics - All pairs have NA values. Please check the .csv file.")
}

#read in tree
tree <- read.tree(treename)
#midpoint root
tree <- midpoint(tree, node.labels="delete")
#extract labels from tree
seq_labels<-unlist(tree["tip.label"])

#SANITY CHECK: plot tree with node labels
#ggtree(tree) + geom_text(aes(label=node), hjust=-.3)

pairs<-strsplit(str_sub(result_tree$Sequences, start=2, end=-2), ";")
pairs<-lapply(pairs, unlist)
pair_1<-unlist(lapply(pairs,function(x) x[1]))
pair_1<-gsub("\\*","_",pair_1)
pair_2<-unlist(lapply(pairs,function(x) x[2]))
pair_2<-gsub("\\*","_",pair_2)
result_tree$seq_1<-pair_1
result_tree$seq_2<-pair_2

#calculate p values 
QS_pv<-2*(1-pnorm(abs(result_tree$Proposed_test),mean=0,sd=1)) #standardnormalverteilt
Bowker_pv<-pchisq(q=result_tree$Bowker_test, df=6, lower.tail=FALSE)#chi square vert mit df=6 #nicht ganz sicher of two sided oder one sided??????
Stuart_pv<-pchisq(q=result_tree$Stuart_test, df=3, lower.tail=FALSE)#chi square vert mit df=3
IS_pv<-pchisq(q=result_tree$Internal_Symmetry, df=3, lower.tail=FALSE)#chi square vert mit df=3

if (s=="true") Sat_cassius1_pv<-pnorm(result_tree$Sat_test_Cassius1, 0, sd=sqrt(3/result_tree$Alignment_Length), lower.tail=FALSE)
if (s=="true") Sat_cassius2_pv<-pnorm(result_tree$Sat_test_Cassius2, 0, sd=sqrt(3/result_tree$Alignment_Length), lower.tail=FALSE)
if (s=="true") Chi_Test_pv<-pchisq(q=result_tree$Chi_test, df=9, lower.tail=FALSE)

if (s=="true") {
  results_pv <- data.frame(Pair=result_tree$Sequences, pair_1, pair_2, Bowker_pv, Bowker_ts=result_tree$Bowker_test,Stuart_pv,Stuart_ts=result_tree$Stuart_test, 
                               IS_pv, IS_ts=result_tree$Internal_Symmetry, QS_pv, QS_ts=result_tree$Proposed_test, 
                               Sat_cassius1_pv, Sat_test_Cassius1=result_tree$Sat_test_Cassius1, Sat_cassius2_pv, Sat_test_Cassius2=result_tree$Sat_test_Cassius2,
                               Chi_Test_pv, Chi_test=result_tree$Chi_test)
} else {
  results_pv<- data.frame(Pair=result_tree$Sequences, pair_1, pair_2, Bowker_pv, Bowker_ts=result_tree$Bowker_test,Stuart_pv,Stuart_ts=result_tree$Stuart_test, 
                               IS_pv, IS_ts=result_tree$Internal_Symmetry, QS_pv, QS_ts=result_tree$Proposed_test)
}

#######----------------FUNCTIONS--------------------
heat_success_bio<-function(pair_1, pair_2, test_pv, reject=TRUE, length,seq_lables=seq_labels){
  pair_1_tmp<-c(pair_1, pair_2)
  pair_2_tmp<-c(pair_2, pair_1)
  test_pv<-rep(test_pv,2)
  
  test<-data.frame(pair_1=factor(pair_1_tmp),pair_2=factor(pair_2_tmp),test_pv=test_pv)
  
  test<-test%>%mutate(test_rej=factor(ifelse(test_pv>=0.05, 0, 1), levels=c(0,1), ordered=TRUE))
  #shape_values<-ifelse(reject==TRUE, c(1,4), c(4,1))
  if (length < 200) {
    test%>%ggplot(aes(x=pair_1, y=pair_2))+geom_tile(aes(fill=test_pv))+
      #geom_raster(aes(fill = test_pv), interpolate=TRUE) +
      scale_fill_gradient(high="red",low="blue", limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE)) +
      #geom_point(aes(shape=test_rej), size=3)+
      labs(x="Sequence", y="Sequence", fill="p-value", shape="H0")+
      #scale_shape_manual(labels=c("keep H0", "reject H0"), values=c(1,4), drop=F) +
      #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  } else {
    test%>%ggplot(aes(x=pair_1, y=pair_2))+geom_tile(aes(fill=test_pv))+
      #geom_raster(aes(fill = test_pv), interpolate=TRUE) +
      scale_fill_gradient(high="red",low="blue", limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE)) +
      #geom_point(aes(shape=test_rej), size=3)+
      labs(x="Sequence", y="Sequence", fill="p-value", shape="H0")+
      #scale_shape_manual(labels=c("keep H0", "reject H0"), values=c(1,4), drop=F) +
      #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
            axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
  }
}

heat_success<-function(pair_1, pair_2, test_pv, reject=TRUE, seq_lables=seq_labels){
  pair_1_tmp<-c(pair_1, pair_2)
  pair_2_tmp<-c(pair_2, pair_1)
  test_pv<-rep(test_pv,2)
  
  test<-data.frame(pair_1=factor(pair_1_tmp),pair_2=factor(pair_2_tmp),test_pv=test_pv)
  
  test<-test%>%mutate(test_rej=factor(ifelse(test_pv>=0.05, 0, 1), levels=c(0,1), ordered=TRUE))
  #shape_values<-ifelse(reject==TRUE, c(1,4), c(4,1))
  test%>%ggplot(aes(x=pair_1, y=pair_2))+geom_tile(aes(fill=test_pv))+
    scale_fill_gradient(high="red",low="blue", limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE)) +
    geom_point(aes(shape=test_rej), size=3)+
    labs(x="",y="sequence", fill="p-value", shape="H0")+
    scale_shape_manual(labels=c("retain", "reject"), values=c(1,4), drop=F)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank() 
    )
}

get_longest_path <- function(tree) {
  Nleaf <- length(unlist(tree['tip.label']))
  tree_root <- getRoot(tree)
  max_path <- 0
  for (i in 1:Nleaf) {
    tip_path <- length(get.path(tree, i, tree_root))
    if (tip_path > max_path) {
      max_path <- tip_path
    }
  }
  return (max_path)
}

edges_rejected_freq <- function(tree, pair_1, pair_2, test_pv, freq_true){
  pair_1_index<-mapply(function(x) which(x==seq_labels), pair_1)
  pair_2_index<-mapply(function(x) which(x==seq_labels), pair_2)
  paths<-mapply(function(x,y) get.path(tree, x,y), pair_1_index, pair_2_index)
  rej_paths<-mapply(function(x) ifelse(x>=0.05,0,1), test_pv )
  edge<-data.frame(tree$edge)
  colnames(edge)=c("parent", "node")
  
  num_edges<-length(edge$parent)
  edges_rej<-c(rep(0,num_edges))
  edges_freq<-c(rep(0,num_edges))
  for(i in 1:num_edges){
    p<-edge$parent[i]
    n<-edge$node[i]
    e<-c(p,n)
    for(j in 1:length(paths)){
      if(p %in% paths[[j]] && n %in% paths[[j]]){
        edges_rej[i]<-edges_rej[i]+rej_paths[j]
        edges_freq[i]<-edges_freq[i]+1
      }
    }
    if(edges_freq[i]!=0 && edges_rej[i]!=0) edges_rej[i]=edges_rej[i]/edges_freq[i]
  }
  edge["edge_rej"]<-round(edges_rej,2)
  return(edge)
}


colored_tree <- function(tree, d1){
  p<-ggtree(tree, size=0.8, layout="slanted") + geom_tiplab(colour="black", size=3.5, align = FALSE)+
    labs(title="H0 rejected (freq) on each edge")
  p%<+% d1 + geom_label(aes(x=branch, label=edge_rej)) + 
    aes(color=edge_rej) + 
    scale_colour_gradient2(low="blue", mid="#FF0099", high="red",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    #scale_colour_gradient2(low="#edd96e", mid="#e89b53", high="#7a002d",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    labs(color="H0 rejected (freq)")
}

colored_tree_clean <- function(tree, d1){
  p<-ggtree(tree, size=0.8, layout="slanted") + #geom_tiplab(colour="black", size=3.5, align = FALSE)+
    labs(title="H0 rejected (freq) on each edge")
  p%<+% d1 + #geom_label(aes(x=branch, label=edge_rej)) + 
    aes(color=edge_rej) + 
    scale_colour_gradient2(low="blue", mid="#FF0099", high="red",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    #scale_colour_gradient2(low="#edd96e", mid="#e89b53", high="#7a002d",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    labs(color="H0 rejected (freq)")
}


compressed_tree <- function(tree, col_tree, freq_threshold, max_root_distance) {
  col_tree<-col_tree + geom_tiplab(colour="black", size=3, align = FALSE) 
  d1 <- data.frame(col_tree[["data"]][["parent"]],col_tree[["data"]][["node"]])
  colnames(d1) <- c("parent", "node")
  d1["edge_rej"] <- col_tree[["data"]][["edge_rej"]]
  Nleaf <- length(unlist(tree['tip.label']))
  Nint <- as.numeric(tree['Nnode'])
  tree_root <- getRoot(tree)
  collapse_list_cleaned <- c()
  collapse_list <- c()

  for (node_tmp in (Nleaf+1):(Nleaf+Nint)) {
    path_to_root <- get.path(tree,tree_root,node_tmp)
    path_to_root_length <- length(path_to_root)
    children <- Children(tree, node_tmp)
    warden <- TRUE
  
    for (nd in collapse_list) {
'      path_between <- length(get.path(tree, nd, node_tmp))
      path_nd <- length(get.path(tree, tree_root, nd))

      if (path_between < path_to_root_length & path_nd < path_to_root_length) {
        warden <- FALSE
      }'
      if (nd %in% path_to_root) {
        warden <- FALSE
        break
      }
      
    }
  
    if (warden) {
      for (child in children) {
        child_freq <- subset(d1, node == child, select = edge_rej)[1,1]
      
        if (path_to_root_length < max_root_distance) {
            warden <- FALSE
            break
        }
        else {
            if (child_freq > freq_threshold) {
              warden <- FALSE
              break
            }
          }
        #if (child_freq > freq_threshold || path_to_root_length < max_root_distance) {
        #  warden <- FALSE
        #  break
        #}
      }
    }
  
    if (warden & node_tmp != tree_root) {
      col_tree <- col_tree %>% collapse(node = node_tmp)
      collapse_list <- append(collapse_list, node_tmp)
    }
  }
  #remove descendants of collapsed nodes
  collapse_list_cleaned <- collapse_list
  for (nd in collapse_list) {
    root_path <- get.path(tree, nd, tree_root)
    if (any(collapse_list[collapse_list!=nd] %in% get.path(tree, nd, tree_root))) {
      collapse_list_cleaned <- collapse_list_cleaned[collapse_list_cleaned!=nd]
    }
  }
  
  return (list(col_tree, collapse_list_cleaned))
}

#tree length
Nleaf <- length(unlist(tree['tip.label']))
p2_col <- c()


######-------------PLOTS-------------------###############
#---------------------------------------BOWKER----------------------------------
if (Nleaf > 50) {
  #threshold for the collapsing
  a <- ceiling(0.28*get_longest_path(tree))
  
  h1<-heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv, reject=TRUE, Nleaf) +
    ggtitle(paste("Bowker Test, omitted pairs: ", toString(omitted), "%"))
  Bowker_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv, TRUE)
  p1<-colored_tree_clean(tree, Bowker_edge_freq)
  p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Bowker Test, Midpoint root")
  p1_labeled <- p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
  
  if (Nleaf < 200) {
    p1<-p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
  }
  
  p2_list <- compressed_tree(tree,p1,0.3,a-1)
  p2 <- p2_list[[1]]
  p2_col <- p2_list[[2]]
  
  if (length(p2_col) > 0) {
    p2 <- p2 + 
      labs(title = "Bowker Test, Compressed Tree, Midpoint root") +
      geom_label(data=subset(p2$data, branch.length > 0.02 & !(node %in% p2_col)), aes(x=branch, label=edge_rej)) +
      geom_point2(data=subset(p2$data,node %in% p2_col),shape=21, size = 7, fill='yellow') + 
      geom_text(data=subset(p2$data,node %in% p2_col), aes(label=node), fontface='bold') 
  
    subtrees <- list()
    for (cl in p2_col) {
      b <- viewClade(p1_labeled, cl) + 
        labs(color="", title = paste("Node ",toString(cl))) + 
        theme(legend.position = "none")
      subtrees <- append(subtrees, list(b))
    }
    p3 <- gridExtra::grid.arrange(grobs = subtrees, top="Collapsed Nodes")
  }
  
  p1 <- p1 + geom_label(data=subset(p1$data, branch.length > 0.005), aes(x=branch, label=edge_rej))
} else {
  Bowker_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv, TRUE)
  p1<-colored_tree(tree, Bowker_edge_freq)
  p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Bowker Test, Midpoint root")
  h1<-heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv, reject=TRUE) +
    ggtitle(paste("Bowker Test, omitted pairs: ", toString(omitted), "%"))
}
pdf(paste(folder, "plot_Bowker_test.pdf", sep="/"),width = 13, height = 14)
h1
p1
if (Nleaf > 50 & length(p2_col)>0) {
  print(p2)
  plot(p3)
} else {
  txt1="No compressed tree."
  txt2="Less than 50 species or no collapsed clades."
  plot.new()
  text(x=.1,y=.8,txt1, cex=0.7)
  text(x=.1, y=.7, txt2, cex=0.7)
}
dev.off()

#---------------------------------------STUART----------------------------------
if (Nleaf > 50) {
  #threshold for the collapsing
  a <- ceiling(0.28*get_longest_path(tree))
  
  h1<-heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv, reject=TRUE, Nleaf) +
    ggtitle(paste("Stuart Test, omitted pairs: ", toString(omitted), "%"))
  Stuart_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv, TRUE)
  p1<-colored_tree_clean(tree, Stuart_edge_freq)
  p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Stuart Test, Midpoint root") 
  p1_labeled <- p1 + geom_tiplab(data=subset(p1$data, edge_rej > 0.5), colour="black", size=3, align = FALSE) 
  if (Nleaf < 200) {
    p1<-p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
  }
  
  p2_list <- compressed_tree(tree,p1,0.3,a-1)
  p2 <- p2_list[[1]]
  p2_col <- p2_list[[2]]
  
  if (length(p2_col) > 0) {
    p2 <- p2 + 
      labs(title = "Stuart Test, Compressed Tree, Midpoint root") +
      geom_label(data=subset(p2$data, branch.length > 0.02 & !(node %in% p2_col)), aes(x=branch, label=edge_rej)) +
      geom_point2(data=subset(p2$data,node %in% p2_col),shape=21, size = 7, fill='yellow') + 
      geom_text(data=subset(p2$data,node %in% p2_col), aes(label=node), fontface='bold') 
  
    subtrees <- list()
    for (cl in p2_col) {
      b <- viewClade(p1_labeled, cl) + 
        labs(color="", title = paste("Node ",toString(cl))) + 
        theme(legend.position = "none")
      subtrees <- append(subtrees, list(b))
    }
    p3 <- gridExtra::grid.arrange(grobs = subtrees, top="Collapsed Nodes")
  }
  
  p1 <- p1 + geom_label(data=subset(p1$data, branch.length > 0.005), aes(x=branch, label=edge_rej))
} else {
  Stuart_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv, TRUE)
  p1<-colored_tree(tree, Stuart_edge_freq)
  p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Stuart Test, Midpoint root")
  h1<-heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv, reject=TRUE) +
    ggtitle(paste("Stuart Test, omitted pairs: ", toString(omitted), "%"))
}
pdf(paste(folder, "plot_Stuart_test.pdf", sep="/"),width = 13, height = 14)
h1
p1
if (Nleaf > 50 & length(p2_col)>0) {
  print(p2)
  plot(p3)
} else {
  txt1="No compressed tree."
  txt2="Less than 50 species or no collapsed clades."
  plot.new()
  text(x=.1,y=.8,txt1, cex=0.7)
  text(x=.1, y=.7, txt2, cex=0.7)
}
dev.off()

#---------------------------------------IS----------------------------------
if (Nleaf > 50) {
  #threshold for the collapsing
  a <- ceiling(0.28*get_longest_path(tree))
  
  h1<-heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv, reject=TRUE, Nleaf) +
    ggtitle(paste("IS Test, omitted pairs: ", toString(omitted), "%"))
  IS_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv, TRUE)
  p1<-colored_tree_clean(tree, IS_edge_freq)
  p1<-p1 + labs(title = "H0 rejected (freq) on each edge, IS Test, Midpoint root") 
  p1_labeled <- p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
  if (Nleaf < 200) {
    p1<-p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
  }
  
  p2_list <- compressed_tree(tree,p1,0.3,a-1)
  p2 <- p2_list[[1]]
  p2_col <- p2_list[[2]]
  
  if (length(p2_col) > 0) {
    p2 <- p2 + 
      labs(title = "IS Test, Compressed Tree, Midpoint root") +
      geom_label(data=subset(p2$data, branch.length > 0.02 & !(node %in% p2_col)), aes(x=branch, label=edge_rej)) +
      geom_point2(data=subset(p2$data,node %in% p2_col),shape=21, size = 7, fill='yellow') + 
      geom_text(data=subset(p2$data,node %in% p2_col), aes(label=node), fontface='bold') 
  
    subtrees <- list()
    for (cl in p2_col) {
      b <- viewClade(p1_labeled, cl) + 
        labs(color="", title = paste("Node ",toString(cl))) + 
        theme(legend.position = "none")
      subtrees <- append(subtrees, list(b))
    }
    p3 <- gridExtra::grid.arrange(grobs = subtrees, top="Collapsed Nodes")
  }
  
  p1 <- p1 + geom_label(data=subset(p1$data, branch.length > 0.005), aes(x=branch, label=edge_rej))
} else {
  IS_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv, TRUE)
  p1<-colored_tree(tree, IS_edge_freq)
  p1<-p1 + labs(title = "H0 rejected (freq) on each edge, IS Test, Midpoint root")
  h1<-heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv, reject=TRUE) +
    ggtitle(paste("IS Test, omitted pairs: ", toString(omitted), "%"))
}
pdf(paste(folder, "plot_IS_test.pdf", sep="/"),width = 13, height = 14)
h1
p1
if (Nleaf > 50 & length(p2_col)>0) {
  print(p2)
  plot(p3)
} else {
  txt1="No compressed tree."
  txt2="Less than 50 species or no collapsed clades."
  plot.new()
  text(x=.1,y=.8,txt1, cex=0.7)
  text(x=.1, y=.7, txt2, cex=0.7)
}
dev.off()

#---------------------------------------QS----------------------------------
if (Nleaf > 50) {
  #threshold for the collapsing
  a <- ceiling(0.28*get_longest_path(tree))
  
  h1<-heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv, reject=TRUE, Nleaf) +
    ggtitle(paste("QS Test, omitted pairs: ", toString(omitted), "%"))
  QS_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv, TRUE)
  p1<-colored_tree_clean(tree, QS_edge_freq)
  p1<-p1 + labs(title = "H0 rejected (freq) on each edge, QS Test, Midpoint root") 
  p1_labeled <- p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
  if (Nleaf < 200) {
    p1<-p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
  }
  
  p2_list <- compressed_tree(tree,p1,0.3,a-1)
  p2 <- p2_list[[1]]
  p2_col <- p2_list[[2]]
  
  if (length(p2_col) > 0) {
    p2 <- p2 + 
      labs(title = "QS Test, Compressed Tree, Midpoint root") +
      geom_label(data=subset(p2$data, branch.length > 0.02 & !(node %in% p2_col)), aes(x=branch, label=edge_rej)) +
      geom_point2(data=subset(p2$data,node %in% p2_col),shape=21, size = 7, fill='yellow') + 
      geom_text(data=subset(p2$data,node %in% p2_col), aes(label=node), fontface='bold') 
  
    subtrees <- list()
    for (cl in p2_col) {
      b <- viewClade(p1_labeled, cl) + 
        labs(color="", title = paste("Node ",toString(cl))) + 
        theme(legend.position = "none")
      subtrees <- append(subtrees, list(b))
    }
    p3 <- gridExtra::grid.arrange(grobs = subtrees, top="Collapsed Nodes")
  }
  
  p1 <- p1 + geom_label(data=subset(p1$data, branch.length > 0.005), aes(x=branch, label=edge_rej))
} else {
  QS_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv, TRUE)
  p1<-colored_tree(tree, QS_edge_freq)
  p1<-p1 + labs(title = "H0 rejected (freq) on each edge, QS Test, Midpoint root")
  h1<-heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv, reject=TRUE) +
    ggtitle(paste("QS Test, omitted pairs: ", toString(omitted), "%"))
}
pdf(paste(folder, "plot_QS_test.pdf", sep="/"),width = 13, height = 14)
h1
p1
if (Nleaf > 50 & length(p2_col)>0) {
  print(p2)
  plot(p3)
} else {
  txt1="No compressed tree."
  txt2="Less than 50 species or no collapsed clades."
  plot.new()
  text(x=.1,y=.8,txt1, cex=0.7)
  text(x=.1, y=.7, txt2, cex=0.7)
}
dev.off()

#---------------------------------------SAT1----------------------------------
if (s=="true") {
  if (Nleaf > 50) {
    #threshold for the collapsing
    a <- ceiling(0.28*get_longest_path(tree))
    
    h1<-heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv, reject=TRUE, Nleaf) +
      ggtitle(paste("Saturation Test Cassius 1, omitted pairs: ", toString(omitted), "%"))
    Sat1_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv, TRUE)
    p1<-colored_tree_clean(tree, Sat1_edge_freq)
    p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Saturation Test Cassius 1, Midpoint root")
    p1_labeled <- p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
    if (Nleaf < 200) {
      p1<-p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
    }
    
    p2_list <- compressed_tree(tree,p1,0.3,a-1)
    p2 <- p2_list[[1]]
    p2_col <- p2_list[[2]]
    
    if (length(p2_col) > 0) {
      p2 <- p2 + 
        labs(title = "Saturation Test Cassius 1, Compressed Tree, Midpoint root") +
        geom_label(data=subset(p2$data, branch.length > 0.02 & !(node %in% p2_col)), aes(x=branch, label=edge_rej)) +
        geom_point2(data=subset(p2$data,node %in% p2_col),shape=21, size = 7, fill='yellow') + 
        geom_text(data=subset(p2$data,node %in% p2_col), aes(label=node), fontface='bold') 
    
      subtrees <- list()
      for (cl in p2_col) {
        b <- viewClade(p1_labeled, cl) + 
          labs(color="", title = paste("Node ",toString(cl))) + 
          theme(legend.position = "none")
        subtrees <- append(subtrees, list(b))
      }
      p3 <- gridExtra::grid.arrange(grobs = subtrees, top="Collapsed Nodes")
    }
    p1 <- p1 + geom_label(data=subset(p1$data, branch.length > 0.005), aes(x=branch, label=edge_rej))
    
  } else {
    Sat1_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv, TRUE)
    p1<-colored_tree(tree, Sat1_edge_freq)
    p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Saturation Test Cassius 1, Midpoint root")
    h1<-heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv, reject=TRUE) +
      ggtitle(paste("Saturation Test Cassius 1, omitted pairs: ", toString(omitted), "%"))
  }
  pdf(paste(folder, "plot_Sat_Cassius1_test.pdf", sep="/"),width = 13, height = 14)
  print(h1)
  print(p1)
  if (Nleaf > 50 & length(p2_col)>0) {
    print(p2)
    plot(p3)
  } else {
    txt1="No compressed tree."
    txt2="Less than 50 species or no collapsed clades."
    plot.new()
    text(x=.1,y=.8,txt1, cex=0.7)
    text(x=.1, y=.7, txt2, cex=0.7)
    }
  
  dev.off()
}

#---------------------------------------SAT2----------------------------------
if(s=="true") {
  if (Nleaf > 50) {
    #threshold for the collapsing
    a <- ceiling(0.28*get_longest_path(tree))
    
    h1<-heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius2_pv, reject=TRUE, Nleaf) +
      ggtitle(paste("Saturation Test Cassius 2, omitted pairs: ", toString(omitted), "%"))
    Sat2_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius2_pv, TRUE)
    p1<-colored_tree_clean(tree, Sat2_edge_freq)
    p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Saturation Test Cassius 2, Midpoint root") 
    p1_labeled <- p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
    if (Nleaf < 200) {
      p1<-p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
    } 
    
    p2_list <- compressed_tree(tree,p1,0.3,a-1)
    p2 <- p2_list[[1]]
    p2_col <- p2_list[[2]]
    
    if (length(p2_col) > 0) {
      p2 <- p2 + 
        labs(title = "Saturation Test Cassius 2, Compressed Tree, Midpoint root") +
        geom_label(data=subset(p2$data, branch.length > 0.02 & !(node %in% p2_col)), aes(x=branch, label=edge_rej)) +
        geom_point2(data=subset(p2$data,node %in% p2_col),shape=21, size = 7, fill='yellow') + 
        geom_text(data=subset(p2$data,node %in% p2_col), aes(label=node), fontface='bold') 
  
      subtrees <- list()
      for (cl in p2_col) {
        b <- viewClade(p1_labeled, cl) + 
          labs(color="", title = paste("Node ",toString(cl))) + 
          theme(legend.position = "none")
        subtrees <- append(subtrees, list(b))
      }
      p3 <- gridExtra::grid.arrange(grobs = subtrees, top="Collapsed Nodes")
    }
    
    p1 <- p1 + geom_label(data=subset(p1$data, branch.length > 0.005), aes(x=branch, label=edge_rej))
    
  } else {
    Sat2_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius2_pv, TRUE)
    p1<-colored_tree(tree, Sat2_edge_freq)
    p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Saturation Test Cassius 2, Midpoint root")
    h1<-heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius2_pv, reject=TRUE) +
      ggtitle(paste("Saturation Test Cassius 2 Test, omitted pairs: ", toString(omitted), "%"))
  }
  pdf(paste(folder, "plot_Sat_Cassius2_test.pdf", sep="/"),width = 13, height = 14)
  print(h1)
  print(p1)
  if (Nleaf > 50 & length(p2_col)>0) {
    print(p2)
    plot(p3)
  } else {
    txt1="No compressed tree."
    txt2="Less than 50 species or no collapsed clades."
    plot.new()
    text(x=.1,y=.8,txt1, cex=0.7)
    text(x=.1, y=.7, txt2, cex=0.7)
  }
  dev.off()
}

#---------------------------------------CHI----------------------------------
if (s=="true") {
  if (Nleaf > 50) {
    #threshold for the collapsing
    a <- ceiling(0.28*get_longest_path(tree))
    
    h1<-heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$Chi_Test_pv, reject=TRUE, Nleaf) +
      ggtitle(paste("Chi Test, omitted pairs: ", toString(omitted), "%"))
    Chi_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Chi_Test_pv, TRUE)
    p1<-colored_tree_clean(tree, Chi_edge_freq)    
    p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Chi Test, Midpoint root") 
    p1_labeled <- p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
    if (Nleaf < 200) {
      p1<-p1 + geom_tiplab(colour="black", size=3, align = FALSE) 
    }
    
    p2_list <- compressed_tree(tree,p1,0.3,a-1)
    p2 <- p2_list[[1]]
    p2_col <- p2_list[[2]]
    
    if (length(p2_col) > 0) {
      p2 <- p2 + 
        labs(title = "Chi Test, Compressed Tree, Midpoint root") +
        geom_label(data=subset(p2$data, branch.length > 0.02 & !(node %in% p2_col)), aes(x=branch, label=edge_rej)) +
        geom_point2(data=subset(p2$data,node %in% p2_col),shape=21, size = 7, fill='yellow') + 
        geom_text(data=subset(p2$data,node %in% p2_col), aes(label=node), fontface='bold') 
    
      subtrees <- list()
      for (cl in p2_col) {
        b <- viewClade(p1_labeled, cl) + 
          labs(color="", title = paste("Node ",toString(cl))) + 
          theme(legend.position = "none")
        subtrees <- append(subtrees, list(b))
      }
      p3 <- gridExtra::grid.arrange(grobs = subtrees, top="Collapsed Nodes")
    }
    
    p1 <- p1 + geom_label(data=subset(p1$data, branch.length > 0.005), aes(x=branch, label=edge_rej))
  } else {
    Chi_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Chi_Test_pv, TRUE)
    p1<-colored_tree(tree, Chi_edge_freq) 
    p1<-p1 + labs(title = "H0 rejected (freq) on each edge, Chi Test, Midpoint root")
    h1<-heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Chi_Test_pv, reject=TRUE) +
      ggtitle(paste("Chi Test, omitted pairs: ", toString(omitted), "%"))
  }
  pdf(paste(folder, "plot_Chi_test.pdf", sep="/"),width = 13, height = 14)
  print(h1)
  print(p1)
  if (Nleaf > 50 & length(p2_col)>0) {
    print(p2)
    plot(p3)
  } else {
    txt1="No compressed tree."
    txt2="Less than 50 species or no collapsed clades."
    plot.new()
    text(x=.1,y=.8,txt1, cex=0.7)
    text(x=.1, y=.7, txt2, cex=0.7)
  }
  dev.off()
}

#---------------------------------------VENN DIAGRAMM----------------------------------
# create dataframe with number of rejects for each pair
Bowker_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,FALSE,TRUE), Bowker_pv)
Stuart_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,FALSE,TRUE), Stuart_pv)
IS_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,FALSE,TRUE), IS_pv)
QS_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,FALSE,TRUE), QS_pv)

venn_data<-data.frame(Pair=result_tree$Sequences, Bowker=Bowker_rej, Stuart=Stuart_rej, IS=IS_rej, QS=QS_rej)

pdf(paste(folder, "venn_diag.pdf", sep="/"),width = 13, height = 14)
ggvenn(venn_data, c("Bowker", "Stuart", "QS"))
dev.off()

#-----------------------------RESULTS TESTS STATS + P VALUES------------------------
results_file<-basename(raw_results_file)
pat <- "(.*?_.*?)_(.*)"
results_file <- sub(pat, "\\2", results_file)
results_file <- paste("results", results_file, sep="_")
write.csv(results_pv, file = paste(folder, results_file, sep="/"), row.names = FALSE)


#------------------------------TESTS RETAIN/REJECT----------------------------
Bowker_success<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,"retain","reject"), Bowker_pv)
Stuart_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Stuart_pv)
IS_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), IS_pv)
QS_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), QS_pv)
Reversibility_test<-data.frame(Pair=results_pv$Pair,
                               Bowker_success, Stuart_success, IS_success, QS_success)

write.csv(Reversibility_test, paste(folder, "results_rev_test.csv", sep="/"),row.names = FALSE)

#------------------------------SAT TESTS RETAIN/REJECT----------------------------
if (s=="true") {
  Sat_cassius1_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Sat_cassius1_pv)
  Sat_cassius2_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Sat_cassius2_pv)
  Chi_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Chi_Test_pv)
  Saturation_test<-data.frame(Pair=results_pv$Pair,
                              Sat_cassius1_success, Sat_cassius2_success, Chi_success)
  write.csv(Saturation_test, paste(folder, "results_sat_test.csv", sep="/"),row.names = FALSE)
}

#---------------------------------------ANNOTATED TREE----------------------------------

#Bowker_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv, TRUE)
#Stuart_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv, TRUE)
#IS_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv, TRUE)
#QS_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv, TRUE)
#if (s=="true") Sat1_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv, TRUE)
#if (s=="true") Sat2_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius2_pv, TRUE)
#if (s=="true") Chi_edge_freq<-edges_rejected_freq(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Chi_Test_pv, TRUE)

if (s=="true") {
  labels_dataframe<-data.frame(parent=Bowker_edge_freq$parent,node=Bowker_edge_freq$node, Bowker_freq=Bowker_edge_freq$edge_rej,
                               Stuart_freq=Stuart_edge_freq$edge_rej, IS_freq=IS_edge_freq$edge_rej,
                               QS_freq=QS_edge_freq$edge_rej, Sat1_req=Sat1_edge_freq$edge_rej,
                               Sat2_freq=Sat2_edge_freq$edge_rej, Chi_rej=Chi_edge_freq$edge_rej)
} else {
  labels_dataframe<-data.frame(parent=Bowker_edge_freq$parent,node=Bowker_edge_freq$node, Bowker_freq=Bowker_edge_freq$edge_rej,
                               Stuart_freq=Stuart_edge_freq$edge_rej, IS_freq=IS_edge_freq$edge_rej,
                               QS_freq=QS_edge_freq$edge_rej)
}

y<-full_join(as.tibble(tree), labels_dataframe, by = c("parent",'node'))


outputfile<-paste(basename(treename), ".tree", sep="")
write.beast(as.treedata(y), file=paste(folder, outputfile, sep="/"))