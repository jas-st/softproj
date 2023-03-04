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
results_file<-args[8]
s<-args[9]

#directory for the results
folder<-strsplit(results_file, "/")[[1]][1]

raw_results_file <- "brown_results.csv"
raw_results_file <- "anderson_results.csv"
raw_results_file <- "faircloth_results.csv"
raw_results_file <- "worobey_results.csv"
treename <- "brown.treefile"
treename <- "anderson.treefile"
treename <- "faircloth.treefile"
treename <- "worobey.treefile"
folder <- "brown_results"
s <- "true"

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

#in case smth failed w the c++ reading in - remove all paris with ; or end

#TODO: stats for omitted rows
omitted <- round(100-nrow(result_tree)/nrow(results_raw_tree)*100,2)

#TODO: remove
result_tree<-result_tree%>%dplyr::rename(Sequences="Sequences compared")

#read in tree
tree <- read.tree(treename)
#midpoint root
midpoint_tree <- midpoint(tree)
#extract labels from tree
seq_labels<-unlist(tree["tip.label"])

#SANITY CHECK tree plot with labels
#p <- ggtree(tree) + geom_tiplab()

pairs<-strsplit(str_sub(result_tree$Sequences, start=2, end=-2), ";")
pairs<-lapply(pairs, unlist)
pair_1<-unlist(lapply(pairs,function(x) x[1]))
pair_1<-gsub("\\*","_",pair_1)
pair_2<-unlist(lapply(pairs,function(x) x[2]))
pair_2<-gsub("\\*","_",pair_1)
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


heat_success_bio<-function(pair_1, pair_2, test_pv, reject=TRUE, seq_lables=seq_labels){
  pair_1_tmp<-c(pair_1, pair_2)
  pair_2_tmp<-c(pair_2, pair_1)
  test_pv<-rep(test_pv,2)
  
  test<-data.frame(pair_1=factor(pair_1_tmp),pair_2=factor(pair_2_tmp),test_pv=test_pv)
  
  test<-test%>%mutate(test_rej=factor(ifelse(test_pv>=0.05, 0, 1), levels=c(0,1), ordered=TRUE))
  #shape_values<-ifelse(reject==TRUE, c(1,4), c(4,1))
  test%>%ggplot(aes(x=pair_1, y=pair_2))+geom_tile(aes(fill=test_pv))+
    #geom_raster(aes(fill = test_pv), interpolate=TRUE) +
    scale_fill_gradient(high="red",low="blue", limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE)) +
    #geom_point(aes(shape=test_rej), size=3)+
    labs(x="Sequence", y="Sequence", fill="p-value", shape="H0")+
    #scale_shape_manual(labels=c("keep H0", "reject H0"), values=c(1,4), drop=F) +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(),axis.ticks.y=element_blank())
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


colored_tree <- function(tree, pair_1, pair_2, test_pv){
  d1<-edges_rejected_freq(tree, pair_1, pair_2, test_pv, TRUE)
  p<-ggtree(tree, size=0.8, layout="slanted") + geom_tiplab(colour="black", size=3.5, align = FALSE)+
    labs(title="H0 rejected (freq) on each edge")
  p%<+% d1 + geom_label(aes(x=branch, label=edge_rej)) + 
    aes(color=edge_rej) + 
    scale_colour_gradient2(low="blue", mid="#FF0099", high="red",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    #scale_colour_gradient2(low="#edd96e", mid="#e89b53", high="#7a002d",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    labs(color="H0 rejected (freq)")
}

colored_tree_clean <- function(tree, pair_1, pair_2, test_pv){
  d1<-edges_rejected_freq(tree, pair_1, pair_2, test_pv, TRUE)
  p<-ggtree(tree, size=0.8, layout="slanted") + #geom_tiplab(colour="black", size=3.5, align = FALSE)+
    labs(title="H0 rejected (freq) on each edge")
  p%<+% d1 + #geom_label(aes(x=branch, label=edge_rej)) + 
    aes(color=edge_rej) + 
    scale_colour_gradient2(low="blue", mid="#FF0099", high="red",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    #scale_colour_gradient2(low="#edd96e", mid="#e89b53", high="#7a002d",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    labs(color="H0 rejected (freq)")
}


compressed_tree <- function(tree, col_tree, freq_threshold, max_root_distance) {
  d1 <- data.frame(col_tree[["data"]][["parent"]],col_tree[["data"]][["node"]])
  colnames(d1) <- c("parent", "node")
  d1["edge_rej"] <- col_tree[["data"]][["edge_rej"]]
  Nleaf <- length(unlist(tree['tip.label']))
  Nint <- as.numeric(tree['Nnode'])
  tree_root <- getRoot(tree)
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
          if (child_freq > 0.05) {
            warden <- FALSE
            break
          } 
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
  return (list(col_tree, collapse_list))
}

a <- ceiling(0.28*get_longest_path(tree))

p1<-colored_tree_clean(tree, results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv)
p1 + labs(title = "H0 rejected (freq) on each edge, QS Test, Arbitrary root")
p2<-colored_tree_clean(midpoint_tree, results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv)
p2 + labs(title = "H0 rejected (freq) on each edge, Bowker Test, Midpoint root")
pqls <- compressed_tree(tree,p1,0.3,a)
pqls <- compressed_tree(midpoint_tree,p2,0.3,a-3)
pq <- pqls[[1]]
pql <- pqls[[2]]

pq + 
  labs(title = "Bowker Test, Compressed Tree, Arbitrary root") +
  #labs(title = "Bowker Test, Compressed Tree, Midpoint root") +
  geom_label(data=subset(pq$data, branch.length > 0.005 & !(node %in% pql)), aes(x=branch, label=edge_rej)) +
  geom_point2(data=subset(pq$data,node %in% pql),shape=21, size = 7, fill='yellow') + 
  geom_text(data=subset(pq$data,node %in% pql), aes(label=node), fontface='bold') 

subtrees <- list()
for (cl in pql) {
  b <- viewClade(p2, cl) + 
    labs(color="", title = paste("Node ",toString(cl))) + 
    theme(legend.position = "none")
  subtrees <- append(subtrees, list(b))
}
gridExtra::grid.arrange(grobs = subtrees)


ggtree(midpoint_tree) + geom_text(aes(label=node), hjust=-.3)





pdf(paste(folder, "plot_Bowker_test.pdf", sep="/"),width = 30, height = 30)

p <- colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv)
heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv, reject=TRUE) +
  ggtitle(paste("Bowker Test, omitted pairs: ", toString(omitted), "%"))
  

#annotate_figure(p,top=text_grob("Bowker Test", face = "bold", size = 14))
dev.off()

pdf(paste(folder, "plot_Stuart_test.pdf", sep="/"),width = 13, height = 14)

p<-colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv)
heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv, reject=TRUE)  +
  ggtitle("Stuart Test")
results_pv_all%>%ggplot(aes(x=Stuart_ts))+
  geom_histogram(bins = 40, aes(y = after_stat(!!str2lang("density"))))+
  stat_function(fun = dchisq, args = list(df=3))

#annotate_figure(p,top=text_grob("Stuart Test", face = "bold", size = 14))
dev.off()

# pdf(paste(folder, "plot_IS_test.pdf", sep="/"),width = 13, height = 14)
# p<-ggarrange(
#   p<-colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv),
heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv, reject=TRUE) +
  ggtitle("IS Test")
#   results_pv_all%>%ggplot(aes(x=IS_ts))+
#     geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#     stat_function(fun = dchisq, args = list(df=3))+
#     facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")+labs(x="test statistic")
# )
# annotate_figure(p,top=text_grob("IS Test", face = "bold", size = 14))
# dev.off()
# 
# pdf(paste(folder, "plot_QS_test.pdf", sep="/"),width = 13, height = 14)
# p<-ggarrange(
#   p<-colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv),
   heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv, reject=TRUE) +
     ggtitle("QS Test")
#   results_pv_all%>%ggplot(aes(x=QS_ts))+
#     geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#     stat_function(fun = dnorm, args = list(mean = 0, sd = 1))+
#     facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# )
# annotate_figure(p,top=text_grob("QS Test", face = "bold", size = 14))
# dev.off()
# 
# pdf(paste(folder, "plot_Sat_Cassius1_test.pdf", sep="/"),width = 13, height = 14)
# p<-ggarrange(
#   colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv),
#   heat_success_bio(results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv, reject=TRUE)
#   results_pv_all%>%ggplot(aes(x=Sat_test_Cassius1))+
#     geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#     stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(3/n)))+
#     facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# )
# annotate_figure(p,top=text_grob("Saturation Test Cassius 1", face = "bold", size = 14))
# dev.off()
# 
# pdf(paste(folder, "plot_Sat_Cassius2_test.pdf", sep="/"),width = 13, height = 14)
# p<-ggarrange(
#   colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius2_pv),
#   heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius2_pv, reject=TRUE),
#   results_pv_all%>%ggplot(aes(x=Sat_test_Cassius2))+
#     geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#     stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(3/n)))+
#     facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# )
# annotate_figure(p,top=text_grob("Saturation Test Cassius 2", face = "bold", size = 14))
# dev.off()
# 
# pdf(paste(folder, "plot_Chi_test.pdf", sep="/"),width = 13, height = 7)
# p<-ggarrange(
#   colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Chi_Test_pv),
#   heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Chi_Test_pv, reject=TRUE),
#   results_pv_all%>%ggplot(aes(x=Chi_test))+
#     geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#     stat_function(fun = dchisq, args = list(df=9))+
#     facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# )
# annotate_figure(p,top=text_grob("Chi Test", face = "bold", size = 14))
# dev.off()


write.csv(results_pv_all, paste(folder, "results_pv_all.csv", sep="/"), row.names = FALSE)


# results_pv_all%>%ggplot(aes(x=Bowker_ts))+
#   geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+  
#   stat_function(fun = dchisq, args = list(df=6))+
#   facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# 
# 
# results_pv_all%>%ggplot(aes(x=Stuart_ts))+
#   geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#   stat_function(fun = dchisq, args = list(df=3))+
#   facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# 
# results_pv_all%>%ggplot(aes(x=IS_ts))+
#   geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#   stat_function(fun = dchisq, args = list(df=3))+
#   facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# 
# results_pv_all%>%ggplot(aes(x=QS_ts))+
#   geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#   stat_function(fun = dnorm, args = list(mean = 0, sd = 1))+
#   facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# 
# results_pv_all%>%ggplot(aes(x=Sat_test_Cassius1))+
#   geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#   stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(3/n)))+
#   facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# 
# results_pv_all%>%ggplot(aes(x=Sat_test_Cassius2))+
#   geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#   stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(3/n)))+
#   facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")
# 
# results_pv_all%>%ggplot(aes(x=Chi_test))+
#   geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7))+
#   stat_function(fun = dchisq, args = list(df=9))+
#   facet_wrap(~Pair, ncol=2) + theme(legend.position = "none")

#test failed/success

Bowker_success<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,"retain","reject"), Bowker_pv)
Stuart_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Stuart_pv)
IS_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), IS_pv)
QS_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), QS_pv)
Reversibility_test<-data.frame(Pair=results_pv_all$Pair,
                               Bowker_success, Stuart_success, IS_success, QS_success)

write.csv(Reversibility_test, paste(folder, "results_rev_test.csv", sep="/"),row.names = FALSE)

Sat_cassius1_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Sat_cassius1_pv)
Sat_cassius2_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Sat_cassius2_pv)
Chi_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Chi_Test_pv)
Saturation_test<-data.frame(Pair=results_pv_all$Pair,
                            Sat_cassius1_success, Sat_cassius2_success, Chi_success)

write.csv(Saturation_test, paste(folder, "results_sat_test.csv", sep="/"),row.names = FALSE)

labels_dataframe<-data.frame(node=d1$node, Bowker_freq=d1$color)
y <- full_join(tree, labels_dataframe, by = 'node')
outputfile<-paste(str_sub(treename, start=0, end=-5), ".tree", sep="")
write.beast(y, "huh")


## Venn diagramm
# create dataframe with number of rejects for each pair

Bowker_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,TRUE,FALSE), Bowker_pv)
Stuart_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,TRUE,FALSE), Stuart_pv)
IS_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,TRUE,FALSE), IS_pv)
QS_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,TRUE,FALSE), QS_pv)

venn_data<-data.frame(Pair=result_tree$Sequences, Bowker=Bowker_rej, Stuart=Stuart_rej, IS=IS_rej, QS=QS_rej)

ggvenn(venn_data, c("Bowker", "Stuart", "QS"))

sink(type="message")
close(zz)