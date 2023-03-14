options(warn=-1)

suppressPackageStartupMessages({
  if(!require("tidyverse")) install.packages("tidyverse", repos = "http://cran.us.r-project.org"); library(tidyverse)
  if(!require("ggpubr")) install.packages("ggpubr", repos = "http://cran.us.r-project.org"); library(ggpubr)
  if(!require("ggvenn")) install.packages("ggvenn", repos = "http://cran.us.r-project.org"); library(ggvenn)
  if(!require("tidytree")) install.packages("tidytree", repos = "http://cran.us.r-project.org"); library(tidytree)
  })


suppressWarnings(suppressMessages(
if(!require("BiocManager", quietly = TRUE))
{
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
    BiocManager::install(c("ggtree","treeio"))
}))

suppressWarnings(suppressMessages(library("ggtree")))
suppressWarnings(suppressMessages(library("treeio")))


#read in command line arguments
args <- commandArgs()
treename <- args[6]
n<- strtoi(args[7])
raw_results_file<- args[8]
results_file<-args[9]
k<-args[10]
s<-args[11]
folder<-strsplit(raw_results_file, "/")[[1]][1]


#read in raw results file
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

#remove column names in between simulations + add column with simulation number
result_tree_rename<-results_raw_tree%>%dplyr::rename(Sequences="Sequences")
result_tree<-result_tree_rename%>%filter(Sequences!="Sequences")
result_tree$simulation<-rep(seq(1,k), each=length(unique(result_tree$Sequences)))

#read in tree file 
tree <- read.tree(treename)
seq_labels<-unlist(tree["tip.label"])
#read sequence labels and add as two columns in results dataframe 
pairs<-strsplit(str_sub(result_tree$Sequences, start=2, end=-2), ";")
pairs<-lapply(pairs, unlist)
pair_1<-unlist(lapply(pairs,function(x) x[1]))
pair_2<-unlist(lapply(pairs,function(x) x[2]))
pair_1<-gsub(" ", "", pair_1, fixed = TRUE)
pair_2<-gsub(" ", "", pair_2, fixed = TRUE)
result_tree$seq_1<-pair_1
result_tree$seq_2<-pair_2

#calculate pvalues
QS_pv<-2*(1-pnorm(abs(result_tree$Proposed_test),mean=0,sd=1)) #standardnormalverteilt
Bowker_pv<-pchisq(q=result_tree$Bowker_test, df=6, lower.tail=FALSE)#chi square vert mit df=6 #nicht ganz sicher of two sided oder one sided??????
Stuart_pv<-pchisq(q=result_tree$Stuart_test, df=3, lower.tail=FALSE)#chi square vert mit df=3
IS_pv<-pchisq(q=result_tree$Internal_Symmetry, df=3, lower.tail=FALSE)#chi square vert mit df=3

if (s=="true") Sat_cassius1_pv<-pnorm(result_tree$Sat_test_Cassius1, 0, sd=sqrt(3/result_tree$Alignment_Length), lower.tail=FALSE)
if (s=="true") Sat_cassius2_pv<-pnorm(result_tree$Sat_test_Cassius2, 0, sd=sqrt(3/result_tree$Alignment_Length), lower.tail=FALSE)
if (s=="true") Chi_Test_pv<-pchisq(q=result_tree$Chi_test, df=9, lower.tail=FALSE)

#dataframe with test statistic and pvalues
if (s=="true") {
  results_pv_all <- data.frame(Pair=result_tree$Sequences, pair_1, pair_2, Bowker_pv, Bowker_ts=result_tree$Bowker_test,Stuart_pv,Stuart_ts=result_tree$Stuart_test, 
             IS_pv, IS_ts=result_tree$Internal_Symmetry, QS_pv, QS_ts=result_tree$Proposed_test, 
             Sat_cassius1_pv, Sat_test_Cassius1=result_tree$Sat_test_Cassius1, Sat_cassius2_pv, Sat_test_Cassius2=result_tree$Sat_test_Cassius2,
             Chi_Test_pv, Chi_test=result_tree$Chi_test, simulation=result_tree$simulation)
} else {
  results_pv_all <- data.frame(Pair=result_tree$Sequences, pair_1, pair_2, Bowker_pv, Bowker_ts=result_tree$Bowker_test,Stuart_pv,Stuart_ts=result_tree$Stuart_test, 
             IS_pv, IS_ts=result_tree$Internal_Symmetry, QS_pv, QS_ts=result_tree$Proposed_test, 
             simulation=result_tree$simulation)
}

#take mean over multiple simulations
results_pv<-results_pv_all %>%group_by(Pair)%>%summarise(across(-c(simulation, pair_1, pair_2),mean, na.rm = TRUE))
results_pv$pair_1<-results_pv_all$pair_1[c(seq(1,length(results_pv$Pair)))]
results_pv$pair_2<-results_pv_all$pair_2[c(seq(1,length(results_pv$Pair)))]


#function for heatmap
heat_success<-function(seq_pair, test_pv){
  pairs<-strsplit(str_sub(seq_pair, start=2, end=-2), ";")
  pairs<-lapply(pairs, unlist)
  pair_1<-unlist(lapply(pairs,function(x) x[1]))
  pair_2<-unlist(lapply(pairs,function(x) x[2]))
  pair_1<-gsub(" ", "", pair_1, fixed = TRUE)
  pair_2<-gsub(" ", "", pair_2, fixed = TRUE)
  pair_1_tmp<-c(pair_1, pair_2)
  pair_2_tmp<-c(pair_2, pair_1)
  test_pv<-rep(test_pv,2)
  
  test<-data.frame(pair_1=factor(pair_1_tmp),pair_2=factor(pair_2_tmp),test_pv=test_pv)
  
  test<-test%>%mutate(test_rej=factor(ifelse(test_pv>=0.05, 0, 1), levels=c(0,1), ordered=TRUE))
  #shape_values<-ifelse(reject==TRUE, c(1,4), c(4,1))
  test%>%ggplot(aes(x=pair_1, y=pair_2))+geom_tile(aes(fill=test_pv))+
    scale_fill_gradient(high="red", low="blue", limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE)) +
    geom_point(aes(shape=test_rej), size=3)+
    labs(title="mean p value of each pair" ,x="",y="sequence", fill="p-value", shape="H0")+
    scale_shape_manual(labels=c("retain", "reject"), values=c(1,4), drop=F)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank() 
    )
}

#function that maps results on tree
edges_rejected_freq <- function(tree, seq_pair, test_pv){
  pairs<-strsplit(str_sub(seq_pair, start=2, end=-2), ";")
  pairs<-lapply(pairs, unlist)
  pair_1<-unlist(lapply(pairs,function(x) x[1]))
  pair_2<-unlist(lapply(pairs,function(x) x[2]))
  pair_1<-gsub(" ", "", pair_1, fixed = TRUE)
  pair_2<-gsub(" ", "", pair_2, fixed = TRUE)
  pair_1_index<-mapply(function(x) which(x==seq_labels), pair_1)
  pair_2_index<-mapply(function(x) which(x==seq_labels), pair_2)
  paths<-mapply(function(x,y) ggtree::get.path(tree, x,y), pair_1_index, pair_2_index)
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
      if(p %in% paths[[j]] && n %in% paths[[j]] && is.na(rej_paths[j])==FALSE){
        edges_rej[i]<-edges_rej[i]+rej_paths[j]
        edges_freq[i]<-edges_freq[i]+1
      }
    }
    if(edges_freq[i]!=0 && edges_rej[i]!=0) edges_rej[i]=edges_rej[i]/edges_freq[i]
  }
  edge["edge_rej"]<-round(edges_rej,2)
  return(edge)
}
  
  
colored_tree <- function(tree, seq_pair, test_pv){
  d1<-edges_rejected_freq(tree, seq_pair, test_pv)
  p<-ggtree(tree, layout="slanted") + ggtree::geom_tiplab(colour="black", size=5)+
    labs(title="H0 rejected (freq) on each edge")
  p %<+% d1 + geom_label(aes(x=branch, label=edge_rej)) + aes(color=edge_rej)+ 
    scale_colour_gradient2(low="blue", mid="#FF0099", high="red",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    labs(color="H0 rejected (freq)")
}

## write result csv tables 
write.csv(results_pv_all, results_file, row.names = FALSE)

#test retain/reject
Bowker_success<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,"retain","reject"), Bowker_pv)
Stuart_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Stuart_pv)
IS_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), IS_pv)
QS_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), QS_pv)
Reversibility_test<-data.frame(Pair=results_pv_all$Pair, Simulation=results_pv_all$simulation,
                               Bowker_success, Stuart_success, IS_success, QS_success)

write.csv(Reversibility_test, paste(folder, "results_rev_test.csv", sep="/"),row.names = FALSE)

#write saturation results table
if (s=="true") {
  Sat_cassius1_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Sat_cassius1_pv)
  Sat_cassius2_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Sat_cassius2_pv)
  Chi_success<-mapply(function(x) ifelse(x>=0.05,"retain","reject"), Chi_Test_pv)
  Saturation_test<-data.frame(Pair=results_pv_all$Pair, Simulation=results_pv_all$simulation,
                              Sat_cassius1_success, Sat_cassius2_success, Chi_success)
  write.csv(Saturation_test, paste(folder, "results_sat_test.csv", sep="/"),row.names = FALSE)
}

## save trees in file with labels for H0 rejection of every test
Bowker_edge_freq<-edges_rejected_freq(tree, results_pv$Pair, results_pv$Bowker_pv)
Stuart_edge_freq<-edges_rejected_freq(tree, results_pv$Pair, results_pv$Stuart_pv)
IS_edge_freq<-edges_rejected_freq(tree, results_pv$Pair, results_pv$IS_pv)
QS_edge_freq<-edges_rejected_freq(tree, results_pv$Pair, results_pv$QS_pv)
if (s=="true") Sat1_edge_freq<-edges_rejected_freq(tree, results_pv$Pair, results_pv$Sat_cassius1_pv)
if (s=="true") Sat2_edge_freq<-edges_rejected_freq(tree, results_pv$Pair, results_pv$Sat_cassius2_pv)
if (s=="true") Chi_edge_freq<-edges_rejected_freq(tree, results_pv$Pair, results_pv$Chi_Test_pv)

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

outputfile<-paste(str_sub(treename, start=0, end=-5), ".tree", sep="")
write.beast(as.treedata(y), file=paste(folder, outputfile, sep="/"))


## all plots
#Bowker test
pdf(paste(folder, "plot_Bowker_test.pdf", sep="/"),width = 13, height = 14)
p<-ggarrange(
  ggarrange(colored_tree(tree, results_pv$Pair, results_pv$Bowker_pv),
  heat_success(results_pv$Pair, results_pv$Bowker_pv), nrow=1),
  results_pv_all%>%ggplot(aes(x=Bowker_ts))+
    geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7), bins=30)+  
    stat_function(fun = dchisq, args = list(df=6))+
    facet_wrap(~Pair) + theme(legend.position = "none")+
    labs(title="distribution of test statistics"), nrow=2
  )
annotate_figure(p,top=text_grob("Bowker Test", face = "bold", size = 14))
dev.off()

#Stuart test
pdf(paste(folder, "plot_Stuart_test.pdf", sep="/"),width = 13, height = 14)
p<-ggarrange(
  ggarrange(colored_tree(tree, results_pv$Pair, results_pv$Stuart_pv),
  heat_success(results_pv$Pair, results_pv$Stuart_pv), nrow=1),
  results_pv_all%>%ggplot(aes(x=Stuart_ts))+
    geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7), bins=30)+
    stat_function(fun = dchisq, args = list(df=3))+
    facet_wrap(~Pair) + theme(legend.position = "none")+
    labs(title="distribution of test statistics"), nrow=2
)
annotate_figure(p,top=text_grob("Stuart Test", face = "bold", size = 14))
dev.off()

#IS test
pdf(paste(folder, "plot_IS_test.pdf", sep="/"),width = 13, height = 14)
p<-ggarrange(
  ggarrange(colored_tree(tree, results_pv$Pair, results_pv$IS_pv),
  heat_success(results_pv$Pair, results_pv$IS_pv), nrow=1),
  results_pv_all%>%ggplot(aes(x=IS_ts))+
    geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7), bins=30)+
    stat_function(fun = dchisq, args = list(df=3))+
    facet_wrap(~Pair) + theme(legend.position = "none")+
    labs(title="distribution of test statistics"), nrow=2
)
annotate_figure(p,top=text_grob("IS Test", face = "bold", size = 14))
dev.off()

#QS test
pdf(paste(folder, "plot_QS_test.pdf", sep="/"),width = 13, height = 14)
p<-ggarrange(
  ggarrange(colored_tree(tree, results_pv$Pair, results_pv$QS_pv),
  heat_success(results_pv$Pair, results_pv$QS_pv), nrow=1),
  results_pv_all%>%ggplot(aes(x=QS_ts))+
    geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7), bins=30)+
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1))+
    facet_wrap(~Pair) + theme(legend.position = "none")+
    labs(title="distribution of test statistics"), nrow=2
)
annotate_figure(p,top=text_grob("QS Test", face = "bold", size = 14))
dev.off()

#Cassius 1
if (s=="true"){
  filename<-paste(folder, "plot_Sat_Cassius1_test.pdf", sep="/")
  p<-ggarrange(
    ggarrange(colored_tree(tree, results_pv$Pair, results_pv$Sat_cassius1_pv),
    heat_success(results_pv$Pair, results_pv$Sat_cassius1_pv), nrow=1),
    results_pv_all%>%ggplot(aes(x=Sat_test_Cassius1))+
      geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7), bins=30)+
      stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(3/n)))+
      facet_wrap(~Pair) + theme(legend.position = "none")+
      labs(title="distribution of test statistics"), nrow=2
  )
  annotate_figure(p,top=text_grob("Saturation Test Cassius 1", face = "bold", size = 14))
  ggsave(filename,width = 13, height = 14)
}

#Cassius 2
if (s=="true"){
  filename<-paste(folder, "plot_Sat_Cassius2_test.pdf", sep="/")
  p<-ggarrange(
    ggarrange(colored_tree(tree, results_pv$Pair, results_pv$Sat_cassius2_pv),
    heat_success(results_pv$Pair, results_pv$Sat_cassius2_pv), nrow=1),
    results_pv_all%>%ggplot(aes(x=Sat_test_Cassius2))+
      geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7), bins=30)+
      stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(3/n)))+
      facet_wrap(~Pair) + theme(legend.position = "none")+
      labs(title="distribution of test statistics"), nrow=2
  )
  annotate_figure(p,top=text_grob("Saturation Test Cassius 2", face = "bold", size = 14))
  ggsave(filename,width = 13, height = 14)
}

#Chi test
if (s=="true"){
  filename<-paste(folder, "plot_Chi_test.pdf", sep="/")
  p<-ggarrange(
    ggarrange(colored_tree(tree, results_pv$Pair, results_pv$Chi_Test_pv),
    heat_success(results_pv$Pair, results_pv$Chi_Test_pv), nrow=1),
    results_pv_all%>%ggplot(aes(x=Chi_test))+
      geom_histogram(aes(y = ..density.., colour=Pair, alpha=0.7), bins=30)+
      stat_function(fun = dchisq, args = list(df=9))+
      facet_wrap(~Pair, nrow=3) + theme(legend.position = "none")+
      labs(title="distribution of test statistics"), nrow=2
  )
  annotate_figure(p,top=text_grob("Chi Test", face = "bold", size = 14))
  ggsave(filename,width = 13, height = 14)
}

## Venn diagramm
# create dataframe with number of rejects for each pair
Bowker_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,FALSE,TRUE), Bowker_pv)
Stuart_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,FALSE,TRUE), Stuart_pv)
IS_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,FALSE,TRUE), IS_pv)
QS_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,FALSE,TRUE), QS_pv)

venn_data<-data.frame(Pair=result_tree$Sequences, Bowker=Bowker_rej, Stuart=Stuart_rej, IS=IS_rej, QS=QS_rej)

for (pair in unique(venn_data$Pair)){
  venn_pair<-venn_data%>%filter(Pair==pair)
  plot(ggvenn(venn_pair, c("Bowker", "Stuart", "QS"), show_percentage = FALSE)+labs(title=pair) +
         theme(plot.title = element_text(hjust = 0.5)))
}
