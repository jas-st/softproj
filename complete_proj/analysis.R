zz <- file("log.Rout", open="wt")
sink(zz, type="message")

library(tidyverse)
library(ape)
library(Biostrings)
library(ggtree)
library(ggpubr)
library(scales)



args <- commandArgs()
treename <- args[6]
n<- strtoi(args[7])
raw_results_file<- args[8]
results_file<-args[9]
k<-args[10]


# treename<-"tree.nwk"
# n<-500
# k<-2
# raw_results<-"results_raw_tree.tsv"
# results<-"results_tree.csv"
  
setwd("/Users/luziaberchtold/Desktop/softproj")
result_tree <- read_delim(raw_results_file, 
                          delim = "\t", escape_double = FALSE, 
                          col_types = cols(...13 = col_skip()), 
                          trim_ws = TRUE)

result_tree<-result_tree%>%filter(Pair!="Pair")
result_tree$simulation<-rep(seq(1,k), each=length(unique(result_tree$Pair)))

tree <- read.tree(treename)
p <- ggtree(tree) + geom_tiplab()
seq_labels<-unlist(tree["tip.label"])

colnames_ts<-c("Sat test Cassius 1", "Sat test Cassius 2", "Chi test", "Bowker ts", "Stuart ts", "IS ts", "QS ts")
colnames_nr<-c("Bowker nr", "Stuart nr", "IS nr", "QS nr")

pairs<-strsplit(str_sub(result_tree$Pair, start=2, end=-2), ",")
pairs<-lapply(pairs, strtoi)
pairs<-lapply(pairs, unlist)
pair_1<-unlist(lapply(pairs,function(x) x[1]))
pair_2<-unlist(lapply(pairs,function(x) x[2]))
result_tree$pair_1<-pair_1
result_tree$pair_2<-pair_2

result_tree<-result_tree%>%mutate_at(colnames_nr, ~factor(.x,levels=c("0", "1")))
result_tree[, colnames_ts] <- sapply(result_tree[, colnames_ts], as.numeric)

result_tree_rename<-result_tree%>%dplyr::rename(Sat_test_Cassius1="Sat test Cassius 1", Sat_test_Cassius2="Sat test Cassius 2", 
                                         Chi_test="Chi test", Bowker_ts="Bowker ts", Stuart_ts="Stuart ts", IS_ts="IS ts", QS_ts="QS ts",
                                         Bowker_nr="Bowker nr", Stuart_nr="Stuart nr", IS_nr="IS nr", QS_nr="QS nr")

QS_pv<-2*(1-pnorm(abs(result_tree_rename$QS_ts),mean=0,sd=1)) #standardnormalverteilt
Bowker_pv<-pchisq(q=result_tree_rename$Bowker_ts, df=6, lower.tail=FALSE)#chi square vert mit df=6 #nicht ganz sicher of two sided oder one sided??????
Stuart_pv<-pchisq(q=result_tree_rename$Stuart_ts, df=3, lower.tail=FALSE)#chi square vert mit df=3
IS_pv<-pchisq(q=result_tree_rename$IS_ts, df=3, lower.tail=FALSE)#chi square vert mit df=3
Sat_cassius1_pv<-pnorm(result_tree_rename$Sat_test_Cassius1, 0, sd=sqrt(3/n), lower.tail=FALSE)
Sat_cassius2_pv<-pnorm(result_tree_rename$Sat_test_Cassius2, 0, sd=sqrt(3/n), lower.tail=FALSE)
Chi_Test_pv<-pchisq(q=result_tree_rename$Chi_test, df=9, lower.tail=FALSE)

results_pv_all<-data.frame(Pair=result_tree_rename$Pair, pair_1, pair_2, Bowker_pv, Bowker_ts=result_tree_rename$Bowker_ts,Stuart_pv,Stuart_ts=result_tree_rename$Stuart_ts, 
                                        IS_pv, IS_ts=result_tree_rename$IS_ts, QS_pv, QS_ts=result_tree_rename$QS_ts, 
                                        Sat_cassius1_pv, Sat_test_Cassius1=result_tree_rename$Sat_test_Cassius1, Sat_cassius2_pv, Sat_test_Cassius2=result_tree_rename$Sat_test_Cassius2,
                                        Chi_Test_pv, Chi_test=result_tree_rename$Chi_test, simulation=result_tree_rename$simulation)

results_pv<-results_pv_all %>%group_by(Pair)%>%summarise(across(-simulation,mean, na.rm = TRUE))

heat_success<-function(pair_1, pair_2, test_pv, reject=TRUE, seq_lables=seq_labels){
  pair_1_tmp<-seq_labels[c(pair_1, pair_2)]
  pair_2_tmp<-seq_labels[c(pair_2, pair_1)]
  test_pv<-rep(test_pv,2)
  
  test<-data.frame(pair_1=factor(pair_1_tmp),pair_2=factor(pair_2_tmp),test_pv=test_pv)
  test<-test%>%mutate(test_rej=factor(ifelse(test_pv>=0.05, 0, 1), levels=c(0,1), ordered=TRUE))
  #shape_values<-ifelse(reject==TRUE, c(1,4), c(4,1))
  test%>%ggplot(aes(x=pair_1, y=pair_2))+geom_tile(aes(fill=test_pv))+
    scale_fill_gradient2(low = "white", high = muted("green"), midpoint=0.05)+
    geom_point(aes(shape=test_rej), size=3)+
    labs(x="sequence", y="sequence", fill="p-value", shape="H0")+
    scale_shape_manual(labels=c("keep H0", "reject H0"), values=c(1,4), drop=F)
}

colored_tree <- function(tree, pair_1, pair_2, test_pv){
  paths<-mapply(function(x,y) get.path(tree, x,y), pair_1, pair_2)
  rej_paths<-mapply(function(x) ifelse(x>=0.05,0,1), test_pv )
  
  edges<-sort(unique(unlist(paths)))
  edges_rej<-c(rep(0,length(edges)))
  freq<-table(unlist(paths))
  for(i in edges){
    for(j in 1:length(paths)){
      if(i %in% paths[j][[1]]){
        edges_rej[i]=edges_rej[i]+rej_paths[j]
      }
      edges_rej[i]=edges_rej[i]/freq[i]
    }
  }
  d1<-data.frame(node=edges, color=edges_rej)
  ggtree(tree, layout="slanted") %<+% d1 + aes(color=edges_rej)+ 
    geom_tiplab(offset=0.2, colour="black", size=5) +
    scale_color_gradient2(low = muted("green"), mid="orange", high = muted("red"), midpoint = 0.5, limits=c(0,1), na.value = NA)+
    labs(color="H0 rejected (freq)")+
    geom_label(aes(label=round(edges_rej,4)))
}


pdf("plots/plot_Bowker_test.pdf",width = 13, height = 7)
p<-ggarrange(
  colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv),
  heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv, reject=TRUE)
  )
annotate_figure(p,top=text_grob("Bowker Test", face = "bold", size = 14))
dev.off()

pdf("plots/plot_Stuart_test.pdf",width = 13, height = 7)
p<-ggarrange(
  colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv),
  heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv, reject=TRUE)
)
annotate_figure(p,top=text_grob("Stuart Test", face = "bold", size = 14))
dev.off()

pdf("plots/plot_IS_test.pdf",width = 13, height = 7)
p<-ggarrange(
  colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv),
  heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv, reject=TRUE)
)
annotate_figure(p,top=text_grob("IS Test", face = "bold", size = 14))
dev.off()

pdf("plots/plot_QS_test.pdf",width = 13, height = 7)
p<-ggarrange(
  colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv),
  heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv, reject=TRUE)
)
annotate_figure(p,top=text_grob("QS Test", face = "bold", size = 14))
dev.off()

pdf("plots/plot_Sat_Cassius1_test.pdf",width = 13, height = 7)
p<-ggarrange(
  colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv),
  heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv, reject=TRUE)
)
annotate_figure(p,top=text_grob("Saturation Test Cassius 1", face = "bold", size = 14))
dev.off()

pdf("plots/plot_Sat_Cassius2_test.pdf",width = 13, height = 7)
p<-ggarrange(
  colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius2_pv),
  heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius2_pv, reject=TRUE)
)
annotate_figure(p,top=text_grob("Saturation Test Cassius 2", face = "bold", size = 14))
dev.off()

pdf("plots/plot_Chi_test.pdf",width = 13, height = 7)
p<-ggarrange(
  colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Chi_Test_pv),
  heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Chi_Test_pv, reject=TRUE)
)
annotate_figure(p,top=text_grob("Chi Test", face = "bold", size = 14))
dev.off()


write.csv(results_pv_all, results_file)
sink(type="message")
close(zz)
