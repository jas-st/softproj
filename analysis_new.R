zz <- file("log.Rout", open="wt")
sink(zz, type="message")


library(scales)
library(tidyverse)
library(ape)
library(Biostrings)
library(ggtree)
library(ggpubr)
library(phangorn)

args <- commandArgs()
treename <- args[6]
n<- strtoi(args[7])
raw_results_file<- args[8]
results_file<-args[9]
#k<-args[10]

folder<-strsplit(results_file, "/")[[1]][1]

# treename<-"tree.nwk"
# n<-10000
# k<-10
# raw_results_file<-"results_tree/results_raw_tree.csv"
# results<-"results_tree.csv"

raw_results_file <- "datasets_eval/worobey_2014a_corrected.csv"
raw_results_file <- "datasets_eval/unmack_2013_new.csv"
treename <- "datasets_eval/worobey_2014a/alignment.nex.treefile"
treename <- "datasets_eval/unmack_2013_new/alignment.nex.treefile"
folder <- "analysis_worobey"
folder <- "datasets_eval/unmack_2013_new"

results_raw_tree_unclean <- read_csv(raw_results_file, 
                                     col_types = cols(Sat_test_Cassius1 = col_number(), 
                                                      Sat_test_Cassius2 = col_number(), 
                                                      Chi_test = col_number(), Bowker_test = col_number(), 
                                                      Stuart_test = col_number(), Internal_Symmetry = col_number(), 
                                                      Proposed_test = col_number(), Alignment_Length = col_number()))
results_raw_tree <- na.omit(results_raw_tree_unclean)

result_tree_rename<-results_raw_tree%>%dplyr::rename(Sequences="Sequences compared")


result_tree<-result_tree_rename%>%filter(Sequences!="Sequences compared")
#result_tree$simulation<-rep(seq(1,k), each=length(unique(result_tree$Sequences)))

tree <- read.tree(treename)
tree <- midpoint(tree)
p <- ggtree(tree) + geom_tiplab()
seq_labels<-unlist(tree["tip.label"])

colnames_ts<-c("Sat test Cassius 1", "Sat test Cassius 2", "Chi test", "Bowker ts", "Stuart ts", "IS ts", "QS ts")

pairs<-strsplit(str_sub(result_tree$Sequences, start=2, end=-2), ";")
pairs<-lapply(pairs, unlist)
pair_1<-unlist(lapply(pairs,function(x) x[1]))
pair_2<-unlist(lapply(pairs,function(x) x[2]))
pair_1<-gsub(" ", "", pair_1, fixed = TRUE)
pair_2<-gsub(" ", "", pair_2, fixed = TRUE)
result_tree$seq_1<-pair_1
result_tree$seq_2<-pair_2

QS_pv<-2*(1-pnorm(abs(result_tree$Proposed_test),mean=0,sd=1)) #standardnormalverteilt
Bowker_pv<-pchisq(q=result_tree$Bowker_test, df=6, lower.tail=FALSE)#chi square vert mit df=6 #nicht ganz sicher of two sided oder one sided??????
Stuart_pv<-pchisq(q=result_tree$Stuart_test, df=3, lower.tail=FALSE)#chi square vert mit df=3
IS_pv<-pchisq(q=result_tree$Internal_Symmetry, df=3, lower.tail=FALSE)#chi square vert mit df=3
Sat_cassius1_pv<-pnorm(result_tree$Sat_test_Cassius1, 0, sd=sqrt(3/result_tree$Alignment_Length), lower.tail=FALSE)
Sat_cassius2_pv<-pnorm(result_tree$Sat_test_Cassius2, 0, sd=sqrt(3/result_tree$Alignment_Length), lower.tail=FALSE)
Chi_Test_pv<-pchisq(q=result_tree$Chi_test, df=9, lower.tail=FALSE)

results_pv_all<-data.frame(Pair=result_tree$Sequences, pair_1, pair_2, Bowker_pv, Bowker_ts=result_tree$Bowker_test,Stuart_pv,Stuart_ts=result_tree$Stuart_test, 
                                        IS_pv, IS_ts=result_tree$Internal_Symmetry, QS_pv, QS_ts=result_tree$Proposed_test, 
                                        Sat_cassius1_pv, Sat_test_Cassius1=result_tree$Sat_test_Cassius1, Sat_cassius2_pv, Sat_test_Cassius2=result_tree$Sat_test_Cassius2,
                                        Chi_Test_pv, Chi_test=result_tree$Chi_test)#, simulation=result_tree$simulation)

results_pv<-results_pv_all #%>%group_by(Pair)%>%summarise(across(-c(simulation, pair_1, pair_2),mean, na.rm = TRUE))
results_pv$pair_1<-results_pv_all$pair_1[c(seq(1,length(results_pv$Pair)))]
results_pv$pair_2<-results_pv_all$pair_2[c(seq(1,length(results_pv$Pair)))]


heat_success<-function(pair_1, pair_2, test_pv, reject=TRUE, seq_lables=seq_labels){
  pair_1_tmp<-c(pair_1, pair_2)
  pair_2_tmp<-c(pair_2, pair_1)
  test_pv<-rep(test_pv,2)
  
  test<-data.frame(pair_1=factor(pair_1_tmp),pair_2=factor(pair_2_tmp),test_pv=test_pv)
  
  test<-test%>%mutate(test_rej=factor(ifelse(test_pv>=0.05, 0, 1), levels=c(0,1), ordered=TRUE))
  #shape_values<-ifelse(reject==TRUE, c(1,4), c(4,1))
  test%>%ggplot(aes(x=pair_1, y=pair_2))+geom_tile(aes(fill=test_pv))+
    #geom_raster(aes(fill = test_pv), interpolate=TRUE) +
    scale_fill_gradient(high="blue",low="red", limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE)) +
    #geom_point(aes(shape=test_rej), size=3)+
    labs(x="sequence 1", y="sequence 2", fill="p-value", shape="H0")+
    #scale_shape_manual(labels=c("keep H0", "reject H0"), values=c(1,4), drop=F) +
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(),axis.ticks.y=element_blank())
}

ggtree(tree) + geom_text(aes(label=node), hjust=-.3)

colored_tree <- function(tree, pair_1, pair_2, test_pv){
  pair_1_index<-mapply(function(x) which(x==seq_labels), pair_1)
  pair_2_index<-mapply(function(x) which(x==seq_labels), pair_2)
  paths<-mapply(function(x,y) get.path(tree, x,y), pair_1_index, pair_2_index)
  rej_paths<-mapply(function(x) ifelse(x>=0.05,0,1), test_pv )

  edges<-sort(unique(unlist(paths)))
  edges_rej<-c(rep(0,length(edges)))
  freq<-table(unlist(paths))
  for(i in edges){
    for(j in 1:length(paths)){
      if(i %in% paths[j][[1]]){
        edges_rej[i]=edges_rej[i]+rej_paths[j]
      }
    }
    edges_rej[i]=edges_rej[i]/freq[i]
  }
  #print(edges_rej), layout="circular"

  d1<-data.frame(node=edges, color=edges_rej)
  ggtree(tree) %<+% d1 + aes(color=edges_rej)+ 
    #geom_tiplab(colour="black", size=3) +
    scale_colour_gradient2(low="blue", mid="#FF0099", high="red",midpoint=0.5, limits=c(0,1), guide = guide_colorbar(draw.ulim = TRUE, draw.llim = TRUE))+
    labs(color="H0 rejected (freq)")#+
    #geom_label(aes(x=branch, label=round(edges_rej,4))) 
}
q <- ggtree(tree) + geom_tiplab()

Nnumber <- as.numeric(tree['Nnode'])+1
tree_root <- getRoot(tree)
collapse_list <- list()
min_length <- 100
for (node in (Nnumber+1):(2*Nnumber-1)) {
  path_to_root_length <- length(get.path(tree,tree_root,node))
  children <- Children(tree, node)
  for (child in children) {
    if (d1[child,2]<0.3 & path_to_root_length > 4 & child > 146) {
      if (path_to_root_length <= min_length) {
        min_length <- path_to_root_length
        p <- p %>% collapse(node = child)
        collapse_list<-append(collapse_list, child)
      }
      
    }
  }
}
for (child in unlist(collapse_list)) {
  p<-p + geom_point2(aes(subset=(node==147)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==275)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==153)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==269)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==154)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==254)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==155)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==157)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==156)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==180)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==179)), shape=21, size = 3, fill='green')
  p<-p + geom_point2(aes(subset=(node==218)), shape=21, size = 3, fill='green')
}



pdf(paste(folder, "plot_Bowker_test.pdf", sep="/"),width = 30, height = 30)

  p <- colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv)
  heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Bowker_pv, reject=TRUE)
  
#annotate_figure(p,top=text_grob("Bowker Test", face = "bold", size = 14))
dev.off()

pdf(paste(folder, "plot_Stuart_test.pdf", sep="/"),width = 13, height = 14)

  p<-colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv)
  heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Stuart_pv, reject=TRUE)
  results_pv_all%>%ggplot(aes(x=Stuart_ts))+
    geom_histogram(bins = 40, aes(y = after_stat(!!str2lang("density"))))+
    stat_function(fun = dchisq, args = list(df=3))

#annotate_figure(p,top=text_grob("Stuart Test", face = "bold", size = 14))
dev.off()

# pdf(paste(folder, "plot_IS_test.pdf", sep="/"),width = 13, height = 14)
# p<-ggarrange(
#   p<-colored_tree(tree, results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv),
#   heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$IS_pv, reject=TRUE),
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
#   heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$QS_pv, reject=TRUE),
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
#   heat_success(results_pv$pair_1, results_pv$pair_2, results_pv$Sat_cassius1_pv, reject=TRUE),
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
library(ggvenn)
Bowker_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,TRUE,FALSE), Bowker_pv)
Stuart_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,TRUE,FALSE), Stuart_pv)
IS_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,TRUE,FALSE), IS_pv)
QS_rej<-mapply(function(x) ifelse(!is.na(x) & x>=0.05,TRUE,FALSE), QS_pv)

venn_data<-data.frame(Pair=result_tree$Sequences, Bowker=Bowker_rej, Stuart=Stuart_rej, IS=IS_rej, QS=QS_rej)

for (pair in unique(venn_data$Pair)){
  venn_pair<-venn_data%>%filter(Pair==pair)
  plot(ggvenn(venn_pair, c("Bowker", "Stuart", "QS"), show_percentage = FALSE)+labs(title=pair) +
         theme(plot.title = element_text(hjust = 0.5)))
}
library(ggVennDiagram)
ggvenn(venn_data, c("Bowker", "Stuart", "QS"),fill_color = c("#EAF1F3", "#FFACB8", "red"))
ggVennDiagram(list(c(0,0,0,0,0),c(1,0,1,1,1)), color="black")
ggVennDiagram(lapply(c(venn_data[2],venn_data[3],venn_data[5]), function(x) which(x == FALSE))) + 
  ggplot2::scale_fill_gradient(low="white",high = "red") +
  scale_color_manual(values = c("black","black","black"))

sink(type="message")
close(zz)