#author: Adela Weber
#email: adela[at]nicoweb.com

#Command line : Rscript --vanilla TopicTerm_Comparison.R -i Name.TopicTerm.tab -o SpecificityProbability_plot -t 1000
#This script plots the specificity of a cluster against the others for the events in decreasing probability, making a plot for each cluster.

#!/usr/bin/env Rscript
library("optparse")
library(ggplot2)
library("reshape2")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Topic-term distribution file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="plot", 
              help= "Name of the output [default = %default -> outputs = plot_cl[X].pdf]", metavar="character"),
  make_option(c("-t", "--topEvents"), type="integer", default = -3, 
              help= "Number of most probable events to plot [default = all events]", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input))
{
  print_help(opt_parser)
  stop("At least one input file must be supplied", call.=FALSE)
}



#This function creates a list of names "Cluster 1", ..., "Cluster k" for the k cluster of the topic-term matrix, to label the rows
#ARGS - TopicTerm matrix
Rownames_clusters = function(TopicTerm)
{
  names = c()
  k = dim(TopicTerm)[1] #number of clusters
  for (cl in 1:k)
  {
    names = c(names, paste("Cluster", cl))
  }
  return(names)
}


#This function creates a dataframe of log2(Prob of a cluster/other cluster) for all non repeated pairs of clusters
#ARGS - TopicTerm matrix
Log2_Ratio = function (TopicTerm)
{
  Log2 <- vector()
  Names <- list()
  for (i in 1:(dim(TopicTerm)[1]-1))
  {
    for (j in (i+1):dim(TopicTerm)[1])
    {
      Log2 <- cbind(Log2, log2(TopicTerm[i,]/TopicTerm[j,]))
      Names <- c(Names, paste("log2_ProbCl",i,"_ProbCl",j, sep = ""))
    }
  }
  colnames(Log2) <- Names
  
  return (data.frame(Log2))
}



#This function returns the probability distribution of a cluster in decreasing order.
#ARGS : TopicTerm - Topic term distribution (clusters x events)
#ARGS : cluster_no - index of the cluster for which we want the decreasing order distrib.
Decreasing_order_cluster_distribution = function (TopicTerm, cluster_no)
{
  return (TopicTerm[cluster_no, order(TopicTerm[cluster_no,], decreasing=TRUE)])
}


#This function plots log2 of cl/other cl in regards of the event ranking (decreasing probability of events for this cluster)
#ARGS - cl_no is the index of the cluster we want
#ARGS - log2_df Dataframe of the all different log2 probability ratios
#ARGS - decr is Decreasing probability of events for cluster cl_no
#ARGS - number of top events to plot
Specificity_probability_graph = function (cl_no, log2_df, decr, n, output_name)
{
  #Find the log2 probability ratio of cluster cl_no
  idx_log2 <- grep(paste("Cl", cl_no, sep = ""), colnames(log2_df)) #Indexes of the log2 ratio that interest us
  log2_df_cl <- log2_df[,idx_log2]
  
  #If log2(otherCluster/cl_no) instead of log2(cl_no/otherCluster), invert it
  idx_not_to_invert <- grep(paste("log2_ProbCl", cl_no, sep = ""), colnames(log2_df)[idx_log2])
  if (length(idx_not_to_invert) == 0) #All log2 to invert
  {
    log2_df_cl <- -log2_df_cl
    if (dim(log2_df)[2] == 1) #If only 2 clusters in total
    {
      colnames(log2_df) <- paste("-", colnames(log2_df), sep="")
    } else #More than 2 clusters
    {
      colnames(log2_df_cl) <- paste("-", colnames(log2_df_cl), sep="")
    }
  }
  else if (length(idx_not_to_invert) < length(colnames(log2_df_cl))) #Some log2 to invert
  {
    log2_df_cl[, -idx_not_to_invert] <- -log2_df_cl[, -idx_not_to_invert]
    colnames(log2_df_cl)[-idx_not_to_invert] <- paste("-", colnames(log2_df_cl)[-idx_not_to_invert], sep="")
  }
  
  #Events in decreasing probability order for cluster no_cl
  Order_events <- names(decr)
  if (dim(log2_df)[2] == 1) #If only 2 clusters
  {
    Data <- data.frame(cbind(rep(colnames(log2_df), n), log2_df[Order_events[1:n],]))
    Data[,2] <- as.numeric(levels(Data[,2]))[Data$Log2]
    colnames(Data) <- c("Clusters", "Log2")
    EventsRanking <- 1:n
  } else #more than 2 clusters
  {
    Data <- melt(log2_df_cl[Order_events[1:n],], variable.name="Clusters", value.name = "Log2") #Take only n top events
    EventsRanking <- rep(1:n, (dim(log2_df_cl)[2]))
  }
  
  #Graphical representation
  ggplot(Data, aes(x= EventsRanking, y=Data$Log2, colour = Data$Clusters)) +
    geom_point(alpha = 0.3) + #COMMENT OUT this line if you don't want the dots
    geom_smooth() +
    xlab(paste("Probability ranking of", n, "top events"))+
    ylab(paste("log2(probabilityCl", cl_no,"/probabilityOtherClusters", sep=""))+
    ggtitle(paste("Specificity of events against\nevent probability ranking for cluster", cl_no))
  ggsave(paste(output_name, "_cl", cl_no, ".pdf", sep = ""), device = "pdf")
}



#####EXTRA FUNCTIONS #######

#This function returns the event names of the n most probable events of a cluster, and the probabilities if wanted.
#ARGS : n - number of events wanted
#ARGS : Decr_cl - probability distribution of a cluster in decreasing order, with names of the events as tags
#ARGS : probabilities - logical, if TRUE gives the probabilies (names as tags), if FALSE gives only the names of the events
n_top_events = function(n, Decr_cl, probabilities)
{
  if (probabilities)
  {
    return(Decr_cl[1:n])
  } else
  {
    return(names(Decr_cl[1:n])) 
  }
}



#This function finds the number of top events of a cluster - they are the number of events needed to have more than x% of your distribution
#ARGS : Decr_proba - decreasing topic term probability of your cluster
#ARGS : x - minimum cumulated probability sum we want to overcome
Find_n_best_events = function (Decr_proba, x)
{
  n = 1
  name <- names(Decr_proba[n])
  while(cumsum(Decr_proba[1:n])[n] < x)
  {
    n <- n+1
    name <- names(Decr_proba[n])
  }
  return (n)
}


#This function returns a list of the skipped exons sizes for all events - to analyze microexon presence
#ARGS : events - ordered (prob ranking) event names
#Reminder of events_ids : gene_id;type_of_event_code;chromosome:end_exon1-start_exon2(skipped_exon):end_exon2-start_exon3:strand(+or-)
ExonSize = function(events)
{
  Size = list()
  for (ev in events)
  {
    Exonlims <- strsplit(strsplit(ev, "-")[[1]][2], ":")[[1]]
    Size = c(Size, as.integer(Exonlims[2])-as.integer(Exonlims[1]))
  }
  return(unlist(Size))
}





##### MAIN #####


#TopicTerm matrices
TopicTerm <- as.matrix(read.table(opt$input, header=TRUE, check.names=FALSE, sep = "\t", as.is=TRUE))
rownames(TopicTerm) = Rownames_clusters(TopicTerm)

if (opt$topEvents > dim(TopicTerm)[2] || opt$topEvents <= 0)
{
  print("Taking default value for number of most probable events to plot")
  opt$topEvents <- dim(TopicTerm)[2]
}

#Calculating the log2 cluster probability ratios
Log2_df <- Log2_Ratio(TopicTerm)

#Event-clusters distribution, in decreasing order
for (i in 1:dim(TopicTerm)[1])
{
  assign(paste("Decr_cl", i, sep = ""), Decreasing_order_cluster_distribution(TopicTerm, i))
}

#Specificity probability graph
for (i in 1:dim(TopicTerm)[1])
{
  Specificity_probability_graph(i, Log2_df, eval(as.name(paste("Decr_cl", i, sep=""))), opt$topEvents, opt$output)
}
