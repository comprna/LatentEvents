#author: Adela Weber
#email: adela[at]nicoweb.com

#Command line : Rscript --vanilla Plot_DocumentTopic.R -i Name.DocTopic.tab -o DocTopic_plot.pdf -f Sample_label.txt
#This script makes the structure histogram plot of the document-topic distribution (sample-cluster distribution), with samples grouped by labels.

#CountClust : Copyright (c) Dey K, Hsiao J and Stephens M (2016)

#!/usr/bin/env Rscript
library("optparse")
library(CountClust)
library(RColorBrewer)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Document-topic distribution file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="plot.pdf", 
              help= "Name of the output plot [default= %default]", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help= "Tab-separated sample-label file with headers, containing (at least) all samples (documents) from the document-topic distribution file.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input) || is.null(opt$file))
{
  print_help(opt_parser)
  stop("One input file and one sample-label file arguments must be supplied", call.=FALSE)
}


#DocTopic matrix extraction
docTopic <- as.matrix(read.table(opt$input, header=FALSE, row.names = 1, sep = "\t", as.is=TRUE))
n_topics <- dim(docTopic)[2] #number of clusters (topics)
colnames(docTopic) <- seq(1, n_topics)


Sample_Label <- read.table(opt$file, header=TRUE, sep = "\t", check.names=FALSE, as.is=TRUE)
Labels <- Sample_Label[,2] #Label names
names(Labels) <- Sample_Label[,1] #Sample names
Labels <- Labels[rownames(docTopic)] #Select and order only the samples from docTopic
y_label <- colnames(Sample_Label)[2]

if (length(na.omit(Labels)) < length(rownames(docTopic)))
{
  stop("Mismatch in names of samples", call.=FALSE)
}

#Annotation data frame
Annotation <- data.frame(sample_id = rownames(docTopic), tissue_label = Labels)

#Structure plot
pdf(NULL)
StructureGGplot(omega = docTopic, 
                annotation = Annotation, 
                palette = colorRampPalette(brewer.pal(12,"Paired"))(n_topics), 
                yaxis_label = y_label, 
                order_sample = TRUE, 
                axis_tick = list(axis_ticks_length = .1, axis_ticks_lwd_y = .1, axis_ticks_lwd_x = .1,  axis_label_size = 7, axis_label_face = "bold")
                )
ggsave(filename = opt$output, device = "pdf")
