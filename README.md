# LatentEvents
Latent Dirichlet Allocation on alternative splicing events.

Here we explain how these tools work.


## Overview
----------------------------

Latent Dirichlet Allocation is a generative statistical model (belonging to probabilistic topic models) usually used on text corpora to find a series of topics. We will be applying this principle to find clusters for a series of samples, mostly from biological data. A great aspect of LDA is that it can allow each sample to belong to multiple clusters, with a probability associated to each of the clusters.

In our code we are using the scikit-learn's implementation of Latent Dirichlet Allocation. From a document-term counts matrix and a series of parameters (number of topics, etc), this function can give us a document-topic distribution and a topic-term distribution (Fig.1).

To adapt LDA to our needs :
- **Counts** have become **PSI measures**
- **Documents** have become **samples**
- **Topics** have become **clusters**
- **Terms** have become **alternative splicing events**

This vocabulary has been used interchangeably throughout the code.


![Figure1.jpg](https://github.com/comprna/LatentEvents/blob/master/Figures/Matrices.jpg)
**Fig.1** - Latent Dirichlet Allocation (LDA) input and output matrices.

Figure 2 shows with more clarity the relations between the distributions and the different words of vocabulary. Each sample will have a cluster distribution and each cluster will have an event distribution.

![Figure2.jpg](https://github.com/comprna/LatentEvents/blob/master/Figures/Distributions.jpg)
**Fig.2** - Intuitions on Latent Dirichlet Allocation. In this imaginary example we have two samples, three clusters and seven events, the possible distributions we could obtain, and how they all connect together.

These tools are ready to use, they can be used directly after download.

Greatly inspired by : Dey KK, Hsiao CJ, Stephens M (2017) Visualizing the structure of RNA-seq expression data using grade of membership models. PLoS Genet 13(3): e1006599. https://doi.org/10.1371/journal.pgen.1006599

## Latent Dirichlet Allocation on splicing events - LatentEvents_LDA.py
----------------------------

Developed on Python 3.4

### Modules used
- scikit-learn version 0.18.2 (0.18 minimum required)
- NumPy version 1.12.1
- needs SUPPA's tools.py (directory lib/)

### Description 
This tool calculates PSIs for a given ioe and transcript expression files and clusters the filtered alternative splicing events with a Latent Dirichlet Allocation model. Events can be filtered for a given event type, to keep micro exons only and/or with a given minimal expression.

### Command
```
python LatentEvents_LDA.py [options]
```
**Command (minimal options)**
```
python LatentEvents_LDA.py -i <ioe-file> -e <expression-file> -o <output-name-for-all-outputs> -l <number-clusters> <number-LDA-max-iterations> <number-feature-selection-iterations> <number-events-to-check-for-feature-selection> -t <number-top-events-in-output> -r <number-random-matrices>
```

### Command options

- **-h** | **--help**: Display the help message describing the different parameters.
- **-i** | **--ioe-file**: Input file with the definition of the events and corresponding transcripts.
- **-e** | **--expression-file**: Transcript expression file containing the abundances of all transcripts.
- **-o** | **--output-file**: Common name for all outputs (later distinguished by extensions: **&lt;output-file&gt;**.DocTopic.tab, **&lt;output-file&gt;**.TopicTerm.tab, **&lt;output-file&gt;**.top **&lt;number&gt;** Events.txt). **Note**: No extension must be given.
- **-N** : Factor that will be multiplied by PSIs to transform them into counts. Default = 1000.
- **-l** | **--four-lda-params**:  Latent Dirichlet Allocation parameters. **MIND THE ORDER**:
    - Number of topics/clusters (sklearn's *LatentDirichletAllocation* function)
    - Number of maximum iterations (sklearn's *LatentDirichletAllocation* function)
    - Number of iterations of feature selection (put to 0 if none wanted)
    - Number of least probable features (events) to check in each iteration of feature selection : this number of events with the worst probabilities will be be checked, if in common between clusters, they will be removed.
- **-t** | **--n-top-events**: Number of most probable events of every cluster to show in output file.
- **-p** | **--cutoff-prob**: Cumulated probability cutoff, so as to find the number of events needed to overcome this cutoff (total number of main events). Default = 0.8.
- **-r** | **--no-random-matrices**: Number of random matrices to create when evaluating clusters.
- **-y** | **--event-type**: Filter to keep only a single type of event. To choose among SE, A5, A3, MX, RI, AF, AL and all. Default = all.
- **-f** | **--total-filter**: Minimum mean for all samples of the total expression of the transcripts involved in the event. It will filter out the events that do not reach this total expression value for the transcripts defining the event (the denominator of the PSI calculation). Default = 1. /!\ Not the same as the one in psiPerEvent, which is a filter per sample.
- **-mi** | **--microexon-enrichment**: if option written in command, events will be filtered to keep only micro exons (<= 50 nt). Default = False.
- **-m** | **--mode**: verbose mode to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL.

**Feature selection** : for a given number of iterations, the last events in terms of probability for each cluster will be checked (number given by arguments). If some events are in common in all clusters, they will be removed. In figure 2, if we check the 3 least probable events for each cluster in the cluster-event distributions, event 4 is in common for all clusters and will be removed before recalculating the distributions.


### Input files

An ioe file and a "transcript expression file" are required as input. See SUPPA’s README for more information on the structure of these files (https://github.com/comprna/SUPPA#input-files-1, https://github.com/comprna/SUPPA#ioe).

### Output files

This tool outputs three files:
- The document-topic (sample-cluster) distribution given by the Latent Dirichlet Allocation model, with sample names as row names, no headers, tab-separated.
- The topic-term (cluster-event id) distribution given by the Latent Dirichlet Allocation model, with event IDs as headers, no row names, tab-separated.
- A file of the <n> most probable events of every cluster. Clusters preceded by a comment line `#Cluster_<cl-number>: <number-of-main-events>/<total-number-of-events> events needed to overcome a probability of <prob-cutoff>`. On each new line there is an event id separated by a tab from its probability in the cluster.

The name of the output is generated as follows:
- **&lt;output-file&gt;**.DocTopic.tab
- **&lt;output-file&gt;**.TopicTerm.tab
- **&lt;output-file&gt;**.top **&lt;n&gt;** Events.txt : **&lt;n&gt;** being the number of top events shown in output, given as an argument in the command.

#### An example of the first lines of a document-topic distribution output:
```
[Samples Cluster1    Cluster2    Cluster3]
GTEX-XXEK-1126-SM-4BRUX	  3.61907773e-07	  9.99989018e-01	  1.06201762e-05
GTEX-S32W-1926-SM-4AD63	  4.80354649e-07	  8.50639711e-01	  1.49359809e-01
GTEX-X4EP-1026-SM-4QAS5	  3.58454581e-07	  9.99999275e-01	  3.66085427e-07
GTEX-WQUQ-1926-SM-4OOSA	  6.03616895e-07	  9.99998640e-01	  7.56013399e-07
GTEX-14E1K-0326-SM-5S2PE	  4.06353591e-07	  9.99999179e-01	  4.14283757e-07
...
```
The line between brackets is only for clarity, it isn't present in the output.

#### A topic-term distribution output:
```
EventID1  EventID2  EventID3 EventID4   ...
<Event1-probability-in-cluster1>    <Event2-probability-in-cluster1>    <Event3-probability-in-cluster1>    <Event4-probability-in-cluster1>    ...
<Event1-probability-in-cluster2>    <Event2-probability-in-cluster2>    <Event3-probability-in-cluster2>    <Event4-probability-in-cluster2>    ...
<Event1-probability-in-cluster3>    <Event2-probability-in-cluster3>    <Event3-probability-in-cluster3>    <Event4-probability-in-cluster3>    ...
```
The event ids are the same as in SUPPA.

#### An example of top 4 events output for 3 clusters and a 0.8 cutoff to find the number of main events:
```
# Cluster_1: 4001/33359 events needed to overcome a probability of 0.8
ENSG00000197971.14;SE:chr18:76984894-76988495:76988527-76988877:-	0.0165711348953
ENSG00000131095.11;SE:chr17:44913430-44913728:44913823-44914028:-	0.00737414082667
ENSG00000131095.11;SE:chr17:44913823-44914028:44914088-44915026:-	0.00734398832255
ENSG00000197971.14;SE:chr18:76980471-76984775:76984894-76988877:-	0.00650691042543
# Cluster_2: 3787/33359 events needed to overcome a probability of 0.8
ENSG00000168484.12;SE:chr8:22161870-22162574:22162732-22163080:+	0.0141654220066
ENSG00000133112.16;SE:chr13:45338776-45339497:45339602-45339994:-	0.010038983434
ENSG00000229117.8;SE:chr12:56116799-56117189:56117211-56117482:+	0.00804762123155
ENSG00000147403.16;SE:chrX:154399396-154399487:154399594-154399803:+	0.00773487101274
# Cluster_3:  1858/33359 events needed to overcome a probability of 0.8
ENSG00000163631.16;SE:chr4:73416353-73417531:73417669-73418088:+	0.0670261311022
ENSG00000163631.16;SE:chr4:73419639-73420254:73420321-73421092:+	0.0638619242146
ENSG00000158874.11;SE:chr1:161223050-161223350:161223425-161223595:-	0.0271666946849
ENSG00000158874.11;SE:chr1:161222522-161222918:161223050-161223350:-	0.0269758417982
```

## Sample-cluster visualisation - Plot_DocumentTopic.R
----------------------------

Developed with R

### Modules used
- CountClust : Copyright (c) Dey K, Hsiao J and Stephens M (2016)
- optparse
- RColorBrewer


### Description 
This script makes the structure histogram plot of the document-topic distribution (sample-cluster distribution), with samples grouped by labels.

### Command
```
Rscript Plot_DocumentTopic.R -i <document-topic-file> -o <output-plot-name> -f <sample-label-file>
```

### Command options
- **-i** | **--input**: Document-topic distribution file, with row names, without headers.
- **-o** | **--output**: Name of the output plot. Default = plot.pdf
- **-f** | **--file**: Tab-separated sample-label file with headers, containing all samples (documents) from the document-topic distribution file.

### Input files
- Document-topic distribution file : tab-separated distribution for each document (sample), with row names being document names. No headers. - *output of LatentEvents_LDA.py*
- Sample-label file : Tab-separated file with headers. On each line we have a sample and its label. At least all samples from the document-topic distribution must be present.

#### Example of sample-label file
```
Samples Tissues
GTEX-XXEK-1126-SM-4BRUX Liver
GTEX-S32W-1926-SM-4AD63	Liver
GTEX-X4EP-1026-SM-4QAS5	Liver
GTEX-WQUQ-1926-SM-4OOSA	Liver
GTEX-14E1K-0326-SM-5S2PE	Liver
GTEX-WFON-1726-SM-4LVMQ	Liver
GTEX-13OW6-2626-SM-5IFF2	Liver
GTEX-O5YT-0826-SM-3TW8N	Liver
GTEX-XBEC-1526-SM-4AT68	Liver
GTEX-13QJC-0726-SM-5RQJK	Liver
GTEX-11TT1-1726-SM-5EQLJ	Liver
GTEX-11UD2-1626-SM-5EQM3	Liver
GTEX-XXEK-0626-SM-4BRWE	Lung
GTEX-111VG-0726-SM-5GIDC	Lung
GTEX-N7MT-0126-SM-26GMB	Lung
GTEX-SN8G-0926-SM-4DM5I	Lung
GTEX-WWYW-0926-SM-3NB2Z	Lung
GTEX-13O3Q-0526-SM-5KM18	Lung
...
```

### Output file
A structure histogram plot by the name **&lt;output-plot&gt;**, grouping all samples by labels. Samples in each batch will have been sorted by the proportional membership of the most representative cluster in that batch.

![StructurePlot.jpg](https://github.com/comprna/LatentEvents/blob/master/Figures/SampleClusterDistrib_plot.jpg)
**Fig.3** - Structure plot of GTEX data. 503 samples : 110 from liver, 288 from lung, 105 from brain. *LatentEvents_LDA.py* run for 150 max iterations, 3 clusters, all events, expression filter of 1, no feature selection.

## Cluster-event comparison - TopicTerm_Comparison.R
----------------------------

Developed with R

### Modules used
- optparse
- ggplot2
- reshape2

### Description
This script plots the specificity of a cluster against the others for events in decreasing probability for said cluster. The specificity is the log2 probability ratio of a cluster against another one. Each cluster will have an output plot with a curve for each of the other clusters.
### Command line
```
Rscript TopicTerm_Comparison.R -i <topic-term-file> -o <output-name> -t <x-limit>
```
**Command example**
```
Rscript --vanilla TopicTerm_Comparison.R -i Example.TopicTerm.tab -o SpecificityProbability_plot -t 1000
```
### Command options  
- **-i** | **--input**: Topic-term distribution file
- **-o** | **--output**: Common name of the output plots, that will differ with the extension. Default = plot -> outputs = plot\_cl **&lt;number&gt;** .pdf. **Note** : No extension must be given.
- **-t** | **--topEvents**: Number of most probable events to plot (limit of the horizontal line of the plot). Default = all events.

### Input file
- Topic-term distribution file : tab-separated distribution for each topic (cluster) of the terms (splicing events), with headers being event ids. No row names.  - *output of LatentEvents_LDA.py*


### Output files
Outputs by the name **&lt;output-name&gt;**\_cl **&lt;cluster-number&gt;** .pdf. The cluster number will be automatically generated, only the output name without extension can be given. Each file will have a specificity-probability plot. 

The specificity is the log2 ratio of the probability of the cluster against the other clusters (there will be a curve for each of the other clusters), for the events in decreasing order of probability of the cluster.
There will therefore be as many outputs as clusters (topics).
```
Example of output plot
```
