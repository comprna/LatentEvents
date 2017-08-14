# -*- coding: utf-8 -*-
"""
@author: Adela Weber
@email: adela[at]nicoweb.com
"""


import sys
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
from lib.tools import *
from sklearn.decomposition import LatentDirichletAllocation #Scikit-learn version 0.18.2
import numpy as np


description = \
    "Description:\n\n" + \
    "This tool reads an ioe file and a transcript expression file and does\n" + \
    "a Latent Dirichlet Allocation on a given event type or all event types for different samples.\n"
    
parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=True)

parser.add_argument("-i", "--ioe-file", help="Input file with the event-transcripts equivalence (.ioe format).", required=True)
parser.add_argument("-e", "--expression-file", required=True,
                    help="Input transcript expression file.")
parser.add_argument("-o", "--output-file", required=True,
                    help="Common output name for LDA matrices and best events per cluster files.")
parser.add_argument("-N", default = 1000, type = int,
                    help='Number to multiply PSIs (for the transformation into counts).\n[Default = 1000]')
parser.add_argument("-l", "--four-lda-params", nargs='+', required = True, type=int, 
                    help="Latent Dirichlet Allocation parameters, MIND THE ORDER:\n\t- Number of topics\n\t- Number of maximum iterations\n\t- Number of iterations of feature selection\n\t- Number of least probable features to check in each iteration\nFor no feature selection, put no. of iterations or/and features to check to 0")
parser.add_argument("-t", "--n-top-events", required = True, type = int,
                    help="Number of most probable events of every cluster to show in output file")
parser.add_argument("-p", "--cutoff-prob", default = 0.8, type = float,
                    help='Cumulated probability cutoff (so as to find the number of events needed to overcome it - total number of main events).\n[Default = 0.8]')
parser.add_argument("-r", "--no-random-matrices", required = True, type = int,
                    help="Number of random matrices to create when evaluating clusters")
parser.add_argument("-y", "--event-type", default = "all", type = str, choices = ["SE", "A5", "A3", "MX", "RI", "AF", "AL", "all"],
                    help='Filter to keep only a single type of event. To choose among SE, A5, A3, MX, RI, AF, AL, all.\n[Default = all]')
parser.add_argument("-f", "--total-filter", type=float, default = 1,
                    help="Minimum mean for all samples of the total expression of the transcripts involved in the event. [Default = 1]")
parser.add_argument("-mi", "--microexon-enrichment", action='store_true',
                    help='Filter to keep microexons.\n[Default = False]')
parser.add_argument("-m", "--mode", default="INFO",
                    help="to choose from DEBUG, INFO, WARNING, ERROR and CRITICAL")




def Distributions_writer (Matrix, OutputName, MatrixType, headers) :
    """
    Creates a tabulated matrix file (.tab) with row names for document-topic matrix
    but column names for topic-term matrix.
    
    Matrix - LDA output
    OutputName - name of the output file
    MatrixType - either "document-topic" or "topic-term", defines distribution matrix content
    headers - list of labels : First line of output file, corresponds to column names
        for topic-term matrix but to row names for document-topic matrix
    """
    if MatrixType == "document-topic": #No headers but row names
        OutputName += ".DocTopic.tab"
        try :
            f = open(OutputName, "w")
            logger.info("Writing document-topic matrix in %s" %OutputName)
            for samp_idx in range(len(Matrix)):
                f.write(headers[samp_idx] + "\t" + np.array2string(Matrix[samp_idx], separator="\t").replace("\n", "\t")[1:-1] + "\n")
            f.close()
            logger.info("%s closed." %OutputName)
        except BaseException:
             logger.error("Error when writing the output : %s" % sys.exc_info()[1].args[0])
             sys.exit(1)
    elif MatrixType == "topic-term":
        OutputName += ".TopicTerm.tab"
        logger.info("Writing topic-term matrix in %s" %OutputName)
        try :
            np.savetxt(OutputName, Matrix, delimiter = "\t", header = ("\t").join(headers), comments = "")
        except BaseException:
             logger.error("Error when writing the output : %s" % sys.exc_info()[1].args[0])
             sys.exit(1)


def PSI_to_counts(psiMatrix, N):
    """
    Transforms psi matrix to 'counts' matrix to run LDA : PSIs multiplied by N and transformed into integers.
    """
    return (psiMatrix * N).astype(int)



def number_main_events_cluster(Vector_TopicTerm, TopicTerm_decr_idx, prob):
    """
    Gives the number of events for the given cluster needed to overcome prob
    cumulative event probability
    
    Vector_TopicTerm - topic term distribution
    TopicTerm_decr_idx - original order indexes of a topic-term distribution
        sorted in decreasing order for a given cluster
    prob - probability cutoff (between 0 and 1)
    """
    
    Nb_mainEvents = 0
    Cum_prob = 0 #Cumulative event probability
    while Cum_prob < prob:
        Cum_prob += Vector_TopicTerm[TopicTerm_decr_idx[Nb_mainEvents]]
        Nb_mainEvents += 1
    return Nb_mainEvents


def n_top_events_writer (TopicTerm, eventsIDs, n, OutputName, prob):
    """
    Gives the n most probable events of every cluster in an output file, with a "#" followed by cluster name followed by the
    number of main events and then in new lines the top n events ids (in decreasing order of relevance) with their probability (separated by a tab)
    
    TopicTerm - Topic-term distribution matrix
    eventsIDs - Names of events
    n - Number of most representative events you want in output
    OutputName - output name of the file containing these results
    prob - probability cutoff for number_main_events_cluster function
    """
    try:
        OutputName += ".top" + str(n) + "Events.txt"
        f = open(OutputName, "w")
        logger.info("Writing best events for each cluster in %s." %OutputName)
        for clust in range(len(TopicTerm)):
            TpcTrm_idx_decr = (TopicTerm[clust]).argsort()[::-1] #Original order indexes of a topic-term distribution sorted in decreasing order for a cluster
            main_events = number_main_events_cluster(TopicTerm[clust], TpcTrm_idx_decr, prob) #Number of relevant events
            Events_max_idx = TpcTrm_idx_decr[:n]
            f.write("# Cluster_" + str(clust + 1) + ": " + str(main_events) + "/" + str(len(TopicTerm[clust])) + " events needed to overcome a probability of " + str(prob) + "\n")
            for idx in Events_max_idx:
                f.write(eventsIDs[idx] + "\t" + str(TopicTerm[clust][idx]) + "\n")
        f.close()
        logger.info("%s closed." %OutputName)
    except BaseException:
        logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
        sys.exit(1)


def LDA (CountsMatrix, K, no_iter) :
    """
    Fits the Latent Dirichlet Allocation model and returns the model, the topic-term distribution and the document-topic distribution
    
    CountsMatrix - Counts matrix (document x terms matrix)
    K - number of topics (clusters) (sklearn's LDA function parameter)
    no_iter - maximum number of iterations (sklearn's LDA function parameter)
    """
    lda = LatentDirichletAllocation(n_topics = K, learning_method='online', max_iter = no_iter).fit(CountsMatrix) #Model fitting
    
    #components_ : topic-term distribution matrix (no. topics x no. events)
    TopicTerm = lda.components_ /lda.components_.sum(axis = 1)[:, np.newaxis] #Normalization to have distributions (sum = 1)
       
    #transform method : document-topic distribution matrix (no. samples x no. topics), normalized lines for scikit-learn versions >= 0.18
    DocTopic = lda.transform(CountsMatrix)
        
    if DocTopic.sum(axis = 1).sum() != len(DocTopic) :
        logger.error("Unnormalized lines in LDA transform method. Check your scikit-learn version (0.18 needed).")
        sys.exit(1)
    return lda, TopicTerm, DocTopic


def EventSelection(Countmatrix, eventIDs, TopicTerm, check_last) :
    """
    Takes the check_last least probable events from each cluster and checks which events are in common. 
    These events will be eliminated from the count matrix and event_id list.
    
    Countmatrix - numpy array. Counts matrix (Document x terms matrix)
    eventIDs - list of the event IDs, in the same order as the matrix (terms = events)
    TopicTerm - Topic term distribution matrix
    check_last - number of least probable events we want to check (> 0)
    """
    Break = False #If no events in common, do not continue iterating (outside function)
    Common_events_idx = (TopicTerm[0]).argsort()[::-1][-check_last:] #check_last last original order indexes of a topic-term distribution sorted in decreasing order for a cluster (first cluster here)
    for clust in range(1, len(TopicTerm)):
        TpcTrm_idx_decr = (TopicTerm[clust]).argsort()[::-1][-check_last:] #check_last last original order indexes of a topic-term distribution sorted in decreasing order for a cluster
        Common_events_idx = np.intersect1d(Common_events_idx, TpcTrm_idx_decr, assume_unique = True) #returns common indexes in original order of topic-term distribution
    if len(Common_events_idx) > 0:
        logger.info("%d common least probable events" %len(Common_events_idx))
        Countmatrix = np.delete(Countmatrix, Common_events_idx, axis = 1)
        events_removed = 0
        for ev in Common_events_idx:
            eventIDs.remove(eventIDs[ev-events_removed])
            events_removed += 1
    else:
        Break = True
    return Break, eventIDs, Countmatrix


"""
#TEST EventSelection
Count = np.array([[ 11., 2., 4,   3.,  1., 0],
                  [ 4.,  2., 7,   5.,  1., 0],
                  [ 7.,  2., 11., 5.,  1., 0],
                  [ 11., 2., 6,   3.,  1., 0],
                  [ 8.,  2., 7,   5.,  1., 5],
                  [ 9.,  2., 11., 12., 1., 7]])

eventsIDS = ['A', 'B', 'C', 'D', 'E', 'F']

TopicTerm = Count[:]
TopicTerm = np.transpose(np.transpose(TopicTerm)/np.sum(TopicTerm, axis = 1))
EventSelection(Count, eventsIDS, TopicTerm, 3)
"""


def RandomizeMatrix(matrix):
    """
    Gives a copy of matrix with shuffled values for each line
    
    matrix - numpy array. Counts matrix (Document x terms matrix)
    """
    RandomMatrix = np.ones(np.shape(matrix))
    for row in range(len(RandomMatrix)):
        RandomMatrix[row] = np.random.permutation(matrix[row]) #Shuffles among values of line
    return RandomMatrix


def RandScore(CountsMatrix, K, no_iter):
    """
    Calculates score for observed data with LDA model fitted to randomized matrix
    
    CountsMatrix - numpy array. Counts matrix (Document x terms matrix) for our real data
    K - number of clusters
    """
    #Randomize CountsMatrix
    RandMatrix = RandomizeMatrix(CountsMatrix)
    #LDA for randomized matrix
    lda_rand = LatentDirichletAllocation(n_topics = K, learning_method='online', max_iter = no_iter).fit(RandMatrix) #Model fitting 
    return lda_rand.score(CountsMatrix)


def EvaluateModel(CountsMatrix, K, lda_obs, no_randMatrices, no_iter):
    """
    Calculates score of observed data matrix with models from randomized matrices against
    score of data matrix with model of matrix, returns the scores and a p-value.
    
    CountsMatrix - numpy array. Counts matrix (Document x terms matrix) for our real data
    K - number of clusters 
    lda_obs - Latent Dirichlet Allocation model fitted to our observed matrix
    no_randMatrices - number of random matrices we want to create from our observed matrix
    no_iter - maximum number of iterations for LDA (LDA function argument)
    """
    Scores_rand = np.zeros((no_randMatrices))
    logger.info("Randomizing matrices and fitting Latent Dirichlet Allocation models...")
    for mat in range(no_randMatrices):
        Scores_rand[mat] = RandScore(CountsMatrix, K, no_iter)
    logger.info("Finished with randomizing matrices.")
    Score_obs = lda_obs.score(CountsMatrix)
    pval = np.count_nonzero(Score_obs < Scores_rand)/float(no_randMatrices)
    return pval, Scores_rand, Score_obs



def main():
    
    args = parser.parse_args()
    
    #Parsing arguments
    mode = "logging." + args.mode

    #Setting logging preferences
    logger = logging.getLogger(__name__)
    logger.setLevel(eval(mode))

    #Setting the level of the loggers in lib
    setToolsLoggerLevel(mode)
    
    #Parameters
    if len(args.four_lda_params) == 4:
        no_topics, no_iter, no_feature_sel, check_last = args.four_lda_params
    else:
        logger.error("Wrong number of LDA parameters")
        sys.exit(1)
    
    if check_last == 0: #No feature selection wanted - no. of feature selection iterations set to 0
        no_feature_sel = 0
    
    output_file = args.output_file
    expression_file = args.expression_file
    event_type = args.event_type
    
    np.seterr(divide='raise') #Raise error in zero numpy divisions

    expression_dictionary = {}  # to store transcript expression [transcript_id] = {[colId] = expression}
    col_ids = []  # to store all the column_id for the expression fields (samples)

    try:

        #BUFFERING EXPRESSION FILE (TPMs)
        factory = FactoryReader()
        r = factory.getReader("expression") #Expression file (TPM) type
        r.openFile(expression_file)
        logger.info("Buffering transcript expression levels.")
        line = r.readLine()
        
        try:
            while True:
                try:
                    arguments = nextel(line)
                    col_ids = arguments["colIds"] #Name of the samples (columns of expression file)
                    expression_dictionary[arguments["transcript_id"]] = {}
                    # Assign for each colId the expression [colId] = expression
                    for key, value in arguments["expression"].items():
                        expression_dictionary[arguments["transcript_id"]][key] = float(value)
                        
                except ValueError:
                    logger.error("%s expression is not a float. Skipping..." % arguments["transcript_id"])
                    continue
                
        except StopIteration:
            if not expression_dictionary:
                logger.error("No expression values have been buffered.")
                sys.exit(1)    

                
            #BUFFERING IOE - to filter events (only above a certain expr.), to build the LDA "count" (psi) matrix
            factory = FactoryReader()
            r = factory.getReader("ioe")
            r.openFile(args.ioe_file)
            logger.info("Filtering events from ioe file : %s events above minimal expression kept." %event_type) 
            line = r.readLine()
            
            events_ids = [] #list of events IDs
            CountsMatrix = [] #Initialisation of count matrix
            CountSkip = 0 #Skipped events counter
            CountFilter = 0 #Unsufficiently expressed events counter
            
            Microexon = args.microexon_enrichment
            if Microexon:
                logger.info("Keeping only events regarding micro exons...")
            try:
                while True: #for each event
                    arguments = nextel(line)
                    
                    if event_type != "all" and ";" + event_type + ":" not in arguments["event_id"]: #Only for given event type
                        continue

                    Total_transcripts_event = np.zeros((len(col_ids), 1)) # Vector of cumulative expession of all transcripts with current event for each sample (column)
                    CountPerEvent = np.zeros((len(col_ids),1)) #Vector of cumulative sum of alternative transcript TPMs of current event per sample

                    if arguments["event_id"] in events_ids :
                        logger.error("Duplicated event %s. Skipping line..." %
                                     arguments["event_id"])
                        continue
                    
                    #FILTERING EVENTS TO ENRICH MATRIX WITH MICROEXON EVENTS
                    if Microexon :
                        #We enrich our events with micro exons - filter
                        Exon_lims = str.split(str.split(arguments["event_id"], "-")[1], ":") #List of (str) exon start and ending positions
                        if abs(int(Exon_lims[1]) - int(Exon_lims[0])) > 50: #skip exons superior to 50 nt
                            continue
                        


                    skip = False  # to avoid checking total_iso if event == "NA"
                    
                    for tr in arguments["alt_iso"].rstrip("\n").split(","): #Transcripts containing the exon of the event
                        try:
                            for x in range(len(col_ids)):
                                CountPerEvent[x,] += expression_dictionary[tr][col_ids[x]] #cumulative expression per sample
                        except KeyError:
                            logger.error(
                                "Transcript %s not found in the \"expression file\" ." %
                                tr)
                            logger.error("Event %s skipped in LDA matrix." %
                                         arguments["event_id"])
                            skip = True
                            break
                        
                    if not skip:
                        for tr in arguments["total_iso"].rstrip("\n").split(","):
                                try:
                                    for x in range(len(col_ids)) :
                                        Total_transcripts_event[x,] += \
                                            expression_dictionary[tr][col_ids[x]]
    
                                except KeyError:
                                    logger.error(
                                        "transcript %s not found in the \"expression file\"." %
                                        tr)
                                    logger.error("Event %s skipped in LDA matrix." %
                                                 arguments["event_id"])
                                    skip = True
                                    break
                                
                        if np.mean(Total_transcripts_event) >= args.total_filter and not skip: #mean global expression per event > filter_cutoff. 
                        #/!\ Not the same filter as in SUPPA's psiCalculator
                            #If it passes the filter and we are not skipping the event
                            #Calculating PSI for current event
                            try:
                                CountsMatrix.append(np.divide(CountPerEvent, Total_transcripts_event)) #adding PSI vector to matrix
                                events_ids.append(arguments["event_id"])
                            except FloatingPointError:
                                logger.error("Zero division")

                        elif skip:
                            CountSkip += 1
                        else:
                            CountFilter += 1
                    
                    else:
                        CountSkip += 1

            except StopIteration:
                logger.info("%d events kept." % len(events_ids))
                logger.info("%d events insufficiently expressed." % CountFilter)
                logger.info("%d events skipped." % CountSkip)
                
                #BUILDING "COUNTS" MATRIX FOR LDA FITTING (keeping integers)
                CountsMatrix = np.hstack(CountsMatrix) #Transforming into numpy array
                N = args.N
                CountsMatrix = PSI_to_counts(CountsMatrix, N)
                    
                #LATENT DIRICHLET ALLOCATION
                #Model fitting
                logger.info("Fitting Latent Dirichlet Allocation model...")
                lda, TopicTerm, DocTopic = LDA (CountsMatrix, no_topics, no_iter)
                
                #Feature selection
                for featureselection in range(no_feature_sel): #Filter least probable events common to all clusters
                    logger.info("Iteration %s : Eliminating least probable events common to all clusters" %(str(featureselection + 1)+"/"+str(no_feature_sel)))
                    Break, events_ids, CountsMatrix = EventSelection(CountsMatrix, events_ids, TopicTerm, check_last) #Filter
                    if Break == True: #If CountsMatrix hasn't changed
                        logger.error("No common events in the least probable events of all clusters for iteration %d" %no_feature_sel)
                        break
                    logger.info("Fitting Latent Dirichlet Allocation model...")
                    lda, TopicTerm, DocTopic = LDA(CountsMatrix, no_topics, no_iter) #Update with new CountsMatrix
                
                logger.info("%d final events." % len(events_ids))
                logger.debug("Shape of doc-topic matrix %s" %str(np.shape(DocTopic)))
                logger.debug("Shape of topic-term matrix %s" %str(np.shape(TopicTerm)))
                
                #WRITING OUTPUT
                Distributions_writer (DocTopic, output_file, "document-topic", col_ids) #Document-topic output
                Distributions_writer (TopicTerm, output_file, "topic-term", events_ids) #Topic-term output


                #n MOST PROBABLE EVENTS PER CLUSTER AND NUMBER OF EVENTS NEEDED TO OVERCOME prob% OF CUMULATED PROBABILITY
                n = args.n_top_events
                prob = args.cutoff_prob #cumulated probability cutoff we want to overcome
                n_top_events_writer (TopicTerm, events_ids, n, output_file, prob)
                
                #CLUSTER EVALUATION
                no_rand_matrix = args.no_random_matrices
                logger.info("Evaluating model...")
                pval, Scores_rand, Score_obs = EvaluateModel(CountsMatrix, no_topics, lda, no_rand_matrix, no_iter)
                logger.debug("Scores of randomized matrices models :\n%s" %np.array2string(Scores_rand))
                logger.debug("Scores of counts matrix model : %f" %Score_obs)
                logger.info("P-value : %f" %pval)
           
                
                
    except BaseException:
        logger.error("Unknown error: %s" % sys.exc_info()[1].args[0])
        sys.exit(1)
    logger.info("Done")

if __name__ == '__main__':
    main()