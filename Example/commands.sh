# ensembl_hg38.SEevents.ioe is and ioe file enriched with SE events and some A3
# Tissues15.tpm gives the TPMs of 15 GTEX lung samples, 15 GTEX brain samples and 15 GTEX liver samples
# NamesTissues15.txt a sample-tissue file


python LatentEvents_LDA.py -i ensembl_hg38.SEevents.ioe -e Tissues15.tpm -o 3Tissues_3cl_Example -t 5000 -r 50 -l 3 150 0 0 -y SE

Rscript --vanilla Analysis/Plot_DocumentTopic.R -i 3Tissues_3cl_Example.DocTopic.tab -o StructurePlot_Example.pdf -f NamesTissues15.txt

Rscript --vanilla Analysis/TopicTerm_Comparison.R -i 3Tissues_3cl_Example.TopicTerm.tab -o SpecificityProbabilityPlot_Example -t 1000
