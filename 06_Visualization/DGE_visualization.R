## Project: Fall 2023 CURE Course 
## Authors: R. Keating Godfrey
## Last updated: 23-11-12

install.packages("dplyr","ggplot2","ggrepel")
library(dplyr)
library(ggplot2)
library(ggrepel)


## create a dataframe limited to ony significant targets
sig.results <-subset(results.annotated,results.annotated$padj <0.05)
write.csv(sig.results, "Sig_Results_Geni_Log2FC.csv")


################ Volcano plots ################

## Pages used to create this script
## https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-counts-to-viz-in-r/tutorial.html
## https://stackoverflow.com/questions/15624656/label-points-in-geom-point
## https://ggrepel.slowkow.com/articles/examples.html
## http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines


## Use threshold of >1 Log2FC and padj <0.05

##----- Female vs. Male genitalia -----##
## subset results table to only fields of interest
viz.results <- sig.results[,c(1,2,3,7,9,10,14,15,18,27)]

## Volcano Plot of all targets 
ggplot(data = sig.results, aes(x = log2FoldChange, y = -log10(padj), 
                                label = geneID))+
  geom_point(colour="purple")+
  xlim(-35,35)+ # this sets the scale for the x-axis
  theme_classic()

## Volcano Plot of all targets with significance designated
ggplot(data = sig.results, aes(x = log2FoldChange, y = -log10(padj), 
                           label = geneID))+
  geom_point(color=dplyr::case_when(sig.results$padj > .0001 ~ "grey",
                                    sig.results$log2FoldChange > 3 ~ "blue",
                                    sig.results$log2FoldChange < -3 ~ "blue"),
             size = 2, alpha = 0.65)+
  geom_vline(xintercept=3, colour="blue", linetype="dashed")+
  geom_vline(xintercept=-3, colour="blue", linetype="dashed")+
  geom_hline(yintercept=4, colour="blue", linetype="dashed")+
  ggtitle("Female vs. Male Genitalia")+
  xlim(-35,35)+
  theme_classic()


## Volcano Plot of sig targets with protein family labels
## here the "label =" parameter is set to PFAMs so that protein
## family names will be associated with data points
## You could use "geneID" or "Preferred_name" for this
ggplot(data = sig.results, aes(x = log2FoldChange, y = -log10(padj), 
                                 label = PFAMs))+
  geom_point(color=dplyr::case_when(sig.results$padj > .0001 ~ "grey",
                                    sig.results$log2FoldChange > 3 ~ "blue",
                                    sig.results$log2FoldChange < -3 ~ "blue"),
             size = 2, alpha = 0.65)+ 
  geom_text_repel(data = subset(sig.results, log2FoldChange > 3
                                & sig.results$padj < .0001),
                  size          = 2.5,
                  box.padding   = 0.5,
                  point.padding = 0.05,
                  force         = 20,
                  segment.size  = 0.1,
                  segment.color = "grey",
                  max.overlaps  = 20)+
  geom_text_repel(data = subset(sig.results, log2FoldChange < -3
                                & sig.results$padj < .0001),
                  size          = 2.5,
                  box.padding   = 0.5,
                  point.padding = 0.05,
                  force         = 20,
                  segment.size  = 0.1,
                  segment.color = "grey",
                  max.overlaps  = 20)+
  ggtitle("Female vs Male Genitalia")+
  theme_classic()

################ End Volcano plots ################




################ Box plots of counts data ###########



################ Heatmaps of counts data ###########




################ KEGG and GO Enrichment ################




##
## https://archetypalecology.wordpress.com/2021/01/27/how-to-perform-kegg-and-go-enrichment-analysis-of-non-model-species-using-r/
## https://rstats101.com/separate-a-collapsed-column-into-multiple-rows/
## https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
