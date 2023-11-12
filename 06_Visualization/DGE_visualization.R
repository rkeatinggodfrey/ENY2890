## Project: Fall 2023 CURE Course 
## Authors: R. Keating Godfrey
## Last updated: 23-11-11

install.packages("dplyr","ggplot2","ggrepel")
library(dplyr)
library(ggplot2)
library(ggrepel)

## Pages used to create this script

## Volcano plots:
## https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-counts-to-viz-in-r/tutorial.html
## https://stackoverflow.com/questions/15624656/label-points-in-geom-point
## https://ggrepel.slowkow.com/articles/examples.html
## http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines


## Use threshold of >1 Log2FC and padj <0.05

##----- Female geni vs female genitalia -----##
## subset results table to only fields of interest
viz.results <- results.annotated[,c(1,2,3,7,9,10,14,18,27)]
## create a dataframe limited to ony significant targets
sig.results <-subset(viz.results,viz.results$padj <0.05)

write.csv(sig.results, "Sig_Results_Geni_Log2FC.csv")

## Volcano Plot of all targets 
ggplot(data = sig.results, aes(x = log2FoldChange, y = -log10(padj), 
                                label = geneID))+
  geom_point(colour="purple")+
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
  theme_classic()


## Volcano Plot of sig targets with protein family labels
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
