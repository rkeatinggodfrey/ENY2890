## Project: Fall 2023 CURE Course 
## Authors: R. Keating Godfrey
## Last updated: 23-11-12

install.packages("dplyr","ggplot2","ggrepel")
library(dplyr)
library(ggplot2)
library(ggrepel)

## if it's not already in your Global Environment, import the
## annotated data frame
results.annotated <- read.csv("Genitalia_Males_v_Females_Results_table_02.csv", 
                              header=T)

## create a dataframe limited to significant targets
results.sig <-subset(results.annotated,results.annotated$padj <0.05)
write.csv(results.sig, "Sig_Results_Geni_Log2FC.csv")


################ Volcano plots ################

## Pages used to create this script
## https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-counts-to-viz-in-r/tutorial.html
## https://stackoverflow.com/questions/15624656/label-points-in-geom-point
## https://ggrepel.slowkow.com/articles/examples.html
## http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines


## Use threshold of >1 Log2FC and padj <0.05

##----- Female vs. Male genitalia -----##
## subset results table to only fields of interest
viz.results <- results.sig[,c(1,2,3,7,9,10,14,15,18,27)]

## Volcano Plot of all targets 
ggplot(data = viz.results, aes(x = log2FoldChange, y = -log10(padj), 
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
  ggtitle("Male vs Female Genitalia")+
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




################ GO Enrichment ################

## Resources:
## https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
## https://archetypalecology.wordpress.com/2021/01/27/how-to-perform-kegg-and-go-enrichment-analysis-of-non-model-species-using-r/
## https://rstats101.com/separate-a-collapsed-column-into-multiple-rows/
## https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/


require(DOSE)
require (clusterProfiler)


## As a first example, let's just look at enrichment terms across all differentially
## expressed genes (this will include male and female, up or down)

## get a list of differential expressed genes
## from the significant results table
deg.genes <-results.sig$geneID

## Now we need to get GO terms associated with those genes
## some genes have mulltiple GOs, so let's makes those show up 
## as separate rows

## define an object called "go" that's really just the results.annoated table
go <- results.annotated
## remove empty entries
go <- subset(kegg,kegg$GOs !="-")
## remove the "GO:" from in front of every GO term and replace with "GO"
go$GOs <- gsub("GO:","GO",as.character(go$GOs))
## now if a gene has multiple GO terms associated with it
## put each on a new row with that gene ID
go <- go %>% separate_rows(GOs)

## create a data frame of GO to gene id to use in enrichment analysis
go2gene <-go[, c("GOs","geneID")]

## run enrichment analysis (in the clusterProfiler package)
enriched <- enricher(deg.genes, TERM2GENE=go2gene)
head(summary(enriched))
barplot(enriched)



## As a second example, let's just look at FEMALE enrichment terms 

## Because we compared Males to Females, genes with
## NEGATIVE Log2FC values were detected at higher expression
## levels in females
## So you can get a list of these by subsetting the significant results to
## those with Log2FC < 0
l2fc.f <-subset(results.sig,results.sig$log2FoldChange < 0)
deg.genes.f <-l2fc.f$geneID

## The GO terms associated with those genes are the same as above,
## so this part of the script to get a go2gene data frame stays the same
## you don't have to run this again if you already did. 

## some genes have mulltiple GOs, so let's makes those show up 
## as separate rows
## define an object called "go" that's really just the results.annoated table
go <- results.annotated
## remove empty entries
go <- subset(kegg,kegg$GOs !="-")
## remove the "GO:" from in front of every GO term and replace with "GO"
go$GOs <- gsub("GO:","GO",as.character(go$GOs))
## now if a gene has multiple GO terms associated with it
## put each on a new row with that gene ID
go <- go %>% separate_rows(GOs)

## create a data frame of GO to gene id to use in enrichment analysis
go2gene <-go[, c("GOs","geneID")]

## FEMALE enrichment analysis (in the clusterProfiler package)
enriched.f <- enricher(deg.genes.f, TERM2GENE=go2gene)
head(summary(enriched.f))
barplot(enriched.f)
