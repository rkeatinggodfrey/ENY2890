## Project: Fall 2023 CURE Course 
## Authors: R. Keating Godfrey
## Last updated: 23-11-12




################ Volcano plots ################

install.packages("dplyr","ggplot2","ggrepel")
library(dplyr)
library(ggplot2)
library(ggrepel)

## Pages used to create this script
## https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-counts-to-viz-in-r/tutorial.html
## https://stackoverflow.com/questions/15624656/label-points-in-geom-point
## https://ggrepel.slowkow.com/articles/examples.html
## http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines

## if it's not already in your Global Environment, import the
## annotated data frame
results.annotated <- read.csv("Genitalia_Males_v_Females_Results_table_02.csv", 
                              header=T)

## create a dataframe limited to significant targets
results.sig <-subset(results.annotated,results.annotated$padj <0.05)
write.csv(results.sig, "Sig_Results_Geni_Log2FC.csv")

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




################ Heatmaps of DESeq2 data ###########

## here are the packages you need (you don't need to reload if loaded)
library(DESeq2)
library(ggplot2)

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

## Resources
## https://youtu.be/S2_FTg9kaZU?feature=shared 


#################### Example heat map of genes with ########################
############################ L2FC > 2 or < -2  #############################

## turn the results from DESeq2 into a dataframe
## here the row names should be the gene IDs just like
## in the annotation file
sig <- as.data.frame(results.gen) # data frame of all results from DESeq2
sig <- subset(sig,sig$padj < 0.05) # subset to significant results

## Normalized counts from DESeq function (dds.gen.de in our example)
matrix.sig <-counts(dds.gen.de, normalized = T)[rownames(sig),]

## get zscore for each row
matrix.z <-t(apply(matrix.sig,1,scale))

## now take sample IDs from your metadata and make
## them the column names
colnames(matrix.z) <- meta$sample

## Ok let's make a heatmap!!
Heatmap(matrix.z, cluster_rows=T, cluster_columns=T, column_labels=colnames(matrix.z),
        name="z score") 

## That is terrifying ^^ ##


#################### Example heat map of genes with ########################
############################ L2FC > 2 or < -2  #############################

## So what should we do?
## How can we reduce the number of genes on our list?

## (1) Use a threshold for log2fc
sig.p2 <-subset(sig,sig$log2FoldChange >2)
sig.n2 <-subset(sig,sig$log2FoldChange < -2)

sig.2 <- rbind(sig.p2,sig.n2)

## Normalized counts from 
matrix.sig.2 <-counts(dds.gen.de, normalized = T)[rownames(sig.2),]

## get zscore for each row
matrix.z <-t(apply(matrix.sig.2,1,scale))

## now take sample IDs from your metadata and make
## them the column names
colnames(matrix.z) <- meta$sample

## Ok let's make a heatmap!!
Heatmap(matrix.z, cluster_rows=T, cluster_columns=T, column_labels=colnames(matrix.z),
        name="z score") 

#################### END example heat map of genes with ########################
############################ L2FC > 2 or < -2  #############################


#################### Example heat map of genes with ########################
###################### "odor" in their description #########################

## (2) Use GOs or other identifier from annotation file to subset results
## read in the annotation file made from coding sequences
## but in this case make the first column (gene IDs) the row names
annot <- read.csv("rna.annot.cds.emapper.annotations.geneID.csv", header=T,
                  row.names = 1)

sig.annot <- transform(merge(sig,annot,by=0), 
                       row.names=Row.names, Row.names=NULL)
# note often this ^^ results in a reduced number of sig genes in our data frame.
# but not necessarily because we want it to... why did it do this?

## Now choose only those genes with a particular search or GO term
## Resource: https://www.statology.org/r-partial-string-match/
sig.odor <- sig.annot[grep("odor", sig.annot$Description),]
## you could select a second term
sig.chem <- sig.annot[grep("chemo", sig.annot$Description),]
## then add the data frames together with those two terms
sig.chem <-rbind(sig.odor,sig.chem)

## Normalized counts from 
matrix.chem <-counts(dds.gen.de, normalized = T)[rownames(sig.chem),]

## get zscore for each row
matrix.z <-t(apply(matrix.chem,1,scale))

## now take sample IDs from your metadata and make
## them the column names
colnames(matrix.z) <- meta$sample

## Ok let's make a heatmap!!
Heatmap(matrix.z, cluster_rows=T, cluster_columns=T, column_labels=colnames(matrix.z),
        name="z score") 

#################### END example heat map of genes with ########################
###################### "odor" in their description #########################

################ END Heatmaps of DESeq2 data ###########




################ GO Enrichment ################

## Genes are be assigned to Gene Ontology terms (GO terms) based on shared function ##
## You can find more information on the geneontology database ##

## If a GO term is enriched in a data set, it means it is seen more frequently in 
## the control or background data set ##

## In our case the control or background data set are all of the genes that reads mapped
## to, regardless of their p-value. So that's our results.annotated file!

## Resources:
## https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
## https://archetypalecology.wordpress.com/2021/01/27/how-to-perform-kegg-and-go-enrichment-analysis-of-non-model-species-using-r/
## https://rstats101.com/separate-a-collapsed-column-into-multiple-rows/
## https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/

## require these two packages for the analysis
require(DOSE)
require (clusterProfiler)


## As a first example, let's just look at enrichment terms across all differentially
## expressed genes (this will include male and female, up or down)

## create a dataframe limited to significant targets
results.sig <-subset(results.annotated,results.annotated$padj <0.05)

## get a list of differential expressed gene id names
## from the significant results table
deg.genes <-results.sig$geneID

## Now we need to get GO terms associated with those genes
## some genes have mulltiple GOs, so let's makes those show up 
## as separate rows

## define an object called "go" that's really just the results.annotated table
go <- results.annotated
## remove empty entries
go <- subset(go,go$GOs !="-")
## remove the "GO:" from in front of every GO term and replace with "GO"
go$GOs <- gsub("GO:","GO",as.character(go$GOs))
## now if a gene has multiple GO terms associated with it
## put each on a new row with that gene ID
go <- go %>% separate_rows(GOs)

## create a data frame of GO term and associated gene ids to use in enrichment analysis
go2gene <-go[, c("GOs","geneID")]

## run enrichment analysis (in the clusterProfiler package)
## this will compare the list of gene ids / GO terms that were differentially expressed
## with all of the gene ids / GO terms we detected
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



################ END GO Enrichment ################





################ Box plots of counts data ###########




