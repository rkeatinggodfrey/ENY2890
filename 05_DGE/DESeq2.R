## Project: Fall 2023 CURE Course 
## Authors: R. Keating Godfrey
## Last updated: 23-11-08

## Resources:
## https://bioconductor.org/packages/release/bioc/html/DESeq2.html
## https://people.duke.edu/~ccc14/duke-hts-2017/Statistics/08032017/DESeq2-Notebook-introduction.html
## https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

library(DESeq2)
library(dplyr)
library(ggplot2)

install.packages("htmltools")
library(htmltools)



### (1) Read in counts data ###

## OPTION A: read in single file with all counts 
## (compiled using Google Sheets or Excel)

## set working directory to source file location
setwd("~/Documents/2023/03_CURE/DESeq2/DGE_cure")

## read in the already-concatenated counts file
counts <- read.csv("counts_example.csv", header=T)


## OPTION B: read in counts from folder and combine into one dataframe ###

## change working director to location of counts files
setwd("/Users/rkeatinggodfrey/Documents/2023/03_CURE/DESeq2/DGE_cure/counts")

## create a list object that contains all of the .tsv files in this directory
filedir <-list.files(pattern = ".tsv$")

## Import all files as a list of dataframes
all_files <- lapply(filedir,function(x) {
  read.table(file = x, 
             sep = '\t', 
             header = F,
             row.names=1)
})

## Combine into single dataframe
counts <-bind_cols(all_files)

## create a string of file names
file.names<-c(filedir[1], filedir)

## why does it duplicate the first file name? Remove it!
file.names <- file.names[2:26]

## now get rid of the ".tsv" in the file name
file.names <- gsub(".tsv", "", file.names)

## name the columns of your counts dataframe to match your metadata file
colnames(counts)<-c(file.names)

## we used the geneIDs as row names for importing and combining files,
## but now we need them to be a column for our analysis
## so we need to turn the gene id row names into a column (variable)
counts <-tibble::rownames_to_column(counts,"gene_id")



### (2) Import metadata file ###

## change directory to the one where your metadata is located
## unless you're already in that directory
setwd("/Users/rkeatinggodfrey/Documents/2023/03_CURE/DESeq2")

## Import metadata (colData) ###
metadata <- read.csv("/Users/rkeatinggodfrey/Documents/2023/03_CURE/DESeq2/DGE_cure/metadata_example.csv")

## **ONLY RUN If you imported the directory to create the counts dataframe,
## add file names to metadata**
metadata$sample.name <-file.names

## create a new data frame that contains only the columns you need
## for the example metadata this is in columns 1, 2, and 3.
meta <- metadata[,c(1,2,3)]
## you can rename the columns to match the rest of the script
colnames(meta) <-c("sample","sex","body_part")

## make sure sex and body part are factors (categorical data)
## for more about factors and data types in R, see:
## https://swcarpentry.github.io/r-novice-inflammation/12-supp-factors.html
meta$sex <-as.factor(meta$sex)
meta$body_part <-as.factor(meta$body_part)

## check that they are factors
is.factor(meta$sex)
is.factor(meta$body_part)


### (3) Set up DESeq2 Model and Run analysis ###
  
### Construct a DESeq data set object ###

## In this example I am comparing genitalia from males and 
## females, so I just need to designate comparisons by sex (~sex)
dds.gen <-DESeqDataSetFromMatrix(countData=counts, # the counts data object
                             colData=meta, # the meta data
                             design=~sex, # the statistical design 
                             tidy=T)
### Run DESeq analyis ###
dds.gen <- DESeq(dds.gen)

## Get the results 
results.gen <-results(dds.gen)

## Look at summary of analysis
summary(results.gen)

## Look at the top of the results file
head(results.gen(dds,tidy=T))

## Look at the top of the results file ordered by adjusted p-value
results.gen <- results.gen[order(results.gen$padj),]
head(results.gen)


## Get results based on two factor levels
results.mvf <-results(dds.gen, contrast=c("sex","Male","Female"))

## Look at summary of analysis 
## (this will be the same as when you didn't use "contrast")
summary(results.mvf)

## Look at the top of the results file ordered by adjusted p-value
results.mvf <- results.mvf[order(results.mvf$padj),]
head(results.mvf)


## (4) Merge results with functional annotation

## The geneIDs are the row names in the original results table
## turn them into a column to merge with functional annotation from eggnog 
library(tibble)
results.table <- as.data.frame(results.mvf)
results.table <- tibble::rownames_to_column(results.table,"geneID")

## read in the annotation file
annot <- read.csv("rna.annot.cds.emapper.annotations.geneID.csv", header=T)

## merge DGE with annotation file
results.annotated <-merge(results.table,annot, by=c("geneID"))

## you can write this out as a CSV
write.csv(results.annotated,"DGE_annotated_test_female_male_geni.csv")

