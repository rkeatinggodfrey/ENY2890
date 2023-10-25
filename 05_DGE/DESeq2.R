## Project: Fall 2023 CURE Course 
## Authors: R. Keating Godfrey
## Last updated: 23-10-18

## Resources:
## https://bioconductor.org/packages/release/bioc/html/DESeq2.html
## https://people.duke.edu/~ccc14/duke-hts-2017/Statistics/08032017/DESeq2-Notebook-introduction.html
## https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/DESeqDataSet-class


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
library(ggplot2)

install.packages("htmltools")
library(htmltools)
