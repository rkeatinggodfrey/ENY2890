## Differential Gene Expression (DGE)  

<p align="center">
<img width="400px" src="ENY2890/Images/Research_InfoGraphic.jpg">
</p> 

Now that we have 
+ (1) gene counts for each set of reads for each sample and 
+ (2) a metadata file that gives us biological conditions (sex, body part) for each sample,  

we can use statistical analyses to ask if particular genes were detected more frequently in certain samples or in certain biological conditions.  

The general structure of this question is:

y ~ sex * body part  

which asks if gene expression (*y*) is determined by (*~*) the interaction of the variables *sex* and *body part*.  
  

### (1) Download R and R studio to run DESeq2 

To analyze DGE, will use a program written in the R language called DESeq2.  

If you want to run this analysis locally on your own computer, you will need to download and install R and the user interface, RStudio.  

+ You can download R here: https://cran.rstudio.com/
+ You can download R studio here: https://rstudio.com/products/rstudio/download/#download

+ [Here are some helpful step-by-step instructions with screenshots from Colorado State University](https://www.stat.colostate.edu/~jah/talks_public_html/isec2020/installRStudio.html) 

This is a useful guide for working in RStudio
+ [Hands-On Programming with R github](https://rstudio-education.github.io/hopr/starting.html)


### (2) Download counts (.tsv) files from HiPer Gator  

At this point everyone needs a copy of all counts files. Download your counts files and share copies with your group members!  

Make a directory on your computer for differential gene analysis (e.g., "DGE_Analysis"). Make a folder within this directory just for your counts.tsv files (e.g., "counts_files"). Put all of your counts files within this folder. **Do not put anything else in this folder**  

Put a copy of your metadata file in .csv format in the parent differential gene expression directory. **Make sure names in the "sample" column of your metadata match the counts file names in your counts directory**

