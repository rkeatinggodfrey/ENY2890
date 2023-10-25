## Differential Gene Expression (DGE)  

Now that we have (1) gene counts for each set of reads for each sample and (2) a metadata file that gives us biological conditions (sex, body part) for each sample, we can use statistical analyses to ask if particular genes were detected more frequently in certain samples or in certain biological conditions. 
  

### DGE DESeq2  

To analyze DGE, will use a program written in the R language called DESeq2.  

If you want to run this analysis locally on your own computer, you will need to download and install R and the user interface, RStudio.  

+ You can download R here: https://cran.rstudio.com/
+ You can download R studio here: https://rstudio.com/products/rstudio/download/#download

+ [Here are some helpful step-by-step instructions with screenshots from Colorado State University](https://www.stat.colostate.edu/~jah/talks_public_html/isec2020/installRStudio.html) 

This is a useful guide for working in RStudio
+ [Hands-On Programming with R github](https://rstudio-education.github.io/hopr/starting.html)
