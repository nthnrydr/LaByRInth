## RUN ONCE: source("http://bioconductor.org/biocLite.R") 
biocLite("VariantAnnotation") #install the package
library("VariantAnnotation") #load the package
## example("readVcf") #optional, test the function by running example codes

a <- c(1,2,0,NA,2,2,1)
b <- c(NA,2,0,NA,1,2,1)

ifelse(a==b, a, NA)

NA[c(1,3,7)]

this is a line to see if i can get to 80 characters and hope that that
will enable
