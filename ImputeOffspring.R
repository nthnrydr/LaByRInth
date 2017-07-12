## Copyright 2017, Jason Vander Woude, Nathan Ryder

## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at

##   http://www.apache.org/licenses/LICENSE-2.0

## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.

## This work is derived from the original LB-Impute method written by
## Christopher Heffelfinger, Christopher Fragoso, Hongyu Zhao, and
## Stephen Dellaporta. The method was originally written in Java and
## the purpose of this project was to port it to R. The original work
## was also licensed under APLv2 and can be found here:
## https://github.com/dellaporta-laboratory/LB-Impute



imputeoffspring <- function(mystates, mystartcol, varcount, samplecount,
                            mygenotypecheck, myreadcheck, mycoveragecheck, myreadqualcheck,
                            myerr, myrecomb, mymarkovorder, myresolve, myminmarkers,
                            mygenotypeerror, mydrp, mykeeporiginal) {

                                        #TODO: not all indices of the colon seperation are guaranteed to exist. watch for this#
    
    
}

##' Find parental genotype calls.
##'
##' Return a matrix with a row per position and a column per parent
##' plus an additional column for heterozygous calls. Each entry is an
##' integer representing the genotype call for that parent according
##' to the vcf specifications or NA if no call was made.
##' @title
##' @param variants matrix of variants
##' @param parents vector of parent names corresponding to column names of variants
##' @return a numeric matrix representing genotype
##' @author Jason Vander Woude
retrieveparentals <- function(variants, parents) {
    ## Matrix with a column for each parent
    retval <- variants[,parents]

    ## Turn into a list of splits
    listOfEntries <- strsplit(retval, split=":")

    ## This is blindly taking the first call before the '/'
    ## TODO: fix this
    ## TODO: verify that the coverage depth is not 0
    listOfEntries <- lapply(listOfEntries, function(entry) {substr(entry[genotypecheck], 1, 1)})

    ## Convert back to a vector and let as.numeric coerce '.' into 'NA'
    ## Turn the vector back into the properly sized matrix
    retval <- matrix(as.numeric(unlist(listOfEntries)), ncol=length(parents))
    
    ## Bind a column for heterozygous information
    retval <- cbind(retval, NA)
}


