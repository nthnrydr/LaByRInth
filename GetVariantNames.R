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


##' Get the names of all sample variants.
##'
##' Returns a list with 2 entries. The first is the line index of the first
##' chromosome and the second is the character vector of variant names.
##' @title
##' @param unfilteredvariants a list of strings represinting lines in a vcf file
##' @param skip deprecated
##' @param namestart index where
##' @return
##' @author Jason Vander Woude
getName.GetVariantNames <- function(unfilteredvariants, skip, namestart) {
    ## browser()
    f = 1                               #int
    startline = 1                       #int
    nameline = 1                        #int
    foundstart = FALSE                  #boolean
    while(!foundstart) {
        variantlineparts = strsplit(unfilteredvariants[f], split="\t")[[1]] #String[]
        if(variantlineparts[1] == "#CHROM") {
            startline <- f + 1
            nameline <- f
            foundstart <- TRUE
        }
        if(f > length(unfilteredvariants)){
            stop("Unable to find line beginning with #CHROM, which is the standard string in vcf files indicating line with variant names. Exiting")
        }
        f <- f+1
    }
    namelineparts <- strsplit(unfilteredvariants[nameline], "\t")[[1]] #String[]
    samplenames <- c()                  #List<String>
    for(i in seq(namestart, length(namelineparts))){
        samplenames <- c(samplenames, namelineparts[i])
    }
    list(startline, samplenames)
}

