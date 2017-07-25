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


imputeoffspring <- function(variants, parents,
                            mystates, mystartcol, varcount, samplecount,
                            mygenotypecheck, myreadcheck, mycoveragecheck, myreadqualcheck,
                            myerr, myrecomb, mymarkovorder, myresolve, myminmarkers,
                            mygenotypeerror, mydrp, mykeeporiginal) {

                                        #TODO: not all indices of the colon seperation are guaranteed to exist. watch for this#
    
    parmap <- retrieveparentals(variants, parents)
    total <- ncol(variants)
                                        #TODO: don't just take first Allele Depth number#
                                        #TODO: use a nice progress bar if possible or a rotating wheel if in terminal#
    ## allsamplesmap <- matrix(NA, )
}



