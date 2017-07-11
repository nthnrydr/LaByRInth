## Copyright 2015, Jason Vander Woude, Nathan Ryder

## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at

##   http://www.apache.org/licenses/LICENSE-2.0

## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.


##' Impute the missing genotype calls.
##'
##' Use the LB-Impute algorithm to impute the missing genotypes of the
##' variants in the specified vcf file then save the results to the
##' specified output file.
##' @title
##' @param vcfin input vcf file
##' @param vcfout output vcf file
##' @param method specified as one of impute, compare, or randomremove
##' @param parents character vector of parent names in vcf input file
##' @return path to saved vcf file
##' @author Jason Vander Woude
LBMethod <- function(vcfin, vcfout, method, parents, cores=1) {
    ## Check validity of arguments
    if (length(parents != 2)) {
        stop("Two parents must be specified")
    }
    
    if (method == "compare") {       #this is the algorithm for accuracy metrics. Requires a file with correct genotypes, erroneous or missing genotypes, and imputed genotypes.
        compare(args)
    }
    else if (method == "randomremove") {      #I don't think this is used anymore.
        remove(args)
    }
    else if (method == "impute") {      #this is the primary algorithm for imputation.
        impute(args)
    }
}


##' Impute the missing genotype calls.
##'
##' Use the LB-Impute algorithm to impute the missing genotypes of the
##' variants in the specified vcf file then save the results to the
##' specified output file.
##' @title 
##' @param vcfin name of the vcf file to be imputed
##' @param vcfout name of the imputed output vcf file
##' @param parentnames character vector of parent variant
##' @param err Default per variant error rate. Modifiable with -err flag. 
##' @param recomb Default distance until probability of transition to a different state (recombination event) becomes 0.5. Modifiable with -recombdist flag. 
##' @param markovorder Number of markers to "peer forward" when determining the probability of a transition. Actual trellis window is markovorder + 2 (0 marker, 1 marker, and then a number of markers equal to markerorder). Changeable with -window flag 
##' @param offspringimpute Can impute offspring or parents. Default is parents. -offspringimpute sets this to true. 
##' @param parentimpute Default imputation of parental genotypes. -parentimpute is left as a flag for clarity, though not necessary. 
##' @param resolveconflicts Resolving conflicts will, for a marker with differing calls on the forward and reverse path, use the state with the higher probability. -resolveconflicts sets this to true. 
##' @param startcol startcol is the first column of the VCF information to contain genotype calls for samples. 
##' @param assumebiallelic  This is the minimum number of samples that must have a successful state call at a given site for a parental imputation to be made. Can be changed with -minsamples flag. 
##' @param minsamples This is the minimum number of samples that must have a successful state call at a given site for a parental imputation to be made. Can be changed with -minsamples flag. 
##' @param minfraction The minimum fraction of genotype calls assigned to a parent for a given marker for that genotype to be assigned to the parent. Can be changed with -minfraction flag. 
##' @param minmarkers The minimum number of non-missing genotypes on a chromosome for imputation to take place for a given sample. If the sample does not have a sufficient number of markers on a given chromosome, all markers on that chromosome will be set to missing. Will default to the markov order in later steps. Can be changed with -minmarkers flag. 
##' @param genotypeerror The expected genotyping error independent of coverage. For instance, misalignments or paralogous artifacts. Maximum emission probability of any state is 1 - genotyping error. Set with -genotypeerr flag. 
##' @param drp
##' @param keeporiginal
##' @return path to saved vcf file
##' @author Jason Vander Woude
LBImpute <- function(vcfin, vcfout, parentnamess,
                     err = 0.05, recomb = 1000000, markovorder = 6,
                     offspringimpute = FALSE, parentimpute = TRUE,
                     resolveconflicts = FALSE, startcol = 10,
                     assumebiallelic = FALSE, minsamples = 5,
                     minfraction = 0.5, minmarkers = 0,
                     genotypeerror = 0.05, drp = TRUE,
                     keeporiginal = FALSE) {
                                        #TODO: add mclapply option depending on installation#
                                        #TODO: make splitstr as strsplit[[1]]#
    
    startTime <- Sys.time()
                                        #TODO: look at error handling for file read#
    variants <- readLines(vcfin)
    isComment <- sapply(variants, function(line){substr(line,1,1) == "#"})
    ## getNameReturn <- getName.GetVariantNames(variants, NULL, startcol)
    ## startline <- getNameReturn[[1]]
    headerLines <- variants[isComment]     #remove the header, but save so it can be restored later
    variants <- variants[!isComment]
    variants <- do.call(rbind, lapply(variants, function(line){strsplit(line, "\t")[[1]]})) #make table

    header <- headerLines[length(headerLines)] #get column headings
    header <- substr(header, 2, nchar(header)) #remove leading '#'
    colnames(variants) <- toupper(strsplit(header, "\t")[[1]]) #[[1]] is because strsplit returns a 1-element list
    
    ## formatExample <- strsplit(variants[1], "\t")[9] #In a vcf file column 9 is the FORMAT column with colon seperated values
    formatExample <- variants[1, "FORMAT"]
    formatFields <- strsplit(formatExample, ":")[[1]]

    ## Determine which positions in the format and data fields mean what
    genotypecheck <- match("GT", formatFields) #GT = genotype
    readcheck <- match("AD", formatFields)     #AD = allele depth
    genotypecheck <- match("DP", formatFields) #DP = depth
    readquality <- match("GQ", formatFields)   #GQ = genotype read quality

    states <- 3 #homozygous parent 1, homozygous parent 2, heterozygous

    ## Group sites to impute by chromosome
    chroms <- unique(variants[,"CHROM"])
    imputationSets <- lapply(chroms, function(chrom) {
        variants[variants[,"CHROM"] == chrom,]
    })

    if (parentimpute) {
       newVariants <- lapply(imputationSets, imputeparents)
    } else if (offspringimpute) {
        newVariants <- lapply(imputationSets, imputeoffspring)
    }

    ## Write imputations to 
    sink(vcfout)
    writeLines(paste(headerLines, collapse="\n"))
    write.table(newVariants, sep="\t", row.names=F, col.names=F, quote=F)
    sink()
    
    runtime <- as.numeric(Sys.time() - startTime)
    writeLines(paste("Runtime:", runtime, "sec"))
}
