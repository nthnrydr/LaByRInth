### impute.prefs
###   $markov.order
###   $parents
###   $resolve.conflicts
###   $read.err
###   $genotype.err
###   $recomb.err
###   $recomb.dist
###   $recomb.double
###   $window.size
###   $min.samples
###   $min.fraction
###   $min.markers


nprobs.probs <- function(probs) {
    nrow(probs)
}

nstates.probs <- function(probs) {
    ncol(probs)
}


FindPath <- function(variants, probs, prefs) {
    ## probs should be a probs object
    nprobs <- nprobs(probs)
    nstates <- nstates(probs)
    trellis <- prefs$markov.order + 2  # viterbi trellis size
    
    ret.val <- rep(NA, nprobs)          # TODO(Jason): check if this is correct size
    class(ret.val) <- "path"

    if (nprobs < min.markers) {  # too few markers are present
        return ret.val  # return a path of NA's
    } else {  # there are sufficient markers for imputation

        for (i in 2:nprobs) {  
            current.trellis <- min(trellis, )  # TODO(Jason): look into this
            current.path <- rep(nstates, current.trellis - 1)  # last path so inc to first

            best.path <- current.path  # keep track of best path so far
            best.prob <- 0  # keep track of highest probability
            
            while (TRUE) {
                current.path <- increment(current.path)  # next path to check
                last.value <- ret.val[i-1]  # previous imputed site

                if (all(current.path == nstates)) {  # this is the last
                    break  # exit the while loop
                }
            }            
        }
    }
}


TransProb <- function(recomb, dists, prefs) {
    if (recomb) {
        0.5 * (1 - exp(-dist/prefs$recomb.dist))
    } else {
        0.5 * (1 + exp(-dist/prefs$recomb.dist))
    }
}


Increment <- function(x, max, min = 1) {
    if (min >= max) {
        stop("max must be greater than min")
    }
    if (any(x > max)) {
        stop("x has elements greater than max")
    }
    if (any(x < min)) {
        stop("x has elements less than min")
    }

    i <- length(x)

    if (all(x == max)) {  # all elements of x are equal to max
        x <- rep(min, i)  # return a list of min's
    } else if (x[i] < max) {  # typical case
        x[i] <- x[i] + 1  # increase the last element
    } else {  # last element is max
        x <- c(increment(x[-i], max, min), min)  # 
    }
    x
}


## Rewrite of "public void getprobabilities2(List<String> variants, int sample)"
## from ImputeOffspring.java.  This calculates the emission probabilities for
## each state based on allelic depth of coverage (using binomial assumption).

## Incoming Variables:

## formatExample <- strsplit(variants[1], "\t")[9]


## offspring <- a 
## num.variants <- the amount of variants (offspring).  Also length of vcf rows
## minus identifiers
## num.markers <- the amount of markers in the chromosome of a variant.  Also
## length of vcf columns
## genot.exists <- boolean vector of length num.markers of whether the genotype call (?/?)
## is empty (./.) and if it is empty, genotypeboolean = F.
##


GetProbabilities <- function(variants) {
    ## In a vcf file column 9 is the FORMAT column with colon separated values
    format.example <- variants[1, "FORMAT"]  # Retrieve an example from the FORMAT column
    format.fields <- strsplit(format.example, ":")[[1]] #[[1]] is because strsplit returns a 1-element list

    ## Determine which positions in the format and data fields mean what
    genotype.check <- match("GT", format.fields) #GT = genotype
    allele.read.check <- match("AD", format.fields)  #AD = allele depth
    depth.check <- match("DP", format.fields)  #DP = depth
    readqual.check <- match("GQ", format.fields)   #GQ = genotype read quality

    states <- 3  # homozygous parent 1, homozygous parent 2, heterozygous
    just.variants <- variants[, (match("FORMAT", header) + 1):ncol(variants)]  # just the columns of variants
    
    prob.path <- matrix(nrow = ncol(just.variants), ncol = states)  # initialized probability path matrix
    
    ## Begin the 
    apply(, 2, function(variant) {
              sapply(variant, function(marker) {
                         marker.fields <- strsplit(marker, ":")[[1]]
                         if(marker.fields[genotype.check]=="./.") {
                             prob.path[
                     })
          })
          
    GenotExists(offspring$Genotype[, 1])
