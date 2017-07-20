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



str.to.num <- function(str, sep) {
    as.numeric(str.split(str, sep))
}


str.split <- function(str, sep) {
    strsplit(str, sep)[[1]]
}


nsites.probs <- function(probs) {
    nrow(probs)
}


nstates.probs <- function(probs) {
    ncol(probs)
}

## Used for creating a 3-D structure
reorder <- function(x, n) {
    as.vector(sapply(1:n, function(i) {
        v <- rep(FALSE, n)
        v[i] <- TRUE
        x[v]
    }))
}


##' Generate a path from a path tracker
##'
##' Given a path tracker (which is a specific type of matrix used in the viterbi
##'     algorithm) and an index representing the final hidden state of the path,
##'     compute the path that ended in that state.
##' @title 
##' @param path.tracker the integer matrix from which the path comes
##' @param index which row does the path end at
##' @return the path
##' @examples
##' path.tracker.a <- matrix(c(3,1,1,1,1,1,3,2,2,2,3,3),byrow=T,nrow=3)
##' path.tracker.a
##' generatePath(path.tracker.a, 1)
##' generatePath(path.tracker.a, 2)
##' generatePath(path.tracker.a, 3)
##' @author Jason Vander Woude
generatePath <- function(path.tracker, index) {
    path <- rep(NA, ncol(path.tracker))
    for (i in length(path):1) {
        path[i] <- index <- path.tracker[index, i]
    }
    path  # return path
}


##' Find the most probable path using the viterbi algorithm
##'
##' See http://homepages.ulb.ac.be/~dgonze/TEACHING/viterbi.pdf for details
##'     about the viterbi algorithm and understanding this implementation
##' @title 
##' @param probs emission probabilities
##' @param dists vector of distances between sites
##' @param prefs imputation preferences object
##' @return the most probable sequence of states
##' @author Jason Vander Woude
viterbi <- function(probs, dists, prefs) {
    nstates <- nstates(probs)
    path.size <- nsites(probs)
    
    paths.tracker <- matrix(NA, nrow=nstates, ncol=path.size)
    probs.tracker <- probs[j, ]  # initialize probabilities

    ## hard code the first column to the vector 1,2,...,nstates
    ## this is what the generatePath function will need
    paths.tracker[, 1] <- 1:nstates
    
    for (site in 2:path.size) {  # for each site in the path
        dist <- dists[i - 1]
        ## for each possible hidden state at this site
        probs.tracker <- sapply(1:nstates, function(state) {
            extension.probs <- sapply(1:nstates, function(i) {
                ## Probability of being at state i before and transitioning to
                ## the 'state' state
                probs.tracker[i] * transProb(i, state, dist, prefs)
            })
            ## which partial path has the highest probability of moving to state
            ## 'state'
            paths.tracker[state, site] <- which.max(extension.probs)
            max(extension.probs) * probs[site, state]  # return new probability
        })
    }
    generatePath(path.tracker, which.max(probs.tracker))  # return best path
}


##' Find the transission probability between hidden states
##'
##' Use the equations specified in the LB-Impute paper to compute the
##'     transmission probability between two sites
##' @title 
##' @param a first site
##' @param b second site
##' @param dist the distance between sites a and b
##' @param prefs a preferences object
##' @return the transmission probability between these hidden states
##' @author Jason Vander Woude
transProb <- function(a, b, dist, prefs) {
    if (a == b) {
        ## If there is no recombination
        0.5 * (1 + exp(-dist/prefs$recomb.dist))
    } else if (a %in% 1:2 && b %in% 1:2 && !prefs$recomb.double) {
        ## Double recomb occurred and has square of probability of single recomb
        ## This only works for biparental
        (0.5 * (1 - exp(-dist/prefs$recomb.dist)))**2
    } else {
        ## Some type of recombination occurred that has regular probability
        0.5 * (1 - exp(-dist/prefs$recomb.dist))
    }
}


##' Increment a vector as if it were the digits of a number
##'
##' Given a vector, try to increment the last index by 1 unless that would cause
##'     it to be greater than max in which case the last index will be set to
##'     min and the previous index incremented in this same manner.
##' @title 
##' @param x the vector
##' @param max the maximum value that any element should have
##' @param min the minimum value that any element should have
##' @return the incremented vector
##' @author Jason Vander Woude
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


VCF <- function(file, full=T) {
    vcf <- list()
    class(vcf) <- "vcf"
    
    vcf$variants <- readLines(file)
    isComment <- sapply(variants, function(line){substr(line,1,1) == "#"})

    vcf$headerLines <- variants[isComment]     #remove the header, but save so it can be restored later
    vcf$variants <- vcf$variants[!isComment]
    vcf$variants <- do.call(rbind, lapply(vcf$variants, function(line){str.split(line, "\t")}))  #make table

    header <- vcf$headerLines[length(vcf$headerLines)]  #get column heading
    header <- substr(header, 2, nchar(header)) #remove leading '#'
    colnames(vcf$variants) <- strsplit(header, "\t")[[1]] #[[1]] is because strsplit returns a 1-element list
    
    formatExample <- vcf$variants[1, "FORMAT"]
    field.names <- str.split(formatExample, ":")

    ## The "FORMAT" column is the last one before the variants start
    n.variants <- ncol(vcf$variants) - match("FORMAT", colnames(vcf$variants))
    n.sites <- nrow(vcf$variants)
    
    ## TODO(Jason): make this more efficient
    vcf$fields <- lapply(field.names, function(field.ID) {
        ## Inside here I have a field string and want to compute the 3-D array
        ## corresponding to that field
        
        ## how many in third dimension? try doing with predetermining
        arr <- array(NA, dim=c(n.sites, n.variants, 
    })
}


Get.vcf <- function(vcf, samples, chromosomes, field) {
    
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


GetProbabilities <- function(variants, sample, prefs) {
    ## In a vcf file column 9 is the FORMAT column with colon separated values
    ## Retrieve an example from the FORMAT column to determine ordering
    format.example <- variants[1, "FORMAT"]
    ## [[1]] is because strsplit returns a 1-element list
    format.fields <- strsplit(format.example, ":")[[1]] 

    ## n.prefix.cols <- match("FORMAT", colnames(variants))
    ## n.samples <- ncol(variants) - n.prefix.cols
    
    ## Determine which positions in the format and data fields mean what
    genotype.index     <- match("GT", format.fields)  #GT = genotype
    allele.count.index <- match("AD", format.fields)  #AD = allele count
    depth.index        <- match("DP", format.fields)  #DP = depth
    readqual.index     <- match("GQ", format.fields)  #GQ = genotype read quality

    states <- 3  # homozygous parent 1, homozygous parent 2, heterozygous

    sample.variant <- variants[, sample]  # just the column associated

    ret.val <- matrix(NA, nrow = nrow(variants), ncol = states)

    ## From a single entry 
    GetProb <- function(entry) {
        ret.val <- NA
        if (entry == "./.") {
            ret.val <- rep(1, states)
        } else {
            info <- strsplit(entry, ":")[[1]]
            if (info[genotype.index] == "./.") {
                ret.val <- ret(1, states)
            } else {
                ## TODO(Jason): check this against the strange entries noted in
                ## the notes file
                genotype <- str.to.num(info[genotype.index], ",")
                allele.count <- str.to.num(info[allele.count.index], ",")
                ## 0 means reference by the vcf standards
                ref.allele.indices <- which(allele.count == 0)
                ## 1 means first alternate by the vcf standards
                alt.allele.indices <- which(allele.count == 1)

                ref.calls <- sum(allele.count[ref.allele.indices])
                alt.calls <- sum(allele.count[alt.allele.indices])

                ## Due so some strange entries in the vcf files it is possible
                ## that both genotypes are the same number thus one of the
                ## alleles will have no indices and summing over those NA
                ## indices will yield NA's
                if (is.na(ref.calls)) {
                    ref.calls <- 0
                }
                if (is.na(alt.calls)) {
                    alt.calls <- 0
                }

                rerr <- prefs$read.err
                max.allowed <- 1 - (2 * prefs$genotype.err)
                min.allowed <- prefs$genotype.err
                
                ## Calculate the emission probabilities for this site
                ref.prob <- (1 - rerr)**ref.calls * (rerr)**alt.calls
                alt.prob <- (1 - rerr)**alt.calls * (rerr)**ref.calls
                hom.prob <- (0.5)**(ref.calls + alt.calls)  # homozygous

                max.prob <- max(ref.prob, alt.prob, hom.prob)

                normalize <- function(x) {
                    x / max.prob * max.allowed + min.allowed
                }
                
                ## TODO(Jason): Correction: max.allowed should be swapped with
                ## (max.allowed - min.allowed)
                for (state in 1:(states - 1)) {
                    if (parent.map[k, state] == 0) {
                        ret.val[state] = normalize(ref.prob)
                    } else if (parent.map[k, state] == 1) {
                        ret.val[state] = normalize(alt.prob)
                    } else {
                        ret.val[state] = normalize(max(alt.prob, ref.prob))
                    }
                }
                ## TODO(Jason): code not finished here
            }
        }
        ret.val  # implicit return
    }
    
    ## prob.path <- matrix(nrow = ncol(just.variants), ncol = states)  # initialized probability path matrix
    
    ## Begin the 
    apply(, 2, function(variant) {
              sapply(variant, function(marker) {
                         marker.fields <- strsplit(marker, ":")[[1]]
                         if(marker.fields[genotype.check]=="./.") {
                             prob.path[]
                     }})
          })
          
    GenotExists(offspring$Genotype[, 1])

}

