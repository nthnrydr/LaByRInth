### prefs
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
###   $states


str.to.num <- function(str, sep) {
    as.numeric(str.split(str, sep))
}


## More convenient strsplit if length of vector is 1
str.split <- function(str, sep) {
    if (length(str) != 1) {
        warning("Only splitting the first element of the vector")
    }
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


##' Extract the information from a vcf file and save it as a vcf object
##'
##' The returned vcf object will have the following variants, header.lines,
##'     variant.names, chrom.names, GT, AD.
##' @title 
##' @param file the path to the vcf file
##' @return the vcf object created from the file
##' @author Jason Vander Woude
VCF <- function(file) {
    ## TODO(Jason): add filtering step to remove non-biallelic calls or talk to
    ## Jesse about how those should be handled
    
    vcf <- list()
    class(vcf) <- "vcf"

    writeLines("Reading vcf file")
    vcf$variants <- readLines(file)
    isComment <- sapply(vcf$variants, function(line){substr(line,1,1) == "#"})

    vcf$header.lines <- vcf$variants[isComment]     #remove the header, but save so it can be restored later
    vcf$variants <- vcf$variants[!isComment]
    vcf$variants <- do.call(rbind, lapply(vcf$variants, function(line){str.split(line, "\t")}))  #make table

    header <- vcf$header.lines[length(vcf$header.lines)]  #get column heading
    header <- substr(header, 2, nchar(header)) #remove leading '#'
    colnames(vcf$variants) <- strsplit(header, "\t")[[1]] #[[1]] is because strsplit returns a 1-element list
    
    formatExample <- vcf$variants[1, "FORMAT"]
    field.names <- str.split(formatExample, ":")

    ## Verify that required fields are in the VCF
    required.fields <- c("GT", "AD")  # could also use GQ
    if (!all(required.fields %in% field.names)) {
        stop(paste("VCF file does not contain all required fields.",
                   "Required fields are", toString(required.fields)))
    }
    field.indices <- match(required.fields, field.names)
    names(field.indices) <- required.fields

    ## The "FORMAT" column is the last one before the variants start
    format.col <- match("FORMAT", colnames(vcf$variants))
    n.variants <- ncol(vcf$variants) - format.col
    n.sites <- nrow(vcf$variants)

    samples <- vcf$variants[ , (format.col + 1):ncol(vcf$variants)]
    rownames(samples) <- paste0(vcf$variants[, "CHROM"], ":", vcf$variants[, "POS"])
    
    vcf$variant.names <- colnames(samples)
    vcf$chrom.names <- unique(vcf$variants[, "CHROM"])



    
    ## TODO(Jason): make sure that errors don't happend with '.' instead of
    ## '3,4' or '1/2/3' the line "flat.mat <- ..." should be changed


    ## GT section
    writeLines("Converting genoytpe data")
    mat <- apply(samples, 1:2, function(sample) {
        str.split(sample, ":")[field.indices["GT"]]
    })
    ## Replace '|' with '/' in genotype matrix
    mat <- gsub("\\|", "/", mat)
    ## introduce NA by conversion of .
    flat.mat <- suppressWarnings(as.numeric(unlist(strsplit(mat, "/"))))
    third.dim <- length(flat.mat) / n.sites / n.variants
    if (third.dim%%1 != 0) {
        stop("Error reading genotypes. One of the sites is likely not biallelic.")
    }
    vcf$GT <- array(reorder(flat.mat, third.dim), dim=c(n.sites,
                                                        n.variants, third.dim))
    colnames(vcf$GT) <- colnames(samples)
    rownames(vcf$GT) <- rownames(samples)

    
    ## AD section
    writeLines("Converting allelic depth data")
    mat <- apply(samples, 1:2, function(sample) {
        str.split(sample, ":")[field.indices["AD"]]
    })
    ## browser()
    ## Code for biallelic check during browser() it is bad if it prints
    ## something
    ## dumby <- apply(mat, 1:2, function(x){if(length(str.split(x, ",")) != 2) {print(x)}})
    
    ## introduce NA by conversion of .
    flat.mat <- suppressWarnings(as.numeric(unlist(strsplit(mat, ","))))
    third.dim <- length(flat.mat) / n.sites / n.variants
    if (third.dim%%1 != 0) {
        stop("Error reading allelic depths. One of the sites is likely not biallelic.")
    }
    vcf$AD <- array(reorder(flat.mat, third.dim), dim=c(n.sites,
                                                        n.variants, third.dim))
    colnames(vcf$AD) <- colnames(samples)
    rownames(vcf$AD) <- rownames(samples)
    
    vcf  # implicit return
}


##' Get a subset of the data from the vcf object
##'
##' Get data pertaining to the speficied field and subset it by samples and
##'     chromosomes.
##' @title 
##' @param vcf an object of class vcf
##' @param field the name of a single field ("GT" or "AD")
##' @param samples a vector of indices or names of variants
##' @param chromosomes a vector of chromosome names to subset by
##' @return a 3-dimensional array representing the subset of the data
##' @author Jason Vander Woude
Get <- function(vcf, field, samples, chromosomes=NULL) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    if (length(field) != 1) {
        stop("Lenght of field must be 1")
    } 
    if (is.null(chromosomes)) {
        chromosomes <- vcf$chrom.names
    } 
    rows <- vcf$variants[, "CHROM"] %in% chromosomes
    vcf[[field]][rows, samples, , drop=F]
}


##' Resolve heterozygous sites in the samples
##'
##' Returns a matrix of genotype calls for the samples such that the entry is 0
##'     if all calls were for the reference allele; 1 if all calls were for the
##'     alternate allele; and NA if there were calls for both the reference and
##'     alternate allele.
##' @title
##' @param vcf an object of class vcf
##' @param samples a vector of indices or names of variants
##' @param prefs a preferences object
##' @return a matrix of genotypes (0, 1, NA)
##' @author Jason Vander Woude
ResolveHomozygotes <- function(vcf, samples) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    genotype <- Get(vcf, "GT", samples)
    allele.counts <- Get(vcf, "AD", samples, vcf$chrom.names)
    ret.val <- genotype[, , 1]  # initilize to first slice in 3rd dim
    for (r in 1:nrow(genotype)) {
        for (c in 1:ncol(genotype)) {
            alleles <- genotype[r, c, ]
            
            if (any(is.na(alleles))) {  # One of the alleles is NA
                ret.val[r, c] <- NA
            } else if (all(alleles == alleles[1])) {  # Check if all are same)
                ret.val[r, c] <- alleles[1]
            } else if (sum(allele.counts[r, c, ] != 0) == 1) {  # Only 1 counted
                ret.val[r, c] <- alleles[as.logical(allele.counts[r, c, ])]
            } else {  # Contradictory calls
                ret.val[r, c] <- NA
            }
        }
    }
    ret.val  # implicit return
}


##' Get matrix of emission probabilities
##'
##' Computes the emission probabilities for each state based on allelic depth of
##'     coverage (using the binomial assumption). Sample must be of length 1 and
##'     the entries of the parent.geno matrix must be either 0, 1, or NA.
##' @title
##' @param vcf an object of class vcf
##' @param sample an index or name of a variant
##' @param parent.geno a matrix of parental genotypes
##' @param prefs a preferences object
##' @return a matrix of posterior probabilities
##' @author Jason Vander Woude
GetProbabilities <- function(vcf, sample, parent.geno, prefs) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    if (length(sample) != 1) {
        stop("Length of sample must be 1")
    }

    gt <- Get(vcf, "GT", sample)
    ad <- Get(vcf, "AD", sample)

    ret.val <- matrix(NA, nrow = nrow(gt), ncol = prefs$states)

    for (row in 1:nrow(ret.val)) {
           
        ## the '1' in [i, 1, ] is because there is only 1 sample
        geno.calls <- gt[row, 1, ]
        allele.counts <- ad[row, 1, ]

        if (all(is.na(geno.calls))) {
            ret.val[row, ] <- rep(1, prefs$states)
        } else {
            ## TODO(Jason): check this against the strange entries noted in
            ## the notes file

            ## 0 means reference by the vcf standards
            ref.indices <- which(geno.calls == 0)
            ## 1 means first alternate by the vcf standards
            alt.indices <- which(geno.calls == 1)

            ref.calls <- sum(allele.counts[ref.indices])
            alt.calls <- sum(allele.counts[alt.indices])

            ## Due some to strange entries in the vcf files it is possible
            ## that both genotypes are the same number thus one of the
            ## alleles will have no indices and summing over those NA
            ## indices will yield NA's
            if (is.na(ref.calls)) {
                ref.calls <- 0
            }
            if (is.na(alt.calls)) {
                alt.calls <- 0
            }

            max.allowed <- 1 - (2 * prefs$genotype.err)
            min.allowed <- prefs$genotype.err
            
            ## Calculate the emission probabilities for this site
            rerr <- prefs$read.err
            ref.prob <- (1 - rerr)**ref.calls * (rerr)**alt.calls
            alt.prob <- (1 - rerr)**alt.calls * (rerr)**ref.calls
            hom.prob <- (0.5)**(ref.calls + alt.calls)  # homozygous

            max.prob <- max(ref.prob, alt.prob, hom.prob)

            normalize <- function(x) {
                ## TODO(Jason): Correction: max.allowed should be swapped with
                ## (max.allowed - min.allowed)
                x / max.prob * max.allowed + min.allowed
            }

            ## Assumes biallelic TODO(Jason): make more general by creating
            ## alt.prob on the fly for each alternate if the site is not biallelic?
            for (state in 1:(prefs$states - 1)) {
                if (is.na(parent.geno[row, state])) {
                    ret.val[row, state] <- normalize(max(alt.prob, ref.prob))
                } else if (parent.geno[row, state] == 0) {
                    ret.val[row, state] <- normalize(ref.prob)
                } else if (parent.geno[row, state] == 1) {
                    ret.val[row, state] <- normalize(alt.prob)
                } else {
                    stop("Parental genotype was not NA, 0, or 1")
                }
            }
            ret.val[row, prefs$states] <- normalize(hom.prob)
        }
    }
    ret.val  # implicit return
}
