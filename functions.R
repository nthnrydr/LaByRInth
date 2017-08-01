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
    nstates <- nstates.probs(probs)
    path.size <- nsites.probs(probs)

    paths.tracker <- matrix(NA_integer_, nrow=nstates, ncol=path.size)
    ## TODO(Jason): make sure this is right. Changed an invalid index j to 1 on
    ## July 26

    ## This will keep track of the overall probabilites of each of the {nstates}
    ## final paths, so that the probability of the final paths do not have to be
    ## computed again. The probabilities are initialized to the emission
    ## probabilities of the first site of the probs matrix which is the first
    ## row
    if (path.size < 1) {
        stop("viterbi requires that probs have > 0 rows")
    }
    probs.tracker <- log(probs[1, ])

    ## Hard code the first column to the vector 1,2,...,nstates as an index
    ## This is what the generatePath function will need
    paths.tracker[, 1] <- 1:nstates
    if (path.size != 1) {  # if the path size is 1 just use the emission probs
        for (site in 2:path.size) {  # for each site in the path
            
            dist <- dists[site - 1]
            
            ## log of the probability for each possible hidden state at this site
            probs.tracker <- sapply(1:nstates, function(state) {
                
                extension.probs <- sapply(1:nstates, function(i) {
                    ## Probability of being at state i before and transitioning to
                    ## the 'state' state
                    probs.tracker[i] + log(transProb(i, state, dist, prefs))
                })
                
                ## which partial path has the highest probability of moving to state
                ## 'state'
                paths.tracker[state, site] <<- which.max(extension.probs)
                max(extension.probs) + log(probs[site, state])  # return new probability
            })
        }
    }

    ## The code above has already computed the optimal path, but that
    ## information is encoded within the paths.tracker matrix and needs to be
    ## extracted. That is what generatePath will do when passed the path.tracker
    ## matrix and the index of the optimal path.
    generatePath(paths.tracker, which.max(probs.tracker))  # return best path
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

    ## Remove the header, but save so it can be restored later
    vcf$header.lines <- vcf$variants[isComment]     
    vcf$variants <- vcf$variants[!isComment]
    ## Make table
    vcf$variants <- do.call(rbind,
                            lapply(vcf$variants,
                                   function(line){str.split(line, "\t")}))

    header <- vcf$header.lines[length(vcf$header.lines)]  #get column heading
    header <- substr(header, 2, nchar(header)) #remove leading '#'
    colnames(vcf$variants) <- str.split(header, "\t")
    
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



    
    ## TODO(Jason): make sure that errors don't happen with '.' instead of
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
    
    ## Introduce NA by conversion of .
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
##' Get data pertaining to the specified field and subset it by samples and
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
        stop("Length of field must be 1")
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
    genotype <- Get(vcf, "GT", samples, vcf$chrom.names)
    allele.counts <- Get(vcf, "AD", samples, vcf$chrom.names) 
    ret.val <- genotype[, , 1]  # initilize to first slice in 3rd dim
    for (r in 1:nrow(genotype)) {  # by chromosome site
        for (c in 1:ncol(genotype)) {  # by variant/sample
            alleles <- genotype[r, c, ]  # grab the observed alleles from a specific site of a variant chromosome
            
            if (any(is.na(alleles))) {  # One of the alleles is NA
                ret.val[r, c] <- NA
            } else if (all(alleles == alleles[1])) {  # Check if all are same
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
GetProbabilities <- function(vcf, sample, chromosomes, parent.geno, prefs) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    if (length(sample) != 1) {
        stop("Length of sample must be 1")
    }
    
    gt <- Get(vcf, "GT", sample, chromosomes)
    ad <- Get(vcf, "AD", sample, chromosomes)

    ret.val <- matrix(NA_integer_, nrow = nrow(gt), ncol = prefs$states)
    rownames(ret.val) <- rownames(gt)
    colnames(ret.val) <- c(colnames(parent.geno), "HETEROZYGOUS")
    class(ret.val) <- "prob"
             
    for (row in 1:nrow(ret.val)) {
           
        ## the '1' in [i, 1, ] is because there is only 1 sample
        geno.calls <- gt[row, 1, ]
        allele.counts <- ad[row, 1, ]

        if (all(is.na(geno.calls))) {
            ret.val[row, ] <- rep(1, prefs$states)
        } else if (any(is.na(geno.calls))) {
            stop("Some but not all genotype calls are NA")
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

            ## probability constraints from LB-Impute original
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


## Determine which rows are real calls
GetRelevantProbabiltiesIndex <- function(vcf, chromosomes, parent.geno, prefs) {
    if (! inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    rows <- vcf$variants[, "CHROM"] %in% chromosomes
    apply(parent.geno[rows, , drop=F], 1, function(calls) {
        !anyNA(calls)
    })
}


LabyrinthImpute <- function(file, parents) {
    require(parallel)
    
    prefs <- InitializePreferences()
    prefs$parents <- parents

    LabyrinthImputeHelper(VCF(file), prefs)
}


## TODO(Jason): Remove LabyrinthImputeHelper
## TODO(Jason): Save rds version of file to impute again??
## TODO(Jason): Instead of printing the large table can I do a progress bar one
## step at a time?
LabyrinthImputeHelper <- function(vcf, prefs) {
    start.time <- Sys.time()
    
    if (!inherits(vcf, "vcf")) {
        stop("vcf must be of class 'vcf'")
    }
    if (!inherits(prefs, "prefs")) {
        stop("prefs must be of class 'prefs'")
    }

    ## Determine whether to run in parallel and how many cores to use
    if (prefs$parallel) {
        require(parallel)
        prefs$lapply <- function(...,
                                 mc.preschedule=F,
                                 mc.cores=prefs$cores) {
            mclapply(..., mc.preschedule=mc.preschedule, mc.cores=mc.cores)
        }
    } else {
        prefs$lapply <- function(...,
                                 mc.preschedule=F,
                                 mc.cores=prefs$cores) {
            lapply(...)
        }
    }
    
    ## TODO(Jason): don't impute the parents

    chroms <- vcf$chrom.names
    variants <- vcf$variant.names
    parent.geno <- ResolveHomozygotes(vcf, prefs$parents)

    ## Console output code
    n.chrom <- length(chroms)
    n.variants <- length(variants)
    n.sites <- nrow(parent.geno)
    prefs$n.jobs <- n.variants * n.chrom

    writeLines("\n")
    writeLines("+---------------------------------------------------------------------+")
    writeLines("| LaByRInth: Low-coverage Biallelic R-package Imputation              |")
    writeLines("| Copyright 2017 Jason Vander Woude, Nathan Ryder                     |")
    writeLines("| Licensed under the Apache License, Version 2.0                      |")
    writeLines("| Source code: github.com/Dordt-Statistics-Research/LaByRInth         |")
    writeLines("| Funding recieved from the National Science Foundation (IOS-1238187) |")
    writeLines("+---------------------------------------------------------------------+")
    writeLines("")
    writeLines(paste0(" *  Running in ", ifelse(prefs$parallel, "parallel", "serial")))
    writeLines(paste0(" *  Imputing ",
                      n.variants, " variants at ",
                      n.chrom, " chromosomes (",
                      n.sites, " sites)"))
    writeLines(paste0(" *  ", prefs$n.jobs,
                      " imputations will run (",
                      n.variants, " x ", n.chrom, ")"))

    ## Code from https://stackoverflow.com/questions/27726134/
    ## how-to-track-progress-in-mclapply-in-r-in-parallel-package
    ## TODO(Jason): don't use prefs$fifo, but instead try to use a fifo variable
    ## in the progress.env environment
    progress.env <- new.env()
    prefs$fifo <- ProgressMonitor(progress.env)
    assign("progress", 0.0, envir=progress.env)
    prefs$prog.env <- progress.env

        ## Actually run the imputation
    result <- do.call(cbind,
                      prefs$lapply(
                          variants, function(variant) {
                              LabyrinthImputeSample(vcf, variant, parent.geno, prefs)
                          }, mc.preschedule=FALSE, mc.cores=prefs$cores))

    close(prefs$fifo)

    colnames(result) <- variants
    rownames(result) <- rownames(parent.geno)

    ## Replace spaces and colons with dashes
    vcf.out <- paste0("../LaByRInth_", gsub("[ :]", "-", date()), "_.vcf")
    writeLines(" *  Results imputed")
    writeLines(paste0(" *  Writing output to ", vcf.out))
    
    ## Write imputations to
    ## TODO(Jason): don't use sink()
    sink(vcf.out)
    writeLines(vcf$header[1])  # Add header
    writeLines("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    names <- vcf$header[length(vcf$header)]
    writeLines(names)  # Add header
    names <- str.split(names, "\t")
    
    prefix <- cbind(vcf$variants[, 1:which(names=="INFO")], "GT")

    for (i in 1:n.sites) {
        cat(paste0(prefix[i, ], collapse="\t"))
        cat("\t")
        for (j in 1:n.variants) {
            call <- result[i, j]
            if (is.na(call)) {
                text <- "./."
            } else if (call %in% 1:2) {
                ## TODO(Jason): Discuss what should be done here. Sometimes it
                ## will impute and say that the site is from parent 1 even
                ## though parent 1 is NA. Consider 1A:1159442 and U6202-088. It
                ## was imputed to match LAKIN, but LAKIN is NA in the
                ## parent.geno because LAKIN was called heterozygous
                ## 1/0:3,7. But U6202-088 was called 1/1:0,5 suggesting that it
                ## is pretty likely homozygous alternate. It seems like we
                ## should be utilizing this information maybe like TIGER (see
                ## LB-Impute paper) does. Let these calls show up in the output
                ## file as NA's like is being done now to see how often they
                ## occur. It seems like there is a lot of room for improvement
                ## here. We will have to look at and understand the emission
                ## probabilities better though. Is it valid to have sites that
                ## look like ./1 to represent that we confidently know what one
                ## of the two alleles is but not the other?

                geno <- parent.geno[i, call]
                if (is.na(geno)) {
                    geno <- "."
                }
                text <- paste0(geno, "/", geno)
            } else if (call == 3) {
                text <- "0/1"
            } else {
                stop(paste("Invalid genotype call:", call))
            }
            cat(text)
            if (j != n.variants) {
                cat("\t")
            }
        }
        if (i != n.sites) {
            cat("\n")
        }
    }
    sink()  # turn off sink

    end.time <- Sys.time()
    time <- difftime(end.time, start.time)
    runtime <- as.numeric(time)
    units <- attr(time, "units")
    writeLines(paste(" *  Completed in", round(runtime, 2), units, "\n"))
    invisible(result)  # implicit return
}


LabyrinthImputeSample <- function(vcf, sample, parent.geno, prefs) {
    ## TODO(Jason): Talk to Tintle about why the imputation looks so bad when
    ## imputing the parents. Shouldn't it be mostly correct? It is only about
    ## 65% correct. Nevermind, the fact that there is no chance of a parent
    ## being heterozygous is probably what the difference is since we cleaned
    ## them up
    if (sample %in% prefs$parents) {
        return(rep(match(sample, prefs$parents), nrow(parent.geno)))
    }
    
    chroms <- vcf$chrom.names

    do.call(c, prefs$lapply(chroms, function(chrom) {
        result <- LabyrinthImputeChrom(vcf, sample, chrom, parent.geno, prefs)
        writeBin(1/prefs$n.jobs, prefs$fifo)
        if (!prefs$parallel) {  # if running in serial mode
            prefs$prog.env$progress <- PrintProgress(prefs$fifo, prefs$prog.env$progress)
        }  # else the forked process handles this
        result  # implicit return
    }, mc.preschedule=FALSE, mc.cores=prefs$cores))
}


LabyrinthImputeChrom <- function(vcf, sample, chrom, parent.geno, prefs) {

    if (length(sample) != 1) {
        stop("Length of sample must be 1")
    }
    if (length(chrom) != 1) {
        stop("Length of chrom must be 1")
    }
    
    emission.probs <- GetProbabilities(vcf, sample, chrom, parent.geno, prefs)
    
    site.pos <- sapply(rownames(emission.probs), function(name) {
        as.numeric(str.split(name, ":")[2])
    })

    relevant.sites <- GetRelevantProbabiltiesIndex(vcf, chrom, parent.geno, prefs)
    ## <- GetRelevantProbabiltiesIndex(emission.probs)

    ## If there are not enough markers according to user preference (or if there
    ## are 0), then do not do the imputation and return a path of NA's of the
    ## correct length
    if (sum(relevant.sites) < max(1, prefs$min.markers)) {  # boolean addition
        full.path <- rep(NA_integer_, length(relevant.sites))
    } else {
        names(relevant.sites) <- NULL  # Makes debugging easier
        
        ## distances between relevant sites
        dists <- diff(site.pos[relevant.sites])

        relevant.probs <- emission.probs[relevant.sites, , drop=F]
        class(relevant.probs) <- "probs"

        path <- viterbi(relevant.probs, dists, prefs)

        ## TODO(Jason): add code to compute reverse path and reconcile them. This is
        ## probably worth putting into a function call since there are quite a few
        ## lines of code dedicated to filling in the NA's again.
        
        full.path <- relevant.sites

        path.index <- 1
        ## The missing calls that were not relevant will be filled back in
        ## to create the full path from the 
        for (i in seq_along(relevant.sites)) {
            if (relevant.sites[i]) {  # if the site was relevant
                full.path[i] <- path[path.index]  # set to next call
                path.index <- path.index + 1  # increment call index
            } else {
                full.path[i] <- NA_integer_
            }
        }        
    }

    ## At this stage full.path has entries of 1, 2, 3, or NA which indicates
    ## respectively if a site is homozygous parent 1, homozygous parent 2,
    ## heterozygous, or unknown. Now the entries will be converted to 0

    full.path  # implicit return
}


InitializePreferences <- function() {
    prefs <- list()
    class(prefs)            <- "prefs"

    ## Algorithm parameters
    prefs$parents           <- c("LAKIN", "FULLER")                 
    prefs$resolve.conflicts <- FALSE                                
    prefs$recomb.double     <- FALSE                                
    prefs$read.err          <- 0.05                                 
    prefs$genotype.err      <- 0.05                                 
    prefs$recomb.err        <- 0.05                                 
    prefs$recomb.dist       <- 100000                               
    prefs$min.markers       <- 1
    prefs$states            <- length(prefs$parents) + 1

    ## Logistics
    prefs$quiet             <- FALSE
    prefs$cores             <- 4
    prefs$parallel          <- FALSE

    prefs  # implicit return
}

##' Checking preferences for the correct variable type
##'
##' 
##' @title 
##' @param prefs 
##' @return 
##' @author Jason Vander Woude and Nathan Ryder
ValidatePreferences <- function(prefs) {
    if (!inherits(prefs, "prefs")) {
        stop("prefs must be of class 'prefs'")
    }
    if (length(parents) != 2) {
        stop("exactly 2 parents must be specified")
    }
    if (!is.logical(resolve.conflicts)) {
        stop("resolve.conflicts must be of type logical")
    }
    if (!is.logical(recomb.double)) {
        stop("recomb.double must be of type logical")
    }
    if (!(0 <= read.err < 1) || !(0 <= genotype.err < 1) || !(0 <= recomb.err < 1)) {
        stop("error values must be between 0 and 1")
    }
    if (!is.numeric(recomb.dist) || !(recomb.dist > 0)) {
        stop("recombination distance must be a number greater than 0")
    }
    if (!is.numeric(min.markers) || !(min.markers >= 1)) {
        stop("min.markers must be a number greater than or equal to 1")
    }
    if (
    ## TODO
}


ProgressMonitor <- function(env) {
    local({
        f <- fifo(tempfile(), open="w+b", blocking=T)
        if (!prefs$parallel) {  # don't fork if running serially
            return(f)
        }
        if (inherits(parallel:::mcfork(), "masterProcess")) {
            progress <- 0.0
            while (progress < 1 && !isIncomplete(f)) {
                progress <- PrintProgress(f, progress)
            } 
            parallel:::mcexit()
        }
        f  # implicit return
    }, envir=env)
}


PrintProgress <- function(f, curr.prog) {
    msg <- readBin(f, "double")
    progress <- curr.prog + as.numeric(msg)
    cat(sprintf(" *  Imputation progress: %.2f%%\r", progress * 100))
    progress  # implicit return
}
