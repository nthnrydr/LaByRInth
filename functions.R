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
###   $drp


nsites.probs <- function(probs) {
    nrow(probs)
}


nstates.probs <- function(probs) {
    ncol(probs)
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
    } else if (a %in% 1:2 && b %in% 1:2 && !prefs$drp) {
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
increment <- function(x, max, min = 1) {
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
