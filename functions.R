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


findPath <- function(variants, probs, prefs) {
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


transProb <- function(recomb, dists, prefs) {
    if (recomb) {
        0.5 * (1 - exp(-dist/prefs$recomb.dist))
    } else {
        0.5 * (1 + exp(-dist/prefs$recomb.dist))
    }
}


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
