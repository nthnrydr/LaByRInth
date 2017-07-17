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


nprobs.probs <- function(probs) {
    nrow(probs)
}

nstates.probs <- function(probs) {
    ncol(probs)
}


findPath <- function(variants, probs, prefs) {
    ## probs should be a probs object
    nsites <- nprobs(probs)
    nstates <- nstates(probs)
    trellis <- prefs$markov.order + 2  # viterbi trellis size
    
    ret.val <- list()
    class(ret.val) <- "path"
    ret.val$call <- rep(NA, nsites)
    ret.val$post.prob <- matrix(NA, nrow = nsites, ncol = nstates)  # todo
    
    if (nsites < min.markers) {  # too few markers are present
        return ret.val  # return a path of NA's
    } else {  # there are sufficient markers for imputation
        ## Impute the first genotype
        
        for (i in 1:nsites) {  
            current.trellis <- min(trellis, )  # TODO(Jason): look into this
            current.path <- rep(nstates, current.trellis -
                                         1)  # last path so inc to first

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


# http://homepages.ulb.ac.be/~dgonze/TEACHING/viterbi.pdf
viterbiFunction <- function(probs, dists, nstates, path.size, prefs) {
    old.partial.paths <- matrix(NA, nrow=path.size, ncol=nstates)
    old.partial.probs <- rep(1, nstates)

    new.partial.paths <- old.partial.paths
    new.partial.probs <- old.partial.probs

    ## Initialize the first states based on emission probabilities
    ## TODO(Jason): decide if the initialization should just be equal
    ## probibilities which are determined by a default or a quick check of how
    ## likely a heterozygous call and homozygous call are
    for (j in 1:nstates) {
        old.partial.paths[1, j] <- probs[1, j]
    }
    
    for (i in 2:path.size) {  # for each site in the path
        dist <- dists[i - 1]
        next.paths <- lapply(1:nstates) function(j) {  # for each possible ending state at the site
            ## find the probability of each existing path ending at state j
            extension.probs <- sapply(1:nstates, function(state) {
                ## multiply prob of prev state and transition prob to state j
                partial.probs[state] *
                    transProb(partial.paths[i, state], j, dist, prefs)
            })

            best.prob <- max(extension.probs)
            best.path <- which.max(extension.probs)

            ## return the best path with state j appended
            ## also append the probability (it will be immediately stripped off)
            c(partial.paths[1:(i - 1), best.path], j, best.prob)




            
            ## TODO(Jason): move best.path to higher scope for efficiency
            best.path <- NA
            best.prob <- -1  # every path will be better than this
            for (k in 1:nstates) {  # try extending each existing partial path
                prev.state <- old.partial.paths[i, k]
                current.prob <- old.partial.probs[k] *
                    transProb(j, prev.state, dist, prefs)
                if (current.prob > best.prob) {
                    best.path <- k
                    best.prob <- current.prob * probs[i, j]
                }
            }
            
        }
    }
}


transProb <- function(a, b, dists, prefs) {
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
