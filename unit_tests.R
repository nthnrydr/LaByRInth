## Set up the testing preferences
test.prefs <- list()
class(test.prefs)             <- "prefs"
test.prefs$markov.order       <- 5
test.prefs$parents            <- c("LAKIN", "FULLER")
test.prefs$resolve.conflicts  <- FALSE
test.prefs$read.err           <- 0.05
test.prefs$genotype.err       <- 0.05
test.prefs$recomb.err <- 0.05
test.prefs$recomb.dist <- 100000
test.prefs$recomb.double <- FALSE
test.prefs$window.size <- test.prefs$markov.order + 2
test.prefs$min.samples <- prefs$window.size
test.prefs$min.fraction <- NULL  # Don't remember what this is
test.prefs$min.markers <- NULL # Don't remember what this is




## Check the generatePath algorithm
path.tracker.a <- matrix(c(3,1,1,1,1,1,3,2,2,2,3,3),byrow=T,nrow=3)
path.tracker.a
generatePath(path.tracker.a, 1)
generatePath(path.tracker.a, 2)
generatePath(path.tracker.a, 3)

probs <- matrix(c(0.5,0.4,0.1,0.3,0.3,0.4,0.7,0.1,0.2,0.5,0.4,0.1),nrow=4,byrow=T)
probs
transProb <- function(a, b, dists, prefs) {}


## Check
entry <- "0/0:1,0:1:66:0,3,36"




## Check vcf file creation
alphabet <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l")
reorder(alphabet, 3)
array(reorder(alphabet, 3), dim=c(2,2,3))
mat <- matrix(c("a|b|c", "d/e/f", "g/h|i", "j|k/l"), nrow=2)
mat
submat <- gsub("\\|", "/", mat)
submat
alphabet <- unlist(strsplit(submat, "/"))
alphabet
array(reorder(alphabet, 3), dim=c(2,2,3))




## vcf
vcfobj <- VCF(file="../LakinFuller_GBSv2_20170509.vcf")
ResolveHomozygotes(vcfobj, test.prefs$parents)

apply(genotype[1:15, , ], 1:2, function(alleles) {
    

}
