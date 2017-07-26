## Set up the testing preferences
test.prefs <- list()
class(test.prefs)             <- "prefs"
test.prefs$parents            <- c("LAKIN", "FULLER")
test.prefs$resolve.conflicts  <- FALSE
test.prefs$read.err           <- 0.05
test.prefs$genotype.err       <- 0.05
test.prefs$recomb.err         <- 0.05
test.prefs$recomb.dist        <- 100000
test.prefs$recomb.double      <- FALSE
test.prefs$window.size        <- test.prefs$markov.order + 2
test.prefs$min.fraction       <- NULL  # Don't remember what this is
test.prefs$min.markers        <- 1
test.prefs$states             <- 3
test.prefs$quiet              <- FALSE

## prefs
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
## end prefs



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
resolved.parents <- ResolveHomozygotes(vcfobj, test.prefs$parents)
parent.6b.depths <- Get(vcfobj, "AD", test.prefs$parents, "6B")
lakin.6b.depths <- Get(vcfobj, "AD", "LAKIN", "6B")
parent.map <- GetProbabilities(vcfobj, "LAKIN", resolved.parents, test.prefs)



## imputation logistics
vcf.file <- "../LakinFuller_004.vcf"
vcf.file.big <- "../LakinFuller_001.vcf"
vcf.obj <- VCF(vcf.file)
vcf.obj.big <- VCF(vcf.file.big)

parent.geno <- ResolveHomozygotes(vcf.obj, test.prefs$parents)
sample <- "U6202-080"
chrom <- "1B"

data <- Get(vcf.obj, "GT", sample, chrom)

test.chrom.impute <- LabyrinthImputeChrom(vcf.obj, sample, chrom, parent.geno, test.prefs)
test.sample <- LabyrinthImputeSample(vcf.obj, sample, parent.geno, test.prefs)j

prefs$parallel <- TRUE
prefs$cores <- 4
prefs$quiet <- TRUE
test.impute <- LabyrinthImputeHelper(vcf.obj.big, prefs)
## test.impute <- LabyrinthImpute(vcf.file.big, c("LAKIN", "FULLER"))

## dummy <- LabyrinthImpute(file="../LakinFuller_004.vcf", c("LAKIN","FULLER"))

finalResult <- local({
  f <- fifo(tempfile(), open="w+b", blocking=T)
  if (inherits(parallel:::mcfork(), "masterProcess")) {
    # Child
    progress <- 0.0
    while (progress < 1 && !isIncomplete(f)) {
      msg <- readBin(f, "double")
      progress <- progress + as.numeric(msg)
      cat(sprintf("Progress: %.2f%%\r", progress * 100))
    } 
    parallel:::mcexit()
  }
  numJobs <- 100
  result <- mclapply(1:numJobs, function(...) {
    # Do something fancy here... For this example, just sleep
    Sys.sleep(0.05)
    # Send progress update
    writeBin(1/numJobs, f)
    # Some arbitrary result
    sample(1000, 1)
  })
  close(f)
  result
})

dispmclapply <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, 
                            mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                            mc.cleanup = TRUE, mc.allow.recursive = TRUE) {
    local({
        f <- fifo(tempfile(), open="w+b", blocking=T)
        if (inherits(parallel:::mcfork(), "masterProcess")) {
            progress <- 0.0
            while (progress < 1 && !isIncomplete(f)) {
                msg <- readBin(f, "double")
                progress <- progress + as.numeric(msg)
                cat(sprintf("Progress: %.2f%%\n", progress * 100))
            } 
            parallel:::mcexit()
        }
        result <- mclapply(X, function(...) {
                                        # Do something fancy here... For this example, just sleep
            Sys.sleep(0.05)
                                        # Send progress update
            writeBin(1/numJobs, f)
                                        # Some arbitrary result
            sample(1000, 1)
        })
        close(f)
        result
    })
}
