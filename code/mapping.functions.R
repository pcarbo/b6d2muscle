# This file contains a bunch of functions I wrote that are useful for
# the QTL mapping using R/qtl and QTLRel.
chromosomes <- c(as.character(1:19),"X")
generations <- c("F8","F9","F10","F11","F12","F13")

# ----------------------------------------------------------------------
# Output the string using 'cat', then move the cursor back to the
# beginning of the string so that subsequent output will overwrite
# this string.
caterase <- function (s) {
  cat(s,rep("\b",nchar(s)),sep="")
}

# ----------------------------------------------------------------------
# Convert a factor to a vector of integers.
factor2integer <- function (x)
  match(x,levels(x))

# ----------------------------------------------------------------------
# Returns a logical vector such that an entry is TRUE if all the
# entries in the corresponding list element (or column of the table)
# are missing---that is, all the entries are set to NA.
all.missing.col <- function (d)
    sapply(d,function(x) all(is.na(x)))

# ----------------------------------------------------------------------
# Returns a data frame containing the phenotype data stored in a CSV
# file. I convert the generation to a factor manually so that I can
# control the order of generations F8 through F13 in the factor.
read.phenotypes <- function (file, generations)
  within(read.csv(file,comment.char="#",as.is = c("id","generation"),
                  check.names = FALSE),
         generation <- factor(generation,levels = generations))

# ----------------------------------------------------------------------
# Returns a data frame containing the genotype data stored in a CSV
# file. I remove the first two columns of the genotype data frame
# containing the ID and generation, and set the row names for the
# genotype data to the mouse IDs.
read.genotypes <- function (file) {
  geno <- read.csv(file,comment.char="#",as.is = "id",check.names = FALSE)
  row.names(geno) <- geno$id
  geno            <- geno[-(1:2)]
  return(geno)
}

# ----------------------------------------------------------------------
# Returns a data frame containing the marker data stored in a CSV file.
# Here I convert the chromosomes to factors manually (as opposed to
# letting the read.csv function handle this) to make sure that the
# chromosomes are ordered properly.
read.map <- function (file, chromosomes)
  within(read.csv(file,comment.char="#",check.names = FALSE,
                  stringsAsFactors = FALSE),
         chr <- factor(chr,chromosomes))

# ----------------------------------------------------------------------
# Convert the genotypes to a matrix of allele counts.
genotypes2counts <- function (d)
  as.matrix(data.frame(lapply(d,factor2integer),
                       row.names = rownames(d),
                       check.names = FALSE))

# ----------------------------------------------------------------------
# Adjust the map positions (i.e. genetic distances) slightly (by the
# amount j) so that no two markers have the same position.
jitter.gendist <- function (map, j) {

  # Get the set of chromosomes.
  chromosomes <- levels(map$chr)

  # Repeat for each chromosome.
  for (i in chromosomes) {

    # Get the markers on the chromosome.
    markers  <- which(map$chr == i)
    n        <- length(markers)

    # Adjust the genetic distances slightly according to scalar j.
    map[markers,"dist"] <- map[markers,"dist"] + cumsum(rep(j,times = n))
  }

  return(map)
}

# ----------------------------------------------------------------------
# Get the genotype probabilities for the specified samples and markers.
subset.genoprob <- function (gp, samples = NULL, markers = NULL) {

  # If the set of markers is the null object, use all the markers.
  if (is.null(markers)) {
    n       <- length(gp$chr)
    markers <- 1:n
  }

  # If the set of samples is the null object, use all the samples.
  if (is.null(samples)) {
    n       <- nrow(gp$pr)
    samples <- 1:n
  }
  
  class(gp) <- c("Pr","list")
  return(within(gp,{
    pr   <- pr[samples,,markers]
    chr  <- chr[markers]
    dist <- dist[markers]
    snp  <- snp[markers]
  }))
}

# ----------------------------------------------------------------------
# Convert the QTL experiment data from the format used in the QTLRel
# library to the format used in qtl library. The return value is a
# 'cross' object that keeps track of all the data in a single QTL
# experiment; for more details, see the help for function read.cross
# in the qtl library.
#
# Here I assume that the alleles are A and B and the genotypes are
# labeled as AA, AB and BB. qtl requires a paternal grandmother
# ("pgm") phenotype, so a column in the table for this phenotype is
# included if it is missing, and the entries of this column are set to
# zero.
#
# IMPORTANT NOTE: this function does not currently work for markers on
# the X chromosome.
rel2qtl <- function (pheno, geno, map) {

  # Get the set of chromosomes.
  chromosomes <- levels(map$chr)
  
  # Convert the genotypes to a matrix of integers.
  geno <- genotypes2counts(geno)
  
  # Add the paternal grandmother ("pgm") phenotype if it does not
  # already exist.
  if (!is.element("pgm",names(pheno)))
    pheno <- cbind(pheno,pgm = 0)
  
  # Initialize the cross object, setting the genotypes to an empty list.
  cross <- list(geno = list(),pheno = pheno)
  class(cross) <- c("f2","cross")

  # Set the alleles.
  attributes(cross)$alleles <- c("A","B")
  
  # Split the markers by chromosome.
  for (i in chromosomes) {

    # Get the markers on the chromosome.
    markers <- which(map$chr == i)
    
    # Get the genotype and map data corresponding to these markers.
    # Note that I need to coerce the genotype data to be a matrix in
    # the exceptional case when there is only a single marker
    # genotyped on the chromosome.
    m                <- map[markers,]
    d                <- list(data = as.matrix(geno[,markers]),map = m$dist)
    colnames(d$data) <- m$snp
    names(d$map)     <- m$snp
    class(d)         <- "A"

    # Store the genotypes for markers on the chromosome.
    cross$geno[[i]] <- d
  }
  
  return(cross)
}

# ----------------------------------------------------------------------
# Convert the genotype probabilities from the format used by QTLRel to
# the format used by qtl, and store the genotype probabilities in the
# cross object. Note that I currently do not properly handle the X
# chromosome.
#
# IMPORTANT NOTE: this function does not currently work for markers on
# the X chromosome.
rel2qtl.genoprob <- function (cross, gp) {

  # Transpose the array.
  prob <- aperm(gp$pr,c(1,3,2))

  # Get the set of chromosomes.
  chromosomes <- unique(gp$chr)

  # Split the genotype probabilities by chromosome.
  for (i in chromosomes) {

    # Get the data for all markers on this chromosome.
    geno <- cross$geno[[i]]
    
    # Get the markers on the chromosome.
    markers <- which(gp$chr == i)

    # Get the genotype probabilities for these markers.
    geno$prob                 <- prob[,markers,]
    attributes(geno$prob)$map <- geno$map
    dimnames(geno$prob)       <- list(NULL,
                                      names(geno$map),
                                      c("AA","AB","BB"))

    # Store the new genotype data.
    cross$geno[[i]] <- geno
  }
  
  return(cross)
}

# ----------------------------------------------------------------------
# Returns the matrix product A*A'.
matrix.square <- function (A)
  return(A %*% t(A))

# ----------------------------------------------------------------------
# Use all the markers to estimate the n x n pairwise relatedness
# matrix, which I define to be 2 times the matrix of kinship
# coefficients (for the definition of kinship coefficients, see Lange,
# 2002). To allow for uncertainty in the genotypes at the markers, I
# compute the expected value of this matrix using the genotype
# probabilities provided as output from the 'genoProb' function in
# QTLRel.
rr.matrix <- function (gp) { 

  # Get the number of samples (n) and the number of markers (p).
  d <- dim(gp$pr)
  n <- d[1]
  p <- d[3]
  
  # Get the genotype probabilities.
  pAA <- gp$pr[,1,]
  pAB <- gp$pr[,2,]
  pBB <- gp$pr[,3,]

  # Get the probability of the genotype AB averaged over all the markers.
  mAB <- matrix(rep(rowMeans(pAB),times = n),n,n)

  # Return the expected value of the pairwise relatedness matrix.
  return(2*(matrix.square(pAA)/p + matrix.square(pBB)/p) 
         + mAB + t(mAB) - matrix.square(pAB)/p)
}

# ----------------------------------------------------------------------
# This function creates a "scanone" object that will be used to store
# QTL mapping results based on the SNP data provided as input. (For
# further info, see the 'scanone' function in the 'qtl' package.) 
empty.scanone <- function (map) {
  d            <- map[c("chr","dist")]
  row.names(d) <- map[["snp"]]
  names(d)     <- c("chr","pos")
  class(d)     <- c("scanone","data.frame")
  return(d) 
}

# ----------------------------------------------------------------------
# Analyze the QTL experiment data using scanOne from the QTLRel
# library. The function returns a data frame containing the QTL
# mapping results, specifically the LOD scores and the estimated
# additive and dominance effects corresponding to individual markers.
map.cross.rel <- function (pheno, geno, map, gp, vc, phenotype, covariates) {
  out    <- scanOne(pheno[,phenotype],pheno[,covariates],
                    geno,gp,vc,test = "None")
  gwscan <- empty.scanone(map)

  # Get the LOD scores and the proportion of variance explained by each SNP. 
  gwscan$lod <- out$p/(2*log(10))
  gwscan$pve <- out$v/100
  
  # Get the additive and dominance effects of all the SNPs.
  gwscan$add <- sapply(out$parameters,"[[","a")
  gwscan$dom <- sapply(out$parameters,"[[","d")
  
  return(gwscan)
}
