 # Map QTLs for muscle weight phenotypes in advanced crosses of
# C57BL/6J (B6) and DBA/2J (D2) inbred strains. In this script, I
# separately analyze data from the F8 cross, data from the F9-F13
# crosses, and data from the combined F8-F13 sample. I compare two
# different methods for mapping QTLs: (1) a simple linear regression
# approach that does not correct for possible confounding due to
# relatedness ('qtl'); and (2) a linear mixed model that uses
# marker-based estimates of pairwise relatedness to correct for
# possible confounding due to unequal relatedness ('QTLRel'). This
# script also allows for the inclusion of a binary interactive
# covariate (with possible values of 0 or 1).
#
# I use this script to run four separate analyses (labeled A through
# D) of the B6 x D2 advanced intercross line:
#
#   (A) Analyze marker data on all chromosomes except the X
#       chromosome, and in all samples.
#
#   (B) Analyze marker data on all chromosomes except the X
#       chromosome, in all samples, in which we model marker effects
#       that differ among males and females (sex interactions).
#
#   (C) Analyze marker data on all chromosomes, including the X
#       chromosome, in males only.
#
#   (D) Analyze marker data on all chromosomes including the X
#       chromosome, in females only.
#
# SCRIPT PARAMETERS
# -----------------
analysis    <- "A"    # Which analysis to run.
num.perm    <- 1000   # Number of replicates for permutation test.
num.markers <- 22000  # Number of markers for QTL mapping.

# Map QTLs separately for these phenotypes.
phenotypes <- c("TA","EDL","gastroc","soleus","tibia")

# Map QTLs in these sets of samples: F8 mice, F9-F13 mice, and all mice.
samples <- list(F8       = function (pheno) pheno$generation == "F8",
                "F9-F13" = function (pheno) pheno$generation != "F8",
                combined = function (pheno) rep(TRUE,nrow(pheno)))

# Set the analysis-specific parameters. Note that the interactive
# covariate should *not* be included as an additive covariate; I will
# include it as an additive covariate below, as required.
cat("Initiating Analysis ",analysis,".\n",sep="")
if (analysis == "A") {

  # (A) QTL MAPPING IN ALL SAMPLES.
  resultsfile   <- "gwscan.combined.RData"
  keep.mice     <- "all"
  drop.X        <- TRUE
  int.covariate <- NULL
  covariates    <- list(F8       = "sex",
                        "F9-F13" = "sex",
                        combined = c("sex","F8"))
} else if (analysis == "B") {

  # (B) QTL MAPPING IN ALL SAMPLES, WITH SEX INTERACTIONS.
  resultsfile   <- "gwscan.combined.sex.RData"
  keep.mice     <- "all"
  drop.X        <- TRUE
  int.covariate <- "sex"
  covariates    <- list(F8       = NULL,
                        "F9-F13" = NULL,
                        combined = "F8")
} else if (analysis == "C") {

  # (C) QTL MAPPING IN MALES ONLY.
  resultsfile   <- "gwscan.combined.males.RData"
  keep.mice     <- "males"
  drop.X        <- FALSE
  int.covariate <- NULL
  covariates    <- list(F8       = NULL,
                        "F9-F13" = NULL,
                        combined = "F8")
} else if (analysis == "D") {

  # (D) QTL MAPPING IN FEMALES ONLY.
  resultsfile   <- "gwscan.combined.females.RData"
  keep.mice     <- "females"
  drop.X        <- FALSE
  int.covariate <- NULL
  covariates    <- list(F8       = NULL,
                        "F9-F13" = NULL,
                        combined = "F8")
}

# Load packages and function definitions.
library(qtl)
library(plyr)
library(abind)
capture.output(library(QTLRel))
source("mapping.functions.R")
source("map.cross.rr.R")

# Initialize the random number generator.
set.seed(7)

# LOAD DATA
# ---------
# Load the phenotype data.
cat("Loading phenotype data.\n")
pheno <- read.phenotypes("../data/pheno.csv",generations)

# Load the marker and genotype data for the F8 cross. I relabel the
# SNPs by their refSNP identifiers (most SNPs correspond to entries in
# the reference database).
cat("Loading genotype data for F8 cross.\n")
geno1        <- read.genotypes("../data/geno.F8.csv")
map1         <- read.map("../data/map.F8.csv",chromosomes)
rows         <- which(!is.na(map1$refSNP))
map1         <- within(map1,snp[rows] <- refSNP[rows])
names(geno1) <- map1$snp

# Load the marker and genotype data for the F9-F13 cohort. When
# possible, I relabel the SNPs by their refSNP identifiers (roughly
# two thirds of the SNPs have corresponding entries in the reference
# database).
cat("Loading genotype data for F9-F13 cohort.\n")
geno2        <- read.genotypes("../data/geno.csv")
map2         <- read.map("../data/map.csv",chromosomes)[1:5]
rows         <- which(!is.na(map2$refSNP))
map2         <- within(map2,snp[rows] <- refSNP[rows])
names(geno2) <- map2$snp

stop()

# Merge the marker data, removing duplicate SNPs, and sorting the SNPs
# by chromosome, then by base pair position along the chromosome.
cat("Merging genotype data.\n")
rows          <- which(!is.element(map1$snp,map2$snp))
map           <- rbind(map1[rows,],map2)
map           <- map[with(map,order(chr,pos)),]
rownames(map) <- NULL

# Merge the genotype data, sorting the SNPs by chromosome, then by
# base pair position along the chromosome.
geno <- rbind.fill(geno1,geno2)
cols <- match(map$snp,colnames(geno))
geno <- geno[cols]

# Remove large data frames that I no longer need (e.g. geno1, geno2).
rm(geno1,geno2,map1,map2)

if (drop.X) {
    
  # Drop the X chromosome from the analysis.
  markers     <- which(map$chr != "X")
  geno        <- geno[,markers]
  map         <- transform(map[markers,],chr = droplevels(chr))
  chromosomes <- levels(map$chr)
} else {

  # Retain the X chromosome, but treat it as an autosome in the QTL
  # mapping.
  chromosomes     <- 1:20
  levels(map$chr) <- chromosomes
}

# I adjust the map positions (i.e. genetic distances) slightly so that
# no two markers have the same position.
map <- jitter.gendist(map,1e-6)

# Select a random subset of the markers, if requested.
if (num.markers < ncol(geno)) {
  markers <- sort(sample(ncol(geno),num.markers))
  geno    <- geno[,markers]
  map     <- map[markers,]
} else
  num.markers <- ncol(geno)

# DISCARD SAMPLES
# ---------------
if (keep.mice == "males") {

  # Retain only samples that correspond to males.
  rows  <- which(pheno$sex == "M")
  pheno <- pheno[rows,]
  geno  <- geno[rows,]
  
} else if (keep.mice == "females") {

  # Retain only samples that correspond to females.
  rows  <- which(pheno$sex == "F")
  pheno <- pheno[rows,]
  geno  <- geno[rows,]
}

# TRANSFORM TRAITS
# ----------------
# Convert sex to a binary covariate so that female = 0 and male = 1.
# Create a binary covariate that is 1 if the individual is from the F8
# cross, and 0 otherwise.
pheno <- transform(pheno,
                   F8  = as.double(generation == "F8"),
                   sex = factor2integer(sex) - 1)

# COMPUTE GENOTYPE PROBABILITIES
# ------------------------------
# Compute the conditional genotype probabilities using QTLRel. To
# accomplish this, I replace the genotypes with allele counts (AA, AB,
# BB become 1, 2, 3, respectively), and I replace any missing values
# with zeros. I compute the genotype probabilities separately for each
# filial generation since they each exhibit different patterns of
# recombination.
cat("Computing conditional genotype probabilities.\n")
G <- genotypes2counts(geno)
G[is.na(G)] <- 0

# Compute the conditional genotype probabilities for the F8 cross.
rows <- which(pheno$generation == "F8")
gp   <- genoProb(G[rows,],map,step = Inf,method = "Haldane",gr = 8)

# Compute the conditional genotype probabilities for the remaining
# generations.
for (i in 9:13) {
  gen  <- paste("F",i,sep = "")
  rows <- which(pheno$generation == gen)
  r    <- genoProb(G[rows,],map,step = Inf,method = "Haldane",gr = i)

  # Combine the genotype probabilities. This involves combining over
  # the arrays of genotype probabilities (the 'pr' list element) along
  # their rows; each row corresponds to a sample.
  gp$pr <- abind(gp$pr,r$pr,along = 1)
}

# Remove some large variables that I no longer need.
rm(G,r)

# QTL MAPPING USING R/QTL
# -----------------------
# Initialize the data structures containing the R/qtl mapping results
# (gwscan.qtl) and the null distribution of the maximum LOD scores
# (perms.qtl).
r          <- vector("list",3)
names(r)   <- names(samples)
gwscan.qtl <- r
perms.qtl  <- r

# Map QTLs using a simple linear regression approach that does not
# correct for possible confounding due to relatedness. I only include
# markers in the QTL mapping that are genotyped in the set of samples.
for (cohort in names(samples)) {

  # Get the table rows corresponding to the set of the samples, and
  # the markers genotyped in the set of samples. Then convert the
  # genotype and phenotype data to the format used by 'qtl'.
  rows    <- which(samples[[cohort]](pheno))
  markers <- which(!all.missing.col(geno[rows,]))
  data    <- rel2qtl(pheno[rows,],geno[rows,markers],map[markers,])
  data    <- rel2qtl.genoprob(data,subset.genoprob(gp,rows,markers))

  # Get the data for the additive covariates. Note that R/qtl requires
  # that the interactive covariate also be included as an additive
  # covariate.
  covars <- c(covariates[[cohort]],int.covariate)
  if (is.null(covars)) {
    data.cov <- NULL
  }
  else data.cov <- data$pheno[covars]

  # Get the data for the interactive covariates. 
  if (is.null(int.covariate)) {
    data.intcov <- NULL
  } else
    data.intcov <- data$pheno[int.covariate]
    
  # Map QTLs using a simple linear regression approach that does not
  # correct for possible confounding to due relatedness.
  cat("Mapping QTLs in",cohort,"cohort using 'qtl'.\n")
  suppressWarnings(gwscan.qtl[[cohort]] <-
      scanone(data,pheno.col = phenotypes,model = "normal",method = "hk",
              addcovar = data.cov,intcovar = data.intcov,use = "all.obs"))
  
  # Perform permutation tests to calculate thresholds for significance,
  # in which we ignore confounding due to relatedness.
  cat("Generating permutations to estimate null distribution in",
      cohort,"cohort.\n")
  suppressWarnings(perms.qtl[[cohort]] <-
      scanone(data,pheno.col = phenotypes,model = "normal",method = "hk",
              addcovar = data.cov,intcovar = data.intcov,use = "all.obs",
              n.perm = num.perm,verbose = FALSE))
}

# Remove some large variables that I no longer need.
rm(data,data.cov,data.intcov)

# QTL MAPPING USING QTLRel
# ------------------------
# Initialize the data structures containing the QTLRel mapping results
# (gwscan.rel) and the parameter estimates for the QTLRel variance 
# components (params).
r          <- vector("list",3)
names(r)   <- names(samples)
gwscan.rel <- r
params     <- r

# Map QTLs using marker-based estimates of pairwise relatedness to
# correct for possible confounding due to unequal relatedness. Again,
# I only consider markers that are genotyped in the set of samples.
for (cohort in names(samples)) {
  cat("Mapping QTLs in",cohort,"cohort using QTLRel:\n")
  cat("trait     n status\n")
 
  # Initialize another data structure containing the QTLRel mapping
  # results for the given sample.
  r                    <- vector("list",length(phenotypes))
  names(r)             <- phenotypes
  gwscan.rel[[cohort]] <- r
  params[[cohort]]     <- r
  
  # Get the rows and markers for the specified set of samples.
  rows    <- which(samples[[cohort]](pheno))
  markers <- which(!all.missing.col(geno[rows,]))

  # Repeat for each phenotype.
  for (phenotype in phenotypes) {

    # Compute LOD scores for all SNPs across the genome using QTLRel.
    out <- map.cross.rr(pheno[rows,],geno[rows,markers],map[markers,],
                        subset.genoprob(gp,rows,markers),phenotype,
                        covariates[[cohort]],int.covariate)
    
    # Get the results of the QTL mapping.
    gwscan.rel[[cohort]][[phenotype]] <- out$gwscan
    params[[cohort]][[phenotype]]     <- out$params
  }

  # Merge the parameter estimates of the variance components for all
  # phenotypes into a single table.
  params[[cohort]] <- do.call(rbind,params[[cohort]])
}

# SAVE RESULTS TO FILE
# --------------------
cat("Saving results to file.\n")
save(file = resultsfile,
     list = c("analysis","phenotypes","samples","covariates","map",
              "int.covariate","num.markers","keep.mice","drop.X",
              "gwscan.qtl","perms.qtl","gwscan.rel","params"))

