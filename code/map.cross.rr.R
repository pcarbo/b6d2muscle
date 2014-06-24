# This is the function that I use in the map.combined.R script to map
# QTLs separately on each chromosome using QTLRel.
map.cross.rr <- function (pheno, geno, map, gp, phenotype, covariates,
                          int.covariate) {
      
  # Only include samples in the analysis for which the phenotype and
  # all covariates are observed. Here, n is the number of samples used
  # for the QTL mapping.
  cols <- c(phenotype,covariates,int.covariate)
  rows <- which(rowSums(is.na(pheno[cols])) == 0)
  n    <- length(rows)

  # Get the interactive covariate data, if an interactive covariate is
  # chosen.
  if (is.null(int.covariate))
    data.intcov <- NULL
  else
    data.intcov <- pheno[rows,int.covariate]
  
  # Compute the (expected) relatedness matrix using *all* the markers,
  # and estimate the variance components using this relatedness
  # matrix. This step is taken so that we can record the fitted
  # parameter estimates for the variance components using the whole
  # genome.
  cat(sprintf("%-7s %3d ",phenotype,n))
  caterase("(relatedness matrix--all markers)")
  R           <- rr.matrix(subset.genoprob(gp,rows))
  ids         <- pheno$id[rows]
  dimnames(R) <- list(ids,ids)

  # Use the relatedness matrix estimated from the marker data to
  # estimate the variance components, and store the parameters
  # corresponding to these variance components.
  caterase("(variance components--all markers)")
  covars <- c(covariates,int.covariate)
  if (is.null(covars))
    params <- estVC(pheno[rows,phenotype],
                    v = list(AA = R,DD = NULL,AD = NULL,HH = NULL,
                             MH = NULL,EE = diag(n)))$par
  else
    params <- estVC(pheno[rows,phenotype],pheno[rows,covars],
                    v = list(AA = R,DD = NULL,AD = NULL,HH = NULL,
                             MH = NULL,EE = diag(n)))$par

  # Initialize a data structure containing the LOD scores for each
  # chromosome.
  chromosomes  <- levels(map$chr)
  scans        <- vector("list",length(chromosomes))
  names(scans) <- chromosomes

  # Map QTLs separately for each chromosome.
  for (i in chromosomes) {
      
    # Get the markers on the chromosome.
    markers <- which(map$chr == i)
    caterase(sprintf("(Mapping %d markers on chromosome %s)  ",
                     length(markers),i))

    # Estimate the matrix of pairwise relatedness using all markers on
    # all (autosomal) chromosomes *except* the current chromosome.
    R <- rr.matrix(subset.genoprob(gp,rows,which(map$chr != i)))
    dimnames(R) <- list(ids,ids)

    # Use the relatedness matrix estimated from the marker data to
    # estimate the variance components.
    if (is.null(covars))
      r <- estVC(pheno[rows,phenotype],
                 v = list(AA = R,DD = NULL,AD = NULL,
                          HH = NULL,MH = NULL,EE = diag(n)))
    else
      r <- estVC(pheno[rows,phenotype],pheno[rows,covars],
                 v = list(AA = R,DD = NULL,AD = NULL,HH = NULL,
                          MH = NULL,EE = diag(n)))
      
    # Once we have the variance components estimated, build a matrix
    # from the variance components. In this case, we only use two
    # variance components: the "AA" component is the n x n relatedness
    # matrix estimated from the marker data (R); and the "EE"
    # component is the n x n identity matrix.
    vc <- r$par["AA"] * r$v[["AA"]] +
          r$par["EE"] * r$v[["EE"]]

    # Map QTLs on the chromosome. Note that QTLRel, unlike R/qtl,
    # automatically includes sex as an additive covariate when we
    # specify it as an interactive covariate, so there is no need to
    # include it along with the additive covariates.
    if (is.null(covariates))
      out <- scanOne(pheno[rows,phenotype],gdat = geno[rows,markers],
                     prdat = subset.genoprob(gp,rows,markers),vc = vc,
                     intcovar = data.intcov,test = "None")
    else
      out <- scanOne(pheno[rows,phenotype],pheno[rows,covariates],
                     geno[rows,markers],subset.genoprob(gp,rows,markers),
                     vc,data.intcov,test = "None")
    
    # Get the LOD scores and the proportion of variance explained by
    # each SNP.
    chr.scan     <- empty.scanone(map[markers,])
    chr.scan$lod <- out$p/(2*log(10))
    chr.scan$pve <- out$v/100
  
    # Get the additive and dominance effects of all the SNPs. If we've
    # included an interactive covariate, then there are four effects
    # (regression coefficients) for each marker.
    if (is.null(int.covariate)) {
      chr.scan$add <- sapply(out$parameters,"[[","a")
      chr.scan$dom <- sapply(out$parameters,"[[","d")
    }
    else {

      # Get the additive and dominance effects when interactive
      # covariate is 0.
      chr.scan$add0 <- sapply(out$parameters,"[[","a")
      chr.scan$dom0 <- sapply(out$parameters,"[[","d")

      # Get the additive and dominance effects when interactive
      # covariate is 1.
      chr.scan$add1 <- sapply(out$parameters,"[[","intcovar.a")
      chr.scan$dom1 <- sapply(out$parameters,"[[","intcovar.d")
    }
    
    # Save the QTL mapping results for the chromosome.
    scans[[i]] <- chr.scan
  }
  
  # Merge the mapping results from all chromosomes into a single table.
  gwscan           <- do.call(rbind,scans)
  rownames(gwscan) <- do.call(c,lapply(scans,rownames))
  cat("\n")
  
  # Output the mapping results and the estimated parameters
  # corresponding to the variance components.
  return(list(gwscan = gwscan,params = params))
}
