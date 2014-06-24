# A small script to view the distribution of the muscle weight
# phenotypes in the F8 and F9-F13 cohorts.
library(lattice)
library(Hmisc)
source("mapping.functions.R")

# SCRIPT PARAMETERS.
panels <- list(TA      = list(xlim = c(25,65)),
               EDL     = list(xlim = c(4,13)),
               gastroc = list(xlim = c(65,165)),
               soleus  = list(xlim = c(3,11)))

# Load the phenotype data.
pheno <- read.phenotypes("../data/pheno.csv",generations)

# Create a new column in the phenotype data that gives the "cohort",
# which is either F8 or F9-F13.
pheno <- cbind(pheno,data.frame(cohort = factor(pheno$generation != "F8")))
levels(pheno$cohort) <- c("F8","F9-F13")

# Set up the graphics device.
trellis.device(height = 9.5,width = 6)
trellis.par.set(list(fontsize = list(text = 9),
                     layout.widths = list(right.padding = -1),
                     layout.heights = list(top.padding = -3,
                                           bottom.padding = 0)))

n <- length(panels)
for (i in 1:n) {
  phenotype <- names(panels)[i]
  
  # Show the phenotype distribution in F8 females using a histogram.
  rows <- with(pheno,which(generation == "F8" & sex == "F"))
  print(histogram(formula(paste("~",phenotype)),pheno[rows,],
                  col = "dodgerblue",border = "dodgerblue",nint = 10,
                  scales = list(x = list(tck = 0.5),y = list(tck = 0.5)),
                  xlab = phenotype,ylab = "F8 females",
                  type = "count",xlim = panels[[phenotype]]$xlim,
                  ylim = c(0,70)),
        split = c(i,1,n,8),
        more = TRUE)

  # Show the phenotype distribution in F9-F13 females using a histogram.
  rows <- with(pheno,which(generation != "F8" & sex == "F"))
  print(histogram(formula(paste("~",phenotype)),pheno[rows,],
                  col = "darkorange",border = "darkorange",nint = 10,
                  scales = list(x = list(tck = 0.5),y = list(tck = 0.5)),
                  xlab = phenotype,ylab = "F9-F13 females",
                  type = "count",xlim = panels[[phenotype]]$xlim,
                  ylim = c(0,70)),
        split = c(i,2,n,8),
        more = TRUE)

  # Show the phenotype distribution in F8 males using a histogram.
  rows <- with(pheno,which(generation == "F8" & sex == "M"))
  print(histogram(formula(paste("~",phenotype)),pheno[rows,],
                  col = "dodgerblue",border = "dodgerblue",nint = 10,
                  scales = list(x = list(tck = 0.5),y = list(tck = 0.5)),
                  xlab = phenotype,ylab = "F8 males",
                  type = "count",xlim = panels[[phenotype]]$xlim,
                  ylim = c(0,70)),
        split = c(i,3,n,8),
        more = TRUE)
               
  # Show the phenotype distribution in F9-F13 males using a histogram.
  rows <- with(pheno,which(generation != "F8" & sex == "M"))
  print(histogram(formula(paste("~",phenotype)),pheno[rows,],
                  col = "darkorange",border = "darkorange",nint = 10,
                  scales = list(x = list(tck = 0.5),y = list(tck = 0.5)),
                  xlab = phenotype,ylab = "F9-F13 males",
                  type = "count",xlim = panels[[phenotype]]$xlim,
                  ylim = c(0,70)),
        split = c(i,4,n,8),
        more = TRUE)

  # Show the distribution of the phenotype in the F8 and F9-F13
  # cohorts using a box-percentile plot.
  print(bwplot(formula(paste("~",phenotype,"| cohort + sex")),pheno,
               probs = c(0.01,0.05,0.125,0.25,0.375),
               panel = panel.bpplot,means = FALSE,nout = 0.003,
               scat1d.opts = list(lwd = 4,tfrac = 0.01,col = "black"),
               xlab = phenotype,ylab = "",layout = c(1,4),
               xlim = panels[[phenotype]]$xlim,
               scales = list(x = list(tck = 0.5))),
        split = c(i,2,n,2),
        more = TRUE)

}
