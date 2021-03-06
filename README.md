# QTL mapping of muscle weight traits in C57BL/6J x DBA/2J AIL mouse study

This repository contains code and data to accompany publication of

Peter Carbonetto, Riyan Cheng, Joseph Gyekis, Clarissa Parker, David Blizard, Abraham Palmer and Arimantas Lionikas (2014). [Discovery and refinement of muscle weight QTLs in B6 x D2 advanced intercross mice.](http://dx.doi.org/10.1152/physiolgenomics.00055.2014) [Physiological Genomics](http://physiolgenomics.physiology.org) 46: 571-582.

### License

Copyright (c) 2014, Peter Carbonetto

The b6d2muscle project repository is free software: you can redistribute
it and/or modify it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html) as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
**without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purpose**. See
[LICENSE](LICENSE) for more details.

### Overview of data files

Here is a brief summary of the files in the [data](data) directory:

+ [pheno.csv](data/pheno.csv) Muscle weight measurements and other
  phenotypes for 891 mice derived from the B6 and D2 inbred strains.

+ [geno.F8.csv](data/geno.F8.csv) Genotype data for 425 F8 mice
  derived from the B6 and D2 inbred strains.

+ [map.F8.csv](data/map.F8.csv) Information about SNPs genotyped in F8
  mice.

+ [geno.csv](data/geno.csv) Genotype data for 435 F9-F13 mice derived
  from the B6 and D2 inbred strains.

+ [map.csv](data/map.csv) Information about SNPs genotyped in F9-F13
  mice.

### Overview of R source code files

Here is a brief summary of the files in the [code](code) directory:

+ [map.qtls.combined.R](code/map.qtls.combined.R) This is the main
  script for mapping muscle weight QTLs in the B6 x D2 advanced
  intercross line. This script separately analyzes data from the F8
  cross, data from the F9-F13 crosses, and data from the combined
  F8-F13 sample, comparing two different methods for mapping QTLs: (1)
  a simple linear regression approach that does not correct for
  possible confounding due to relatedness ('qtl'); and (2) a linear
  mixed model that uses marker-based estimates of pairwise relatedness
  to correct for possible confounding due to unequal relatedness
  ('QTLRel'). This script also includes an option for modeling
  sex-specific genetic effects on muscle weights.

+ [plot.gwscan.combined.R](code/plot.gwscan.combined.R) Plots the QTL
  mapping results generated by the map.qtls.combined.R script.

+ [examine.pheno.R](code/examine.pheno.R) A small script to view the
  distribution of the muscle weight phenotypes in the F8 and F9-F13
  cohorts.

+ [mapping.functions.R](code/mapping.functions.R) Functions for
  QTL mapping and analyzing the experimental cross data.

+ [map.cross.rr.R](code/map.cross.rr.R) Function used to map QTLs
  separately on each chromosome using QTLRel.

### Credits

The R code implementing the analysis procedures was developed by:<br>
[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago<br> 
June 2014
