# QTL mapping of muscle weight traits in C57BL/6J x DBA/2J AIL mouse study

Short description of this repository goes here.

###License

Copyright (c) 2014, Peter Carbonetto

The lgsmfear project repository is free software: you can redistribute
it and/or modify it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html) as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
**without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purpose**. See
[LICENSE](LICENSE) for more details.

###Overview of data files

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

###Overview of R source code files

Here is a brief summary of the files in the [code](code) directory:

+ [examine.pheno.R](code/examine.pheno.R) A small script to view the
  distribution of the muscle weight phenotypes in the F8 and F9-F13
  cohorts.

+ [mapping.functions.R](code/mapping.functions.R) Functions for
  QTL mapping and analyzing the experimental cross data.

###Credits

The R code implementing the analysis procedures was developed by:<br>
[Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago<br> 
June 2014
