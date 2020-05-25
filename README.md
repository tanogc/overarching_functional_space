# Building an overarching functional space

Code to re-create analysis for "A trait space at an overarching-scale yields more conclusive macroecological patterns of functional diversity". This code calculates FD metrics (functional richness, dispersion and evenness) based on FSs covering three decreasing grain-sizes of trait variation: an overarching European scale (541 taxa), a continental context-dependent scale that covers six regions (180 taxa) and six regional-scale FSs (63-108 taxa). It also show the latitudinal response of taxonomic richness and FD metrics. Finally, it provides code to perform null models to assess if FD patterns were independent of taxonomic richness variation and to identify community assembly mechanisms (filtering or overdispersion) along the latitudinal gradient.

## Original article:

Please, use this citation to reference the code:

```
Múrria, C., Iturrarte, G., Gutiérrez-Cánovas, C., 2020. A trait space at an overarching-scale yields 
more conclusive macroecological patterns of functional diversity. Global Ecology and Biogeography
```

## R files description:

* 0_FD_functions.R: R script to estimate Functional Diversity (FD) metrics
* 0_quality_funct_space_fromdist.R: R function for computing the quality of functional dendrogramm and multidimensional functional spaces. This function is a simplified version of the Appendix S1 associated to Maire et al. 2015 (Global Ecol. and Biogeogr.)
* 1_FS script_all.R: Main code to reproduce the results presented in the paper

## Original data
* traits.txt: fuzzy coding traits for 541 river macroinvertebrate taxa
* regional_taxa.txt: regional occurrence of river macroinvertebrate taxa
* site_taxa.txt: site occurrences of river macroinvertebrate taxa

## Dependencies
To run the code and functions from this repository, you need to install the following packages: 'ade4','ape','clue', 'cluster', 'FD', 'geometry',
'gtools','plyr', 'MuMIn', 'vegan'. We can install them by running the following line:

```
install.packages(c("ade4","ape","clue", "cluster", "FD", "geometry",
                   "gtools","plyr", "MuMIn", "vegan"))

```

Please, send questions or problems related with the use of this code to Cayetano Gutiérrez Cánovas (cayeguti@um.es) or Cesc Múrraia (cmurria@gmail.com).

