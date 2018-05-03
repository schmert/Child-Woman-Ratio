# Child-Woman-Ratio

### 03 May 2018 

This repository contains data and code for replicating **Bayesian estimation of total fertility from a population's age-sex structure** by Carl Schmertmann and Mathew Hauer. Preprint at https://osf.io/xqq4g/ , in press at *Statistical Modelling*. 

### R programs include...

* **Stan Kanamari - final version.R** 
(data and model for the Kanamari do Rio Juruá Brazilian indigenous populationç makes Fig 4)

* **Stan GA 2010 counties.R**
(data and hierarchical model for 159 counties in the state of Georgia)

* **Figures GA 2010 counties.R** (makes Figs 5,6,7 and several supplementary plots that are not included in the paper)

* **compare Bayes vs Oasis.R** (makes Fig 8)

* **big test of alternative priors.R** (code for the sensitivity analysis in Appendix B; makes Fig 9)

### Input Data includes..

* **GA Oasis data.csv** (county-level TFR and q_5 estimates from Georgia's Oasis state data system)

* **BirthRates2006-2010.csv** (county-level ASFRs for 2006-2010 from Georgia's Oasis state data system)

* **bltper_5x5/ directory** (holds text files containing HMD life tables for a variety of countries; used in Appendix B)

### Output Data includes...

* **Kanamari fit.RData** (*stanfit* object with results for the small indigenous population in Brazil)

* **Stan GA 2010 counties fit Fri 11Aug17 1410.RData** (*stanfit* object with result for GA counties)

* **big test of alternative priors.csv** (results from sensitivity experiments in Appendix B)

* many graphs in .eps, .pdf, .png format





