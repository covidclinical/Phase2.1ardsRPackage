# Phase2.1ardsRPackage
R code to run, validate, and submit the analysis for the ards project.

To install this package in R:

```
devtools::install_github("https://github.com/covidclinical/Phase2.1ardsRPackage", subdir="FourCePhase2.1ards", upgrade=FALSE)
```

The main function runAnalysis() has 2 required arguments.

obfusquation => TRUE/ FALSE
obfuscationThreeshord => integer 

```
library(Phase2.1ardsRPackage)

FourCePhase2.1ards::runAnalysis(obfusquation = TRUE, obfuscationThreeshord =3)
```
