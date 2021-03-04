# Phase2.1ardsRPackage
R code to run, validate, and submit the analysis for the ards project.


before installing the package you need to install the 2 following packages 

```
install.packages("caret")
install.packages("icd.data")
library(caret)
library(icd.data)

```

To install this package in R:

```
devtools::install_github("https://github.com/covidclinical/Phase2.1ardsRPackage", subdir="FourCePhase2.1ards", upgrade=FALSE)
```

The main function runAnalysis() has 2 required arguments.

obfuscation => TRUE/ FALSE
obfuscationThreshord => integer 

```
library(FourCePhase2.1ards)

FourCePhase2.1ards::runAnalysis(obfuscation = TRUE, obfuscationThreshord =3)
```


Finally, please submit the results to https://github.com/covidclinical/Phase2.1ardsRSummariesPublic:

```
library(FourCePhase2.1ards)

FourCePhase2.1ards::submitAnalysis()
```

Share with @bertrandmoal your GitHub handle via direct message or the #ards_young Slack channel so you can be added as contributor to the repository.
Note that you would need to use a token to access private repos, see here.
https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token

Briefly, to generate a new token, go to your GitHub settings -> Developer settings -> Personal access tokens -> Generate.

If somehow submitAnalysis() didn’t allow you to upload the results to Phase2.1ardsRSummariesPublic, you can share the results file with  @bertrandmoal via the #ards_young Slack channel.
