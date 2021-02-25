library(dplyr)

### load pheno_code with ICD10
pheno10_ICD10=read.csv("data-raw/phecode_icd10.csv")
pheno10_ICD10$ICD10=str_replace_all(pheno10_ICD10$ICD10, "\\.", "")

### load pheno_code with ICD9
pheno_ICD9=read.csv("data-raw/phecode_icd9_rolled.csv", sep = ",")
pheno_ICD9 <- pheno_ICD9%>%
  dplyr::select(-c("Rollup","Leaf","Ignore.Bool"))%>%
  rename(ICD.String=ICD9.String,
         ICD=ICD9)

pheno10_ICD10 <- pheno10_ICD10%>%
  rename(ICD.String=ICD10.String,
         ICD=ICD10)

pheno_ICD = rbind(pheno_ICD9,pheno10_ICD10)

### load complication class
comp_class <- read.csv(paste0('data-raw/complication_class.csv'),sep = ";")

## load procedure
sev_proc_icd10 <- read.csv('data-raw/sev_proc_icd10.csv')

## load lab
lab_mapping <- read.csv('data-raw/4CE_loinc.csv')

## load medication
med_code= read.csv("data-raw/4CE_medication.csv")

usethis::use_data(pheno_ICD,comp_class,sev_proc_icd10,lab_mapping, med_code, overwrite = TRUE)
