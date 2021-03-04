
#' Runs the analytic workflow for the ards project
#'
#' @keywords 4CE
#'
#' @param obfusquation
#' @param obfuscationThreeshord
#'
#' @return
#' @export
#' @import dplyr tidyr stringr icd caret DT tidyverse icd.data
#' @examples
runAnalysis <- function(obfuscation = TRUE, obfuscationThreshord =3) {

    ## make sure this instance has the latest version of the quality control and data wrangling code available
    devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)

    ## get the site identifier assocaited with the files stored in the /4ceData/Input directory that
    ## is mounted to the container
    currSiteId = FourCePhase2.1Data::getSiteId()

    ## run the quality control
    FourCePhase2.1Data::runQC(currSiteId)


    ## DO NOT CHANGE ANYTHING ABOVE THIS LINE

    ## ========================================
    ## PART 1: Read in Data Tables
    ## ========================================
    message("\nPlease ensure that your working directory is set to /4ceData")
    message("Reading in Input/LocalPatientSummary.csv and Input/LocalPatientObservations.csv...")

    LocalPatientClinicalCourse=read.csv(paste0(FourCePhase2.1Data::getInputDataDirectoryName(),"/LocalPatientClinicalCourse.csv"), na.strings = c('1900-01-01', '1/1/1900'))
    LocalPatientObservations=read.csv(paste0(FourCePhase2.1Data::getInputDataDirectoryName(),"/LocalPatientObservations.csv"),na.strings = '-999')
    LocalPatientSummary=read.csv(paste0(FourCePhase2.1Data::getInputDataDirectoryName(),"/LocalPatientSummary.csv"), na.strings = c('1900-01-01', '1/1/1900'))


    ## ========================================
    ## PART 2: load doc and reformat
    ## ========================================

    ## obfusquation
    obfuscationValue <- -99

    ### reformat ICD10  : keep only the 3 first digit
    LocalPatientObservations <- LocalPatientObservations %>%
        dplyr::mutate(concept_code = as.character(concept_code))%>%
        dplyr::mutate( concept_code = ifelse(concept_type %in% c("DIAG-ICD10","DIAG-ICD9"),substr(concept_code, 0, 3),concept_code))

    ### load ICD10
    ICD10=icd.data::icd10cm2016 %>%
        dplyr::mutate(code=as.factor(code))

    ### load ICD9
    ICD9=icd.data::icd9cm_hierarchy %>%
        dplyr::mutate(code=as.factor(code))

    # path_pack="/RDevelopment/Phase2.1ardsRPackage/FourCePhase2.1ards"
    #
    # ### load pheno_code with ICD10
    # pheno10_ICD10=read.csv(paste0(path_pack,"/doc/phecode_icd10.csv"))
    # pheno10_ICD10$ICD10=str_replace_all(pheno10_ICD10$ICD10, "\\.", "")
    #
    # ### load pheno_code with ICD9
    # pheno_ICD9=read.csv(paste0(path_pack,"/doc/phecode_icd9_rolled.csv"), sep = ",")
    # pheno_ICD9 <- pheno_ICD9%>%
    #     dplyr::select(-c("Rollup","Leaf","Ignore.Bool"))%>%
    #     rename(ICD.String=ICD9.String,
    #            ICD=ICD9)
    #
    # pheno10_ICD10 <- pheno10_ICD10%>%
    #     rename(ICD.String=ICD10.String,
    #            ICD=ICD10)
    #
    # pheno_ICD = rbind(pheno_ICD9,pheno10_ICD10)
    #
    # ### load complication class
    # comp_class <- read.csv(paste0(path_pack,'/doc/complication_class.csv'),sep = ";")
    #
    # ## load procedure
    # sev_proc_icd10 <- read.csv(paste0(path_pack,'/doc/sev_proc_icd10.csv'))
    #
    # ## load lab
    # lab_mapping <- read.csv(paste0(path_pack,'/doc/4CE_loinc.csv'))
    #
    # ## load medication
    # med_code= read.csv(paste0(path_pack,"/doc/4CE_medication.csv"))

    data(pheno_ICD)
    data(comp_class)
    data(sev_proc_icd10)
    data(lab_mapping)
    data(med_code)

    ## input ###
    med_severe = c("SIANES","SICARDIAC")
    loinc_pao2 = "2703-7"
    ARDS=c('J80','518.82')


    # P1 : january - july
    # P2 : august - december

    start_date_p1=as.POSIXct(as.Date("2020-01-01"), tz = Sys.timezone())
    end_date_p1=as.POSIXct(as.Date("2020-07-01"), tz = Sys.timezone())
    start_date_p2=as.POSIXct(as.Date("2020-12-31"), tz = Sys.timezone())

    # date of inclusion
    last_date_inclusion=as.POSIXct(as.Date("2020-12-01"), tz = Sys.timezone())


    ### limit before / after
    limit_d_befor= -8
    limit_d_after= 90


    ## ========================================
    ## PART 3: Group selection
    ## ========================================


    #remove patient entry after the last date of inclusion
    LocalPatientSummary <- LocalPatientSummary %>%
        dplyr::mutate(admission_date = as.POSIXct(admission_date)) %>%
        dplyr::filter( admission_date <= last_date_inclusion)

    LocalPatientClinicalCourse <- LocalPatientClinicalCourse %>%
        dplyr::filter( patient_num %in% unique(LocalPatientSummary$patient_num))

    LocalPatientObservations <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% unique(LocalPatientSummary$patient_num))


    ## select patient by age
    LocalPatientSummary <- LocalPatientSummary %>%
        dplyr::filter( age_group %in% c("18to25","26to49","50to69","70to79","80plus")) %>%
        dplyr::mutate( age_group_spe = dplyr::case_when(
            age_group %in% c("18to25","26to49")~ "18_49",
            age_group %in% c("50to69","70to79","80plus")~ "sup_49"),
            periode_group = dplyr::case_when(
                admission_date >= start_date_p1 &  admission_date < end_date_p1 ~ "P1",
                admission_date >= end_date_p1 &  admission_date <= start_date_p2 ~ "P2"))


    ### patients with ARDS
    pat_ARDS <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::filter(grepl(paste(ARDS,collapse="|"), concept_code) & days_since_admission >=0) %>%
        dplyr::select( patient_num)%>%
        unique()


    ### patients with ventilation
    pat_vent <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::filter(concept_type == "PROC-ICD10"
                      & days_since_admission >=0
                      & concept_code %in% sev_proc_icd10$code) %>%
        dplyr::select( patient_num)%>%
        unique()

    ### patients with severe medication
    pat_med_severe <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::filter(concept_code %in% c("SIANES","SICARDIAC") & days_since_admission >=0) %>%
        dplyr::select( patient_num)%>%
        unique()

    ### add groups to file

    LocalPatientSummary <- LocalPatientSummary %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::mutate(ARDS = if_else(patient_num %in% pat_ARDS$patient_num ,1,0),
                      MED_SEVERE = if_else(patient_num %in% pat_med_severe$patient_num ,1,0),
                      PROC_SEVERE = if_else(patient_num %in% pat_vent$patient_num ,1,0),
                      GROUP = dplyr::case_when(
                          ARDS==1 & age_group_spe== "18_49" ~ "ARDS_18_49",
                          ARDS==1 & age_group_spe== "sup_49" ~ "ARDS_sup_49",
                          ARDS==0 & MED_SEVERE==0 & PROC_SEVERE==0 & age_group_spe== "18_49" ~ "NO_ARDS_18_49",
                          ARDS==0 & MED_SEVERE==0 & PROC_SEVERE==0 & age_group_spe== "sup_49" ~ "NO_ARDS_sup_49",
                          ARDS==0 & (MED_SEVERE==1 | PROC_SEVERE==1) & age_group_spe== "18_49" ~ "OTHERS_18_49",
                          ARDS==0 & (MED_SEVERE==1 | PROC_SEVERE==1) & age_group_spe== "sup_49" ~ "OTHERS_sup_49"))


    LocalPatientObservations <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::left_join(LocalPatientSummary[,c("patient_num","periode_group","GROUP","ARDS","PROC_SEVERE","MED_SEVERE")])

    LocalPatientClinicalCourse <- LocalPatientClinicalCourse %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::left_join(LocalPatientSummary[,c("patient_num","periode_group","GROUP","ARDS","PROC_SEVERE","MED_SEVERE")])%>%
        dplyr::mutate(calendar_date = as.POSIXct(calendar_date))


    ## ========================================
    ## PART 4 : Analysis
    ## ========================================


    ## patient by groups
    num_pat <- LocalPatientSummary %>%
        dplyr::group_by( GROUP) %>%
        dplyr::summarise(npat=n_distinct(patient_num))

    ###output
    output<- num_pat %>%
        dplyr::rename(value=npat)%>%
        dplyr::mutate(concept="patient",
                      variable="number",
                      periode_group="all")%>%
        data.frame()

    output_p <- LocalPatientSummary %>%
        dplyr::group_by( GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::mutate(concept="patient",variable="number")%>%
        data.frame()

    output_gen=rbind(output,output_p)


    ## patient by groups
    output_day <- LocalPatientClinicalCourse %>%
        dplyr::filter( in_hospital == 1) %>%
        dplyr::group_by(GROUP,days_since_admission) %>%
        dplyr::summarise(value=n_distinct(patient_num)) %>%
        dplyr::mutate(concept="patient",
                      variable="number",
                      periode_group="all")%>%
        data.frame()

    output_day_p <- LocalPatientClinicalCourse %>%
        dplyr::filter( in_hospital == 1) %>%
        dplyr::group_by( GROUP,periode_group,days_since_admission) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::mutate(concept="patient",
                      variable="number")%>%
        data.frame()

    output_day=rbind(output_day,output_day_p)


    ## quartil of number of associated ICD10 code/Procedure/medication/labs

    el_freq<- LocalPatientObservations %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(patient_num,GROUP,concept_type)%>%
        dplyr::summarise(freqi= n())%>%
        unique()%>%
        dplyr::group_by(GROUP,concept_type)%>%
        dplyr::summarise(mean_freq=mean(freqi,na.rm = TRUE),
                         std_freq = sd(freqi,na.rm = TRUE))%>%
        pivot_longer(c(mean_freq,std_freq),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::rename(concept=concept_type)%>%
        data.frame()


    el_freq_p<- LocalPatientObservations %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(patient_num,GROUP,concept_type,periode_group)%>%
        dplyr::summarise(freqi= n())%>%
        unique()%>%
        dplyr::group_by(GROUP,concept_type,periode_group)%>%
        dplyr::summarise(mean_freq=mean(freqi,na.rm = TRUE),
                         std_freq = sd(freqi,na.rm = TRUE))%>%
        pivot_longer(c(mean_freq,std_freq),names_to = "variable", values_to = "value")%>%
        dplyr::rename(concept=concept_type)%>%
        data.frame()

    output_gen=rbind(output_gen,el_freq,el_freq_p)


    ## quartil of number of associated ICD10 code/Procedure/medication/labs

    # duration
    duration<-  LocalPatientClinicalCourse %>%
        dplyr::group_by(patient_num,periode_group,GROUP) %>%
        dplyr::arrange(patient_num,days_since_admission)%>%
        dplyr::mutate(temp=in_hospital-lag(in_hospital, default = first(in_hospital)))%>%
        dplyr::summarise(d_hosp= sum(in_hospital),
                         length_duration_1=ifelse(min(days_since_admission[temp==-1]) == Inf,max(days_since_admission),min(days_since_admission[temp==-1])),
                         still_at_hospital= min(days_since_admission[temp==-1]) == Inf,
                         number_rehospit= sum(temp==1),
                         delay_1rehospit_1admin= min(days_since_admission[temp==1]),
                         delay_hospi_rehosp= delay_1rehospit_1admin-length_duration_1)


    ## Length of stay per patient

    los<-  duration %>%
        dplyr::group_by(GROUP)%>%
        dplyr::summarise(mean_los=mean(length_duration_1,na.rm = TRUE),
                         std_los = sd(length_duration_1,na.rm = TRUE))%>%
        pivot_longer(c(mean_los,std_los),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()


    los_p<-  duration %>%
        dplyr::group_by(GROUP, periode_group)%>%
        dplyr::mutate(length_duration_1= ifelse(is.infinite(length_duration_1),NA,length_duration_1))%>%
        dplyr::summarise(mean_los=mean(length_duration_1,na.rm = TRUE),
                         std_los = sd(length_duration_1,na.rm = TRUE))%>%
        pivot_longer(c(mean_los,std_los),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,los,los_p)


    ##  Severity

    sev<-  LocalPatientSummary %>%
        dplyr::mutate(severe = if_else(severe == 1, "ever_severe","never_severe") ) %>%
        dplyr::group_by(GROUP,severe) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=severe) %>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()


    sev_p<-  LocalPatientSummary %>%
        dplyr::mutate(severe = if_else(severe == 1, "ever_severe","never_severe") ) %>%
        dplyr::group_by(GROUP,severe,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=severe) %>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,sev,sev_p)


    ### Diagnostic: Percentage of patient without any diagnostic


    diag_no_day <- LocalPatientObservations %>%
        dplyr::filter(concept_type %in%  c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::group_by(GROUP) %>%
        dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
        dplyr::right_join(output[,c("value","GROUP")],by =c("GROUP"))%>%
        dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
        dplyr::mutate(nodiag_npat=value-count_pat) %>%
        dplyr::select(- c(value,count_pat)) %>%
        dplyr::rename(value=nodiag_npat) %>%
        dplyr::mutate(concept="DIAG-ICD",
                      variable="npat_nodiag",
                      periode_group="all")%>%
        data.frame()

    diag_no_day_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::group_by(GROUP,periode_group) %>%
        dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
        dplyr::right_join(output[,c("value","GROUP")],by =c("GROUP"))%>%
        dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
        dplyr::mutate(nodiag_npat=value-count_pat) %>%
        dplyr::select(- c(value,count_pat)) %>%
        dplyr::rename(value=nodiag_npat) %>%
        dplyr::mutate(concept="DIAG-ICD",
                      variable="npat_nodiag")%>%
        data.frame()

    output_gen=rbind(output_gen,diag_no_day,diag_no_day_p)

    ## Procedure:

    ### Average number of procedure per day

    proc_day <- LocalPatientObservations %>%
        dplyr::filter(concept_type %in% c("PROC-ICD10","PROC-ICD9")) %>%
        dplyr::group_by(patient_num,GROUP,days_since_admission) %>%
        dplyr::summarise(count_proc=n()) %>%
        dplyr::right_join(LocalPatientClinicalCourse[,c("patient_num","GROUP","days_since_admission")],by =c("patient_num","GROUP","days_since_admission"))%>%
        dplyr::mutate(count_proc=ifelse(is.na(count_proc),0,count_proc))%>%
        dplyr::group_by(GROUP,days_since_admission) %>%
        dplyr::summarise(value=mean(count_proc,na.rm = TRUE))%>%
        dplyr::mutate(concept="PROC-ICD",
                      variable="average",
                      periode_group="all")%>%
        data.frame()

    proc_day_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9")) %>%
        dplyr::group_by(patient_num,GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(count_proc=n()) %>%
        dplyr::right_join(LocalPatientClinicalCourse[,c("patient_num","GROUP","days_since_admission")],by =c("patient_num","GROUP","days_since_admission"))%>%
        dplyr::mutate(count_proc=ifelse(is.na(count_proc),0,count_proc))%>%
        dplyr::group_by(GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(value=mean(count_proc,na.rm = TRUE))%>%
        dplyr::mutate(concept="PROC-ICD",
                      variable="average")%>%
        data.frame()

    output_day=rbind(output_day,proc_day,output_day_p)


    ###  Percentage of patient without any procedure

    npat_day_all<- output_day%>%
        dplyr::filter( periode_group =="all")%>%
        dplyr::filter(concept == "patient")

    proc_no_day <- LocalPatientObservations %>%
        dplyr::filter(concept_type %in% c("PROC-ICD10","PROC-ICD9")) %>%
        dplyr::group_by(GROUP,days_since_admission) %>%
        dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
        dplyr::right_join(npat_day_all[,c("value","GROUP","days_since_admission")],by =c("GROUP","days_since_admission"))%>%
        dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
        dplyr::mutate(nodiag_npat=value-count_pat) %>%
        dplyr::select(- c(value,count_pat)) %>%
        dplyr::rename(value=nodiag_npat) %>%
        dplyr::mutate(concept="PROC-ICD",
                      variable="npat_noproc",
                      periode_group="all")%>%
        data.frame()

    proc_no_day_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9")) %>%
        dplyr::group_by(GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
        dplyr::right_join(npat_day_all[,c("value","GROUP","days_since_admission")],by =c("GROUP","days_since_admission"))%>%
        dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
        dplyr::mutate(nodiag_npat=value-count_pat) %>%
        dplyr::select(- c(value,count_pat)) %>%
        dplyr::rename(value=nodiag_npat) %>%
        dplyr::mutate(concept="PROC-ICD",
                      variable="npat_noproc")%>%
        data.frame()

    output_day=rbind(output_day,proc_no_day,proc_no_day_p)


    ## Laboratory

    ### Average number of laboratory per day


    lab_day <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="LAB-LOINC") %>%
        dplyr::group_by(patient_num,GROUP,days_since_admission) %>%
        dplyr::summarise(count_proc=n()) %>%
        dplyr::right_join(LocalPatientClinicalCourse[,c("patient_num","GROUP","days_since_admission")],by =c("patient_num","GROUP","days_since_admission"))%>%
        dplyr::mutate(count_proc=ifelse(is.na(count_proc),0,count_proc))%>%
        dplyr::group_by(GROUP,days_since_admission) %>%
        dplyr::summarise(value=mean(count_proc,na.rm = TRUE))%>%
        dplyr::mutate(concept="LAB-LOINC",
                      variable="average",
                      periode_group="all")%>%
        data.frame()

    lab_day_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="LAB-LOINC") %>%
        dplyr::group_by(patient_num,GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(count_proc=n()) %>%
        dplyr::right_join(LocalPatientClinicalCourse[,c("patient_num","GROUP","days_since_admission")],by =c("patient_num","GROUP","days_since_admission"))%>%
        dplyr::mutate(count_proc=ifelse(is.na(count_proc),0,count_proc))%>%
        dplyr::group_by(GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(value=mean(count_proc,na.rm = TRUE))%>%
        dplyr::mutate(concept="LAB-LOINC",
                      variable="average")%>%
        data.frame()

    output_day=rbind(output_day,lab_day,lab_day_p)


    ### Percentage of patient without any labs

    lab_no_day <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="LAB-LOINC") %>%
        dplyr::group_by(GROUP,days_since_admission) %>%
        dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
        dplyr::right_join(npat_day_all[,c("value","GROUP","days_since_admission")],by =c("GROUP","days_since_admission"))%>%
        dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
        dplyr::mutate(nodiag_npat=value-count_pat) %>%
        dplyr::select(- c(value,count_pat)) %>%
        dplyr::rename(value=nodiag_npat) %>%
        dplyr::mutate(concept="LAB-LOINC",
                      variable="npat_nolab",
                      periode_group="all")%>%
        data.frame()

    lab_no_day_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="LAB-LOINC") %>%
        dplyr::group_by(GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
        dplyr::right_join(npat_day_all[,c("value","GROUP","days_since_admission")],by =c("GROUP","days_since_admission"))%>%
        dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
        dplyr::mutate(nodiag_npat=value-count_pat) %>%
        dplyr::select(- c(value,count_pat)) %>%
        dplyr::rename(value=nodiag_npat) %>%
        dplyr::mutate(concept="LAB-LOINC",
                      variable="npat_nolab")%>%
        data.frame()

    output_day=rbind(output_day,lab_no_day,lab_no_day_p)


    ## Medication:


    ###Average number of medication per day

    med_day <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::group_by(patient_num,GROUP,days_since_admission) %>%
        dplyr::summarise(count_proc=n()) %>%
        dplyr::right_join(LocalPatientClinicalCourse[,c("patient_num","GROUP","days_since_admission")],by =c("patient_num","GROUP","days_since_admission"))%>%
        dplyr::mutate(count_proc=ifelse(is.na(count_proc),0,count_proc))%>%
        dplyr::filter( days_since_admission >= 0 & days_since_admission <15)%>%
        dplyr::group_by(GROUP,days_since_admission) %>%
        dplyr::summarise(value=mean(count_proc,na.rm = TRUE)) %>%
        dplyr::mutate(concept="MED-CLASS",
                      variable="average",
                      periode_group="all")%>%
        data.frame()

    med_day_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::group_by(patient_num,GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(count_proc=n()) %>%
        dplyr::right_join(LocalPatientClinicalCourse[,c("patient_num","GROUP","days_since_admission")],by =c("patient_num","GROUP","days_since_admission"))%>%
        dplyr::mutate(count_proc=ifelse(is.na(count_proc),0,count_proc))%>%
        dplyr::filter( days_since_admission >= 0 & days_since_admission <15)%>%
        dplyr::group_by(GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(value=mean(count_proc,na.rm = TRUE)) %>%
        dplyr::mutate(concept="LAB-LOINC",
                      variable="average")%>%
        data.frame()

    output_day=rbind(output_day,med_day,med_day_p)


    ### Percentage of patient without any medication

    med_no_day <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::group_by(GROUP,days_since_admission) %>%
        dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
        dplyr::right_join(npat_day_all[,c("value","GROUP","days_since_admission")],by =c("GROUP","days_since_admission"))%>%
        dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
        dplyr::mutate(nodiag_npat=value-count_pat) %>%
        dplyr::select(- c(value,count_pat)) %>%
        dplyr::rename(value=nodiag_npat) %>%
        dplyr::mutate(concept="MED-CLASS",
                      variable="npat_nomed",
                      periode_group="all")%>%
        data.frame()

    med_no_day_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::group_by(GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
        dplyr::right_join(npat_day_all[,c("value","GROUP","days_since_admission")],by =c("GROUP","days_since_admission"))%>%
        dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
        dplyr::mutate(nodiag_npat=value-count_pat) %>%
        dplyr::select(- c(value,count_pat)) %>%
        dplyr::rename(value=nodiag_npat) %>%
        dplyr::mutate(concept="MED-CLASS",
                      variable="npat_nomed")%>%
        data.frame()

    output_day=rbind(output_day,med_no_day,med_no_day_p)


    ##  PaO2 - Distribution of number of Pa02 samples  per groups **

    pao2_day <- LocalPatientObservations %>%
        dplyr::mutate( pao2= if_else(concept_code == loinc_pao2,1,0))  %>%
        dplyr::group_by(patient_num,GROUP,days_since_admission) %>%
        dplyr::summarise(numvalue = sum(pao2))%>%
        dplyr::right_join(LocalPatientClinicalCourse[,c("patient_num","GROUP","days_since_admission")],by =c("patient_num","GROUP","days_since_admission"))%>%
        dplyr::mutate(numvalue=ifelse(is.na(numvalue),0,numvalue))%>%
        dplyr::filter( days_since_admission >= 0 & days_since_admission <30)%>%
        dplyr::group_by(GROUP,days_since_admission) %>%
        dplyr::summarise(value=mean(numvalue,na.rm = TRUE)) %>%
        dplyr::mutate(concept="PA02",
                      variable="average",
                      periode_group="all")%>%
        data.frame()

    pao2_day_p <- LocalPatientObservations %>%
        dplyr::mutate( pao2= if_else(concept_code == loinc_pao2,1,0))  %>%
        dplyr::group_by(patient_num,GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(numvalue = sum(pao2))%>%
        dplyr::right_join(LocalPatientClinicalCourse[,c("patient_num","GROUP","days_since_admission")],by =c("patient_num","GROUP","days_since_admission"))%>%
        dplyr::mutate(numvalue=ifelse(is.na(numvalue),0,numvalue))%>%
        dplyr::filter( days_since_admission >= 0 & days_since_admission <30)%>%
        dplyr::group_by(GROUP,days_since_admission,periode_group) %>%
        dplyr::summarise(value=mean(numvalue,na.rm = TRUE)) %>%
        dplyr::mutate(concept="PA02",
                      variable="average")%>%
        data.frame()

    output_day=rbind(output_day,pao2_day,pao2_day_p)

    ## ARDS (ICD-10) versus More than 1 Pa02 between 5 and 20 days after the admission

    ### For the entire population

    ### check if there Pa02
    test_ADRS_Pa02 <- LocalPatientObservations %>%
        dplyr::filter( days_since_admission >= 0 & concept_code ==loinc_pao2)

    if (nrow(pat_ARDS)==0 | nrow(test_ADRS_Pa02)==0) {

        # create output
        output_sens <- data.frame(matrix(ncol = 5, nrow = 1))
        x <- c("popu", "sensibility", "ppv", "npv","specificity")
        colnames(output_sens) <- x

        output_sens$popu="all"
        output_sens$sensibility=-999
        output_sens$ppv=-999
        output_sens$npv=-999
        output_sens$specificity=-999

    }else {

        ADRS_Pa02 <- LocalPatientObservations %>%
            dplyr::filter( days_since_admission >= 0 & concept_code ==loinc_pao2) %>%
            dplyr:: mutate( ARDS =  if_else(GROUP %in% c("ARDS_sup_49","ARDS_18_49"),1,0),
                            PAO2= if_else(concept_code == loinc_pao2 & ( days_since_admission >= 5 & days_since_admission <= 20),1,0))%>%
            dplyr::group_by(patient_num,ARDS,GROUP)%>%
            dplyr::summarise(PAO2sup1 =if_else(sum(PAO2)>2,1,0))%>%
            unique()

        #on the entire population less than 50
        xtab_ADRS_Pa02 <- table( ADRS_Pa02$PAO2sup1, ADRS_Pa02$ARDS,dnn = c( "more than 1 Pa02","ARDS"))
        colnames(xtab_ADRS_Pa02)=c("NO ARDS",	"ARDS")
        rownames(xtab_ADRS_Pa02)=c("No Pa02 sample between 5-20 days ",	"At least one Pa02 sample between 5-20 days")

        SEN_PaO2_ARDS=caret::sensitivity(as.factor(ADRS_Pa02$PAO2sup1), as.factor(ADRS_Pa02$ARDS), "1")
        PPV_PaO2_ARDS=caret::posPredValue(as.factor(ADRS_Pa02$PAO2sup1), as.factor(ADRS_Pa02$ARDS), "1")
        NPV_PaO2_ARDS=caret::negPredValue(as.factor(ADRS_Pa02$PAO2sup1), as.factor(ADRS_Pa02$ARDS), "0")
        SPE_PaO2_ARDS=caret::specificity(as.factor(ADRS_Pa02$PAO2sup1), as.factor(ADRS_Pa02$ARDS), "0")

        # create output
        output_sens <- data.frame(matrix(ncol = 5, nrow = 1))
        x <- c("popu", "sensibility", "ppv", "npv","specificity")
        colnames(output_sens) <- x

        output_sens$popu="all"
        output_sens$sensibility=SEN_PaO2_ARDS
        output_sens$ppv=PPV_PaO2_ARDS
        output_sens$npv=NPV_PaO2_ARDS
        output_sens$specificity=SPE_PaO2_ARDS


    }


    ### For the  population less than 49 year old

    n_ards_y=num_pat$npat[num_pat$GROUP=="ARDS_18_49"]

    if (length(n_ards_y)==0 | nrow(test_ADRS_Pa02)==0) {
        # create output
        output_sens_y <- data.frame(matrix(ncol = 5, nrow = 1))
        colnames(output_sens_y) <- x

        output_sens_y$popu="less_50"
        output_sens_y$sensibility=-999
        output_sens_y$ppv=-999
        output_sens_y$npv=-999
        output_sens_y$specificity=-999

    }else {

        ADRS_Pa02_y<-ADRS_Pa02  %>% dplyr::filter( GROUP %in% c("ARDS_18_49", "NO_ARDS_18_49", "OTHERS_18_49"))
        xtab_ADRS_Pa02_y <- table(ADRS_Pa02_y$ARDS, ADRS_Pa02_y$PAO2sup1, dnn = c("ARDS", "more than 1 Pa02"))

        xtab_ADRS_Pa02_y  <- table( ADRS_Pa02_y$PAO2sup1,ADRS_Pa02_y$ARDS,dnn = c( "more than 1 Pa02","ARDS"))
        colnames(xtab_ADRS_Pa02_y)=c("NO ARDS",	"ARDS")
        rownames(xtab_ADRS_Pa02_y)=c("No Pa02 sample between 5-20 days ",	"At least one Pa02 sample between 5-20 days")

        SEN_PaO2_ARDS_y=caret::sensitivity(as.factor(ADRS_Pa02_y$PAO2sup1), as.factor(ADRS_Pa02_y$ARDS), "1")
        PPV_PaO2_ARDS_y=caret::posPredValue(as.factor(ADRS_Pa02_y$PAO2sup1), as.factor(ADRS_Pa02_y$ARDS), "1")
        NPV_PaO2_ARDS_y=caret::negPredValue(as.factor(ADRS_Pa02_y$PAO2sup1), as.factor(ADRS_Pa02_y$ARDS), "0")
        SPE_PaO2_ARDS_y=caret::specificity(as.factor(ADRS_Pa02_y$PAO2sup1), as.factor(ADRS_Pa02_y$ARDS), "0")

        # create output
        output_sens_y <- data.frame(matrix(ncol = 5, nrow = 1))
        colnames(output_sens_y) <- x

        output_sens_y$popu="less_50"
        output_sens_y$sensibility=SEN_PaO2_ARDS_y
        output_sens_y$ppv=PPV_PaO2_ARDS_y
        output_sens_y$npv=NPV_PaO2_ARDS_y
        output_sens_y$specificity=SPE_PaO2_ARDS_y

    }



    output_sens=rbind(output_sens,output_sens_y)

    ### age

    age<-  LocalPatientSummary %>%
        dplyr::group_by(GROUP,age_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=age_group) %>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    age_p<-  LocalPatientSummary %>%
        dplyr::group_by(GROUP,age_group,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=age_group) %>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,age,age_p)

    ### sex

    sex<-  LocalPatientSummary %>%
        dplyr::group_by(GROUP,sex) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=sex) %>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    sex_p<-  LocalPatientSummary %>%
        dplyr::group_by(GROUP,sex,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=sex) %>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,sex,sex_p)

    ### output

    pre_hospit<-  LocalPatientObservations %>%
        dplyr::filter(concept_type %in% c("DIAG-ICD10","DIAG-ICD9") & days_since_admission <= limit_d_befor ) %>%
        dplyr::mutate(previous_hospi = "YES") %>%
        dplyr::select( patient_num, GROUP,previous_hospi) %>%
        unique()%>%
        dplyr::right_join(LocalPatientSummary[,c("patient_num","GROUP")],by =c("patient_num","GROUP"))%>%
        dplyr::mutate(previous_hospi = if_else(is.na(previous_hospi),"no_previous_hospi","previous_hospi")) %>%
        dplyr::group_by(GROUP,previous_hospi) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=previous_hospi) %>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    pre_hospit_p<-  LocalPatientObservations %>%
        dplyr::filter(concept_type %in% c("DIAG-ICD10","DIAG-ICD9") & days_since_admission <= limit_d_befor ) %>%
        dplyr::mutate(previous_hospi = "YES") %>%
        dplyr::select( patient_num, GROUP,previous_hospi) %>%
        unique()%>%
        dplyr::right_join(LocalPatientSummary[,c("patient_num","GROUP","periode_group")],by =c("patient_num","GROUP"))%>%
        dplyr::mutate(previous_hospi = if_else(is.na(previous_hospi),"no_previous_hospi","previous_hospi")) %>%
        dplyr::group_by(GROUP,previous_hospi,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=previous_hospi) %>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,pre_hospit,pre_hospit_p)

    ### output

    rehosp<-  duration %>%
        dplyr::mutate(rehosp= if_else(number_rehospit > 0, "rehospit", "no_rehospit"))%>%
        dplyr::group_by(GROUP,rehosp) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=rehosp) %>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    rehosp_p<-  duration %>%
        dplyr::mutate(rehosp= if_else(number_rehospit > 0, "rehospit", "no_rehospit"))%>%
        dplyr::group_by(GROUP,rehosp,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=rehosp) %>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,rehosp,rehosp_p)

    ### elixhauser_score


    comorb_names_elix <- get_quan_elix_names()
    comorbs_elix <- as.vector(comorb_names_elix$Abbreviation)

    comorb_elix_min7 <- map_char_elix_codes(
        df = LocalPatientObservations,
        comorb_names = comorb_names_elix,
        t1 = -365,
        t2 = limit_d_befor,
        map_type = 'elixhauser',
        truncate = TRUE
    )

    comorb_elix_min7 <- comorb_elix_min7$index_scores %>%
        rename('elixhauser_score' = van_walraven_score)%>%
        mutate(patient_num = as.integer(patient_num))

    output_elix_min7 <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_min7, by = "patient_num")%>%
        dplyr::group_by(GROUP) %>%
        dplyr::summarise(elix_mean_min7=mean(elixhauser_score, na.rm = TRUE),
                         elix_std_min7=sd(elixhauser_score, na.rm = TRUE))%>%
        pivot_longer(c(elix_mean_min7,elix_std_min7),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="patient")%>%
        dplyr::mutate(periode_group="all")%>%
        data.frame()

    output_elix_min7_p <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_min7, by = "patient_num")%>%
        dplyr::group_by(GROUP,periode_group) %>%
        dplyr::summarise(elix_mean_min7=mean(elixhauser_score, na.rm = TRUE),
                         elix_std_min7=sd(elixhauser_score, na.rm = TRUE))%>%
        pivot_longer(c(elix_mean_min7,elix_std_min7),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,output_elix_min7,output_elix_min7_p)

    comorb_elix_90 <- map_char_elix_codes(
        df = LocalPatientObservations,
        comorb_names = comorb_names_elix,
        t1 = -365,
        t2 = 90,
        map_type = 'elixhauser',
        truncate = TRUE
    )

    comorb_elix_90 <- comorb_elix_90$index_scores %>%
        rename('elixhauser_score' = van_walraven_score)%>%
        mutate(patient_num = as.integer(patient_num))

    output_elix_90 <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_90, by = "patient_num")%>%
        dplyr::group_by(GROUP) %>%
        dplyr::summarise(elix_mean_90=mean(elixhauser_score, na.rm = TRUE),
                         elix_std_90=sd(elixhauser_score, na.rm = TRUE))%>%
        pivot_longer(c(elix_mean_90,elix_std_90),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="patient")%>%
        dplyr::mutate(periode_group="all")%>%
        data.frame()

    output_elix_90_p <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_90, by = "patient_num")%>%
        dplyr::group_by(GROUP,periode_group) %>%
        dplyr::summarise(elix_mean_90=mean(elixhauser_score, na.rm = TRUE),
                         elix_std_90=sd(elixhauser_score, na.rm = TRUE))%>%
        pivot_longer(c(elix_mean_90,elix_std_90),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,output_elix_90,output_elix_90_p)


    ### elixhauser_score by commorbdities

    output_elix_min7_score <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_min7, by = "patient_num")%>%
        dplyr::group_by(GROUP) %>%
        dplyr::summarise(across(comorbs_elix, ~sum(.x, na.rm = TRUE)))%>%
        pivot_longer(comorbs_elix,names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="elix",
                      periode_group="all",
                      time="before")%>%
        data.frame()

    output_elix_min7_score_p <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_min7, by = "patient_num")%>%
        dplyr::group_by(GROUP,periode_group) %>%
        dplyr::summarise(across(comorbs_elix, ~sum(.x, na.rm = TRUE)))%>%
        pivot_longer(comorbs_elix,names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="elix",
                      time="before")%>%
        data.frame()

    output_elix_90_score <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_90, by = "patient_num")%>%
        dplyr::group_by(GROUP) %>%
        dplyr::summarise(across(comorbs_elix, ~sum(.x, na.rm = TRUE)))%>%
        pivot_longer(comorbs_elix,names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="elix",
                      periode_group="all",
                      time="all")%>%
        data.frame()

    output_elix_90_score_p <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_90, by = "patient_num")%>%
        dplyr::group_by(GROUP,periode_group) %>%
        dplyr::summarise(across(comorbs_elix, ~sum(.x, na.rm = TRUE)))%>%
        pivot_longer(comorbs_elix,names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="elix",
                      time="all")%>%
        data.frame()


    out_elix=rbind(output_elix_min7_score,output_elix_min7_score_p,output_elix_90_score,output_elix_90_score_p)

    ### manage ICD10 and 9

    ICD=rbind(ICD10[,c("code","long_desc","major","sub_chapter","chapter")],
              ICD9[,c("code","long_desc","major","sub_chapter","chapter")])


    diagnosis<-LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("DIAG-ICD10","DIAG-ICD9" )) %>%
        dplyr::left_join( ICD10[,c("code","long_desc","major","sub_chapter","chapter")],by =c( "concept_code" = "code"))%>%
        data.frame()


    # all time
    out_iCD_all <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="all")%>%
        data.frame()

    out_iCD_all_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type  %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="all")%>%
        data.frame()


    # before time
    out_iCD_before <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::filter( days_since_admission <= limit_d_befor) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="before")%>%
        data.frame()

    out_iCD_before_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type  %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::filter( days_since_admission <= limit_d_befor) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="before")%>%
        data.frame()

    # after time
    out_iCD_after <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="after")%>%
        data.frame()

    out_iCD_after_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type  %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="after")%>%
        data.frame()

    # after 90 days after admission time
    out_iCD_after90 <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::filter( days_since_admission >= 90) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="after_90")%>%
        data.frame()

    out_iCD_after90_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type  %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
        dplyr::filter( days_since_admission >= 90) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="after_90")%>%
        data.frame()

    out_iCD=rbind(out_iCD_all,out_iCD_before,out_iCD_after,out_iCD_after90,out_iCD_all_p,out_iCD_before_p,out_iCD_after_p,out_iCD_after90_p)


    ### add phenotype

    diag_pheno<-LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("DIAG-ICD10","DIAG-ICD9" )) %>%
        dplyr::left_join( pheno_ICD[,c("ICD","ICD.String","PheCode","Phenotype")],by =c( "concept_code" = "ICD"))%>%
        data.frame()


    # all time
    out_phe_all <-diag_pheno %>%
        dplyr::group_by(PheCode, GROUP) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=PheCode )%>%
        dplyr::mutate(concept="PheCode",
                      periode_group="all",
                      time="all")%>%
        data.frame()

    out_phe_all_p <-diag_pheno %>%
        dplyr::group_by(PheCode, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=PheCode )%>%
        dplyr::mutate(concept="PheCode",
                      time="all")%>%
        data.frame()


    # before time
    out_phe_before <-diag_pheno %>%
        dplyr::group_by(PheCode, GROUP) %>%
        dplyr::filter( days_since_admission <= limit_d_befor) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=PheCode )%>%
        dplyr::mutate(concept="PheCode",
                      periode_group="all",
                      time="before")%>%
        data.frame()

    out_phe_before_p <-diag_pheno %>%
        dplyr::group_by(PheCode, GROUP,periode_group) %>%
        dplyr::filter( days_since_admission <= limit_d_befor) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=PheCode )%>%
        dplyr::mutate(concept="PheCode",
                      time="before")%>%
        data.frame()

    # after time
    out_phe_after <-diag_pheno %>%
        dplyr::group_by(PheCode, GROUP) %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=PheCode )%>%
        dplyr::mutate(concept="PheCode",
                      periode_group="all",
                      time="after")%>%
        data.frame()

    out_phe_after_p <-diag_pheno %>%
        dplyr::group_by(PheCode, GROUP,periode_group) %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=PheCode )%>%
        dplyr::mutate(concept="PheCode",
                      time="after")%>%
        data.frame()

    # after 90 days after admission time
    out_phe_after90 <-diag_pheno %>%
        dplyr::group_by(PheCode, GROUP) %>%
        dplyr::filter( days_since_admission >= 90) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=PheCode )%>%
        dplyr::mutate(concept="PheCode",
                      periode_group="all",
                      time="after_90")%>%
        data.frame()

    out_phe_after90_p <-diag_pheno %>%
        dplyr::group_by(PheCode, GROUP,periode_group) %>%
        dplyr::filter( days_since_admission >= 90) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=PheCode )%>%
        dplyr::mutate(concept="PheCode",
                      time="after_90")%>%
        data.frame()


    out_phe=rbind(out_phe_all,out_phe_all_p,
                  out_phe_before,out_phe_before_p,
                  out_phe_after,out_phe_after_p,
                  out_phe_after90,out_phe_after90_p)

    ### add phenotype

    diag_comp<-LocalPatientObservations %>%
        dplyr::filter( concept_type == "DIAG-ICD10") %>%
        dplyr::left_join( comp_class,by =c( "concept_code" = "ICD10"))%>%
        data.frame()


    # after time
    out_comp_after <-diag_comp %>%
        dplyr::group_by(complication_class, GROUP) %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=complication_class )%>%
        dplyr::mutate(concept="complication_class",
                      periode_group="all",
                      time="after")%>%
        data.frame()

    out_comp_after_p <-diag_comp %>%
        dplyr::group_by(complication_class, GROUP,periode_group) %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=complication_class )%>%
        dplyr::mutate(concept="complication_class",
                      time="after")%>%
        data.frame()

    # after 90 days after admission time
    out_comp_after90 <-diag_comp %>%
        dplyr::group_by(complication_class, GROUP) %>%
        dplyr::filter( days_since_admission >= 90) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=complication_class )%>%
        dplyr::mutate(concept="complication_class",
                      periode_group="all",
                      time="after_90")%>%
        data.frame()

    out_comp_after90_p <-diag_comp %>%
        dplyr::group_by(complication_class, GROUP,periode_group) %>%
        dplyr::filter( days_since_admission >= 90) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=complication_class )%>%
        dplyr::mutate(concept="complication_class",
                      time="after_90")%>%
        data.frame()


    out_comp=rbind(out_comp_after,out_comp_after_p,
                   out_comp_after90,out_comp_after90_p)


    ### medication
    # all time
    out_med_all <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="all")%>%
        data.frame()

    out_med_all_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="all")%>%
        data.frame()

    # before time
    out_med_before <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::filter( days_since_admission < 0) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="before")%>%
        data.frame()

    out_med_before_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::filter( days_since_admission < 0) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="before")%>%
        data.frame()

    # after time
    out_med_after <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="after")%>%
        data.frame()

    out_med_after_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="after")%>%
        data.frame()

    # after 90 days after admission time
    out_med_after90 <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::filter( days_since_admission >= 90) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="after_90")%>%
        data.frame()

    out_med_after90_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="MED-CLASS") %>%
        dplyr::filter( days_since_admission >= 90) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="after_90")%>%
        data.frame()

    out_med=rbind(out_med_all,out_med_before,out_med_after,out_med_after90,out_med_all_p,out_med_before_p,out_med_after_p,out_med_after90_p)

    ### procedure

    # all time
    out_proc_all <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="PROC-ICD10") %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="all")%>%
        data.frame()

    out_proc_all_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="PROC-ICD10")  %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="all")%>%
        data.frame()

    # before time
    out_proc_before <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="PROC-ICD10") %>%
        dplyr::filter( days_since_admission < 0) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="before")%>%
        data.frame()

    out_proc_before_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="PROC-ICD10")  %>%
        dplyr::filter( days_since_admission < 0) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="before")%>%
        data.frame()

    # after time
    out_proc_after <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="PROC-ICD10")  %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="after")%>%
        data.frame()

    out_proc_after_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="PROC-ICD10") %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="after")%>%
        data.frame()

    # after 90 days after admission time
    out_proc_after90 <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="PROC-ICD10")  %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="after_90")%>%
        data.frame()

    out_proc_after90_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="PROC-ICD10")  %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="after_90")%>%
        data.frame()

    out_proc=rbind(out_proc_all,out_proc_before,out_proc_after,out_proc_after90,out_proc_all_p,out_proc_before_p,out_proc_after_p,out_proc_after90_p)

    out_proc_diag_med=rbind(out_iCD,out_proc,out_med,out_elix,out_phe,out_comp)

    ###  Lab -- output

    output_lab_all <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="LAB-LOINC") %>%
        dplyr::group_by(GROUP, days_since_admission, concept_code) %>%
        dplyr::summarise(npat=n_distinct(patient_num),
                         mean_value=mean(value, na.rm = TRUE),
                         std_value=sd(value, na.rm = TRUE),
                         mean_log_value = mean(log(value + 0.5), na.rm = TRUE),
                         std_log_value = sd(log(value + 0.5), na.rm = TRUE))%>%
        dplyr::mutate(periode_group="all")%>%
        data.frame()

    output_lab_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type =="LAB-LOINC") %>%
        dplyr::group_by(GROUP, days_since_admission, concept_code,periode_group) %>%
        dplyr::summarise(npat=n_distinct(patient_num),
                         mean_value=mean(value, na.rm = TRUE),
                         std_value=sd(value, na.rm = TRUE),
                         mean_log_value = mean(log(value + 0.5), na.rm = TRUE),
                         std_log_value = sd(log(value + 0.5), na.rm = TRUE))%>%
        data.frame()

    output_lab=rbind(output_lab_all,output_lab_p)

    ### outcomes

    out_out <- LocalPatientSummary %>%
        dplyr::inner_join(duration[,c("patient_num","length_duration_1")], by ="patient_num") %>%
        dplyr::mutate(dead_28=ifelse(deceased == 1 & floor(as.numeric(difftime(death_date, admission_date, units = "days")))<28,1,0),
                      dead_90=ifelse(deceased == 1 & floor(as.numeric(difftime(death_date, admission_date, units = "days")))<90,1,0),
                      out_hospit_28=ifelse(length_duration_1<28,1,0),
                      out_hospit_90=ifelse(length_duration_1<90,1,0))%>%
        dplyr::group_by(GROUP)%>%
        dplyr::summarise(dead=sum(deceased),
                         dead_less_28=sum(dead_28),
                         dead_less_90=sum(dead_90),
                         out_hospit_28=sum(out_hospit_28),
                         out_hospit_90=sum(out_hospit_90))%>%
        dplyr::right_join(num_pat,by = "GROUP")%>%
        dplyr::select(GROUP,dead,dead_less_28,dead_less_90,out_hospit_28,out_hospit_90)%>%
        pivot_longer(c(dead,dead_less_28,dead_less_90,out_hospit_28,out_hospit_90),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="patient",
                      periode_group="all")%>%
        data.frame()

    out_out_p <- LocalPatientSummary %>%
        dplyr::inner_join(duration[,c("patient_num","length_duration_1")], by ="patient_num") %>%
        dplyr::mutate(dead_28=ifelse(deceased == 1 & floor(as.numeric(difftime(death_date, admission_date, units = "days")))<28,1,0),
                      dead_90=ifelse(deceased == 1 & floor(as.numeric(difftime(death_date, admission_date, units = "days")))<90,1,0),
                      out_hospit_28=ifelse(length_duration_1<28,1,0),
                      out_hospit_90=ifelse(length_duration_1<90,1,0))%>%
        dplyr::group_by(GROUP,periode_group)%>%
        dplyr::summarise(dead=sum(deceased),
                         dead_less_28=sum(dead_28),
                         dead_less_90=sum(dead_90),
                         out_hospit_28=sum(out_hospit_28),
                         out_hospit_90=sum(out_hospit_90))%>%
        dplyr::right_join(num_pat,by = "GROUP")%>%
        dplyr::select(GROUP,periode_group,dead,dead_less_28,dead_less_90,out_hospit_28,out_hospit_90)%>%
        pivot_longer(c(dead,dead_less_28,dead_less_90,out_hospit_28,out_hospit_90),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,out_out,out_out_p)

    ### status

    patnum_day=unique(LocalPatientSummary[,c("siteid","patient_num","GROUP","periode_group")])

    patnum_day=patnum_day[rep(1:nrow(patnum_day),times = 3),]
    patnum_day$days_since_admission <- rep(c(6,27,89), each=nrow(patnum_day)/3)

    status7 <- LocalPatientClinicalCourse %>%
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","siteid","GROUP","periode_group"))%>%
        dplyr::filter(days_since_admission == "6" )%>%
        dplyr::mutate(status_at_d = ifelse(is.na(deceased),"no_data",
                                           ifelse(deceased==1,"death",
                                                  ifelse(in_hospital==1,"in_hospital","out_hospital"))))%>%
        dplyr::group_by(GROUP,status_at_d) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=status_at_d) %>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::mutate(concept="status_at_7")%>%
        data.frame()

    status28 <- LocalPatientClinicalCourse %>%
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","siteid","GROUP","periode_group"))%>%
        dplyr::filter(days_since_admission == "27" )%>%
        dplyr::mutate(status_at_d = ifelse(is.na(deceased),"no_data",
                                           ifelse(deceased==1,"death",
                                                  ifelse(in_hospital==1,"in_hospital","out_hospital"))))%>%
        dplyr::group_by(GROUP,status_at_d) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=status_at_d) %>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::mutate(concept="status_at_28")%>%
        data.frame()

    status90 <- LocalPatientClinicalCourse %>%
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","siteid","GROUP","periode_group"))%>%
        dplyr::filter(days_since_admission == "89" )%>%
        dplyr::mutate(status_at_d = ifelse(is.na(deceased),"no_data",
                                           ifelse(deceased==1,"death",
                                                  ifelse(in_hospital==1,"in_hospital","out_hospital"))))%>%
        dplyr::group_by(GROUP,status_at_d) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=status_at_d) %>%
        dplyr::mutate(periode_group="all")%>%
        dplyr::mutate(concept="status_at_90")%>%
        data.frame()


    status7_p <- LocalPatientClinicalCourse %>%
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","siteid","GROUP","periode_group"))%>%
        dplyr::filter(days_since_admission == "6" )%>%
        dplyr::mutate(status_at_d = ifelse(is.na(deceased),"no_data",
                                           ifelse(deceased==1,"death",
                                                  ifelse(in_hospital==1,"in_hospital","out_hospital"))))%>%
        dplyr::group_by(GROUP,status_at_d,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=status_at_d) %>%
        dplyr::mutate(concept="status_at_7")%>%
        data.frame()

    status28_p <- LocalPatientClinicalCourse %>%
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","siteid","GROUP","periode_group"))%>%
        dplyr::filter(days_since_admission == "27" )%>%
        dplyr::mutate(status_at_d = ifelse(is.na(deceased),"no_data",
                                           ifelse(deceased==1,"death",
                                                  ifelse(in_hospital==1,"in_hospital","out_hospital"))))%>%
        dplyr::group_by(GROUP,status_at_d,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=status_at_d) %>%
        dplyr::mutate(concept="status_at_28")%>%
        data.frame()

    status90_p <- LocalPatientClinicalCourse %>%
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","siteid","GROUP","periode_group"))%>%
        dplyr::filter(days_since_admission == "89" )%>%
        dplyr::mutate(status_at_d = ifelse(is.na(deceased),"no_data",
                                           ifelse(deceased==1,"death",
                                                  ifelse(in_hospital==1,"in_hospital","out_hospital"))))%>%
        dplyr::group_by(GROUP,status_at_d,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=status_at_d) %>%
        dplyr::mutate(concept="status_at_90")%>%
        data.frame()


    output_gen=rbind(output_gen,status7,status7_p,status28,status28_p,status90,status90_p)


    ## ========================================
    ## PART 4 : Analysis
    ## ========================================


    siteid=unique(LocalPatientSummary$siteid)

    output_day$siteid=siteid
    output_gen$siteid=siteid
    output_lab$siteid=siteid
    out_proc_diag_med$siteid=siteid
    output_sens$siteid=siteid

    ## manage std with 1 patient
    output_lab <- output_lab %>%
        dplyr::mutate(std_value = ifelse(npat == 1, 0, std_value),
                      std_log_value = ifelse(npat == 1, 0, std_log_value))

    var_gen_num = c("number","n_severe" ,"26to49","50to69","70to79","80plus","18to25","female","male","no_previous_hospi","previous_hospi","no_rehospit","rehospit","dead","dead_less_28","dead_less_90","out_hospit_28","out_hospit_90","in_hospital", "death", "out_hospital","no_data")

    var_gen_ag =c("mean_freq","std_freq", "mean_los", "std_los")

    ## manage obfuscation

    ## add obfuscation

    obfusc=data.frame("siteid" = siteid, "truefalse" = obfuscation, "value" = obfuscationThreshord)

    if(obfusquation){

        group_less_obf <- output_gen %>%
            dplyr::filter(variable =="number" &  value< obfuscationThreshord)%>%
            dplyr::select(GROUP,periode_group)

        output_gen <- output_gen %>%
            dplyr::mutate(
                value = ifelse(
                    test = variable %in% var_gen_num &  value < obfuscationThreshord & value != 0,
                    yes = obfuscationValue, no = value),
                value = ifelse(
                    test = (GROUP %in% group_less_obf$GROUP) & (periode_group %in% group_less_obf$periode_group) & (variable %in% var_gen_ag),
                    yes = obfuscationValue, no = value))

        output_day <- output_day %>%
            dplyr::mutate(
                value = ifelse(
                    test = value < obfuscationThreshord & value != 0 & variable == "number",
                    yes = obfuscationValue, no = value))

        out_proc_diag_med <- out_proc_diag_med %>%
            dplyr::mutate(
                value = ifelse(
                    test = value < obfuscationThreshord & value != 0,
                    yes = obfuscationValue, no = value))

        output_lab <- output_lab %>%
            dplyr::mutate(
                npat = ifelse(
                    test = npat < obfuscationThreshord & npat != 0,
                    yes = obfuscationValue, no = npat),
                mean_value = ifelse(
                    test = npat < obfuscationThreshord & npat != 0,
                    yes = obfuscationValue, no = mean_value),
                std_value= ifelse(
                    test = npat < obfuscationThreshord & npat != 0,
                    yes = obfuscationValue, no = std_value),
                mean_log_value= ifelse(
                    test = npat < obfuscationThreshord & npat != 0,
                    yes = obfuscationValue, no = mean_log_value),
                std_log_value= ifelse(
                    test = npat < obfuscationThreshord & npat != 0,
                    yes = obfuscationValue, no = std_log_value))

    }


    write.csv(output_gen, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_ards_young_gen",".csv")), row.names = FALSE, na = "")
    write.csv(output_day, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_ards_young_day",".csv")), row.names = FALSE, na = "")
    write.csv(output_lab, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_ards_young_lab",".csv")), row.names = FALSE, na = "")
    write.csv(out_proc_diag_med, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_ards_young_proc_diag_med",".csv")), row.names = FALSE, na = "")
    write.csv(obfusc, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_obfusc",".csv")), row.names = FALSE, na = "")
    write.csv(output_sens, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_sens",".csv")), row.names = FALSE, na = "")

}

