
#' Runs the analytic workflow for the ards project
#'
#' @keywords 4CE
#' @return
#' @export
#' @import dplyr tidyr stringr icd caret DT tidyverse icd.data metafor

runAnalysis <- function() {

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

    ## obfuscation
    obfuscationValue <- -99

    obfuscationThreshord= as.numeric(FourCePhase2.1Data::getObfuscation(toupper(currSiteId)))

    if(obfuscationThreshord==0){
      obfuscation = FALSE
    } else {
      obfuscation = TRUE
    }

    ### missing value
    miss_value= -999

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

    # ## ICD
    # ICD_comor_comp <- read_csv2(paste0(path_pack,'/data-raw/ICD_comor_comp.csv'))
    #
    # ## load elix
    #  elix_score= read_csv2(paste0(path_pack,'/data-raw/elix_score.csv'))

    # pheno_ICD=FourCePhase2.1ards::pheno_ICD
    # comp_class=FourCePhase2.1ards::comp_class
    # sev_proc_icd10=FourCePhase2.1ards::sev_proc_icd10
    # lab_mapping=FourCePhase2.1ards::lab_mapping
    # med_code= FourCePhase2.1ards::med_code
    # data(pheno_ICD)
    # data(comp_class)
    # data(sev_proc_icd10)
    # data(lab_mapping)
    # data(med_code)

    ## input ###
    med_severe = c("SIANES","SICARDIAC")
    loinc_pao2 = "2703-7"
    ARDS=c('J80','518.82','518')
    proc_cpt= c("CPT4:31500","CPT4:94656","CPT4:94657","CPT4:94002")

    # P1 : january - july
    # P2 : august - december

    start_date_p1=as.POSIXct(as.Date("2020-01-01"), tz = Sys.timezone())
    end_date_p1=as.POSIXct(as.Date("2020-07-01"), tz = Sys.timezone())
    start_date_p2=as.POSIXct(as.Date("2020-12-31"), tz = Sys.timezone())

    # date of inclusion
    last_date_inclusion=as.POSIXct(as.Date("2020-12-31"), tz = Sys.timezone())


    ### limit before / after
    limit_d_befor= -14
    limit_d_after= 90

    message("load doc and reformat=> OK")
    ## ========================================
    ## PART 3: Group selection
    ## ========================================

    #remove patient admitted after the last date of inclusion
    LocalPatientSummary <- LocalPatientSummary %>%
        dplyr::mutate(admission_date = as.POSIXct((as.Date(admission_date)))) %>%
        dplyr::filter( admission_date <= last_date_inclusion)%>%
        data.frame()

    LocalPatientClinicalCourse <- LocalPatientClinicalCourse %>%
        dplyr::filter( patient_num %in% unique(LocalPatientSummary$patient_num))%>%
        data.frame()

    LocalPatientObservations <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% unique(LocalPatientSummary$patient_num))%>%
        data.frame()

    ## select patient by age
    LocalPatientSummary <- LocalPatientSummary %>%
        dplyr::filter( age_group %in% c("18to25","26to49","50to69","70to79","80plus")) %>%
        dplyr::mutate( age_group_spe = dplyr::case_when(
            age_group %in% c("18to25","26to49")~ "18_49",
            age_group %in% c("50to69","70to79","80plus")~ "sup_49"),
            periode_group = dplyr::case_when(
                admission_date >= start_date_p1 &  admission_date < end_date_p1 ~ "P1",
                admission_date >= end_date_p1 &  admission_date <= start_date_p2 ~ "P2"))%>%
        data.frame()

    ### patients with ARDS
    pat_ARDS <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::filter(grepl(paste(ARDS,collapse="|"), concept_code) & days_since_admission >=0) %>%
        dplyr::select( patient_num)%>%
        dplyr::distinct()%>%
        data.frame()

    ### patients with ventilation
    pat_vent <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::filter((concept_type == "PROC-ICD10"
                      & days_since_admission >=0
                      & concept_code %in% FourCePhase2.1ards::sev_proc_icd10$code)
                      |(concept_type == "PROC-CPT"
                        & days_since_admission >=0
                        & concept_code %in% proc_cpt))%>%
        dplyr::select(patient_num)%>%
        dplyr::distinct()%>%
        data.frame()

    ### patients with severe medication
    pat_med_severe <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::filter(concept_code %in% c("SIANES","SICARDIAC") & days_since_admission >=0) %>%
        dplyr::select( patient_num)%>%
        dplyr::distinct()%>%
        data.frame()

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
                          ARDS==0 & (MED_SEVERE==1 | PROC_SEVERE==1) & age_group_spe== "sup_49" ~ "OTHERS_sup_49"))%>%
        data.frame()

    LocalPatientObservations <- LocalPatientObservations %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::left_join(LocalPatientSummary[,c("patient_num","periode_group","GROUP","ARDS","PROC_SEVERE","MED_SEVERE")])%>%
        data.frame()

    LocalPatientClinicalCourse <- LocalPatientClinicalCourse %>%
        dplyr::filter( patient_num %in% LocalPatientSummary$patient_num) %>%
        dplyr::left_join(LocalPatientSummary[,c("patient_num","periode_group","GROUP","ARDS","PROC_SEVERE","MED_SEVERE")])%>%
        dplyr::mutate(calendar_date = as.POSIXct((as.Date(calendar_date))))%>%
        data.frame()


    #remove patient without any diagnostics


    ### number patient per groups
    output_temp <- LocalPatientSummary %>%
      dplyr::group_by( GROUP) %>%
      dplyr::summarise(npat=n_distinct(patient_num)) %>%
      dplyr::rename(value=npat)%>%
      dplyr::mutate(concept="patient",
                    variable="number",
                    periode_group="all")%>%
      data.frame()

    ### number patient without diagnostics
    diag_no_day <- LocalPatientObservations %>%
      dplyr::filter(concept_type %in%  c("DIAG-ICD10","DIAG-ICD9")) %>%
      dplyr::group_by(GROUP) %>%
      dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
      dplyr::right_join(output_temp[,c("value","GROUP")],by =c("GROUP"))%>%
      dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
      dplyr::mutate(nodiag_npat=value-count_pat) %>%
      dplyr::select(- c(value,count_pat)) %>%
      dplyr::rename(value=nodiag_npat) %>%
      dplyr::mutate(concept="DIAG-ICD",
                    variable="npat_nodiag",
                    periode_group="all")%>%
      data.frame()

    ### remove patients without diag from LocalPatientObservations ,LocalPatientSummary, LocalPatientClinicalCourse
    pat_withday <- LocalPatientObservations %>%
      dplyr::filter(concept_type %in%  c("DIAG-ICD10","DIAG-ICD9"))%>%
      dplyr::select(patient_num)%>%
      unique()%>%
      data.frame()

    pat_nodiag=setdiff(LocalPatientSummary$patient_num,pat_withday$patient_num)

   if (!length(pat_nodiag) == 0) {

     LocalPatientSummary <- LocalPatientSummary %>%
       dplyr::filter(! patient_num %in% pat_nodiag )%>%
       data.frame()

     LocalPatientObservations <- LocalPatientObservations %>%
       dplyr::filter(! patient_num %in% pat_nodiag )%>%
       data.frame()

     LocalPatientClinicalCourse <- LocalPatientClinicalCourse %>%
       dplyr::filter(! patient_num %in% pat_nodiag )%>%
       data.frame()

   }

    message( paste("ncol LocalPatientSummary = ",ncol(LocalPatientSummary)))
    message( paste("nrow LocalPatientSummary = ",nrow(LocalPatientSummary)))

    print(table(LocalPatientSummary$GROUP))
    print(table(LocalPatientSummary$GROUP,LocalPatientSummary$periode_group))

    message("Group selection => OK")


    ## ========================================
    ## PART 4 : Analysis
    ## ========================================

    ## patient by groups
    num_pat <- LocalPatientSummary %>%
        dplyr::group_by( GROUP) %>%
        dplyr::summarise(npat=n_distinct(patient_num))%>%
        data.frame()

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
        dplyr::distinct()%>%
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
                         delay_hospi_rehosp= delay_1rehospit_1admin-length_duration_1)%>%
        data.frame()


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

    # diag_no_day <- LocalPatientObservations %>%
    #     dplyr::filter(concept_type %in%  c("DIAG-ICD10","DIAG-ICD9")) %>%
    #     dplyr::group_by(GROUP) %>%
    #     dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
    #     dplyr::right_join(output[,c("value","GROUP")],by =c("GROUP"))%>%
    #     dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
    #     dplyr::mutate(nodiag_npat=value-count_pat) %>%
    #     dplyr::select(- c(value,count_pat)) %>%
    #     dplyr::rename(value=nodiag_npat) %>%
    #     dplyr::mutate(concept="DIAG-ICD",
    #                   variable="npat_nodiag",
    #                   periode_group="all")%>%
    #     data.frame()
    #
    # diag_no_day_p <- LocalPatientObservations %>%
    #     dplyr::filter( concept_type %in% c("DIAG-ICD10","DIAG-ICD9")) %>%
    #     dplyr::group_by(GROUP,periode_group) %>%
    #     dplyr::summarise(count_pat=n_distinct(patient_num)) %>%
    #     dplyr::right_join(output[,c("value","GROUP")],by =c("GROUP"))%>%
    #     dplyr::mutate(count_pat=ifelse(is.na(count_pat),0,count_pat))%>%
    #     dplyr::mutate(nodiag_npat=value-count_pat) %>%
    #     dplyr::select(- c(value,count_pat)) %>%
    #     dplyr::rename(value=nodiag_npat) %>%
    #     dplyr::mutate(concept="DIAG-ICD",
    #                   variable="npat_nodiag")%>%
    #     data.frame()

    # output_gen=rbind(output_gen,diag_no_day,diag_no_day_p)
    output_gen=rbind(output_gen,diag_no_day)

    ## Procedure:

    ### Average number of procedure per day

    proc_day <- LocalPatientObservations %>%
        dplyr::filter(concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT")) %>%
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
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT")) %>%
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
        dplyr::filter(concept == "patient")%>%
        data.frame()

    proc_no_day <- LocalPatientObservations %>%
        dplyr::filter(concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT")) %>%
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
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT")) %>%
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

    ### check if there is enough patient for sensibility and specificity analysis

    ADRS_Pa02 <- LocalPatientObservations %>%
        dplyr::filter( days_since_admission >= 0 & concept_code ==loinc_pao2) %>%
        dplyr:: mutate( ARDS =  if_else(GROUP %in% c("ARDS_sup_49","ARDS_18_49"),1,0),
                        PAO2= if_else(concept_code == loinc_pao2 & ( days_since_admission >= 5 & days_since_admission <= 20),1,0))%>%
        dplyr::group_by(patient_num,ARDS,GROUP)%>%
        dplyr::summarise(PAO2sup1 =if_else(sum(PAO2)>2,1,0))%>%
        dplyr::distinct()%>%
        data.frame()

    if (nrow(pat_ARDS)==0 |length(unique(ADRS_Pa02$PAO2sup1))<2) {

        # create output
        output_sens <- data.frame(matrix(ncol = 5, nrow = 1))
        colnames(output_sens)  <- c("popu", "sensibility", "ppv", "npv","specificity")

        output_sens$popu="all"
        output_sens$sensibility=-999
        output_sens$ppv=-999
        output_sens$npv=-999
        output_sens$specificity=-999

    }else {

        #on the entire population less than 50
        xtab_ADRS_Pa02 <- table(ADRS_Pa02$PAO2sup1, ADRS_Pa02$ARDS,dnn = c( "more than 1 Pa02","ARDS"))
        colnames(xtab_ADRS_Pa02)=c("NO ARDS","ARDS")
        rownames(xtab_ADRS_Pa02)=c("No Pa02 sample between 5-20 days ",	"At least one Pa02 sample between 5-20 days")

        SEN_PaO2_ARDS=caret::sensitivity(as.factor(ADRS_Pa02$PAO2sup1), as.factor(ADRS_Pa02$ARDS), "1")
        PPV_PaO2_ARDS=caret::posPredValue(as.factor(ADRS_Pa02$PAO2sup1), as.factor(ADRS_Pa02$ARDS), "1")
        NPV_PaO2_ARDS=caret::negPredValue(as.factor(ADRS_Pa02$PAO2sup1), as.factor(ADRS_Pa02$ARDS), "0")
        SPE_PaO2_ARDS=caret::specificity(as.factor(ADRS_Pa02$PAO2sup1), as.factor(ADRS_Pa02$ARDS), "0")

        # create output
        output_sens <- data.frame(matrix(ncol = 5, nrow = 1))
        colnames(output_sens) <- c("popu", "sensibility", "ppv", "npv","specificity")

        output_sens$popu="all"
        output_sens$sensibility=SEN_PaO2_ARDS
        output_sens$ppv=PPV_PaO2_ARDS
        output_sens$npv=NPV_PaO2_ARDS
        output_sens$specificity=SPE_PaO2_ARDS

    }


    ### For the  population less than 49 year old

    n_ards_y=num_pat$npat[num_pat$GROUP=="ARDS_18_49"]

    ADRS_Pa02_y<-ADRS_Pa02  %>%
        dplyr::filter( GROUP %in% c("ARDS_18_49", "NO_ARDS_18_49", "OTHERS_18_49"))%>%
        dplyr::distinct()%>%
        data.frame()

    if (length(n_ards_y)==0 | length(unique(ADRS_Pa02_y$PAO2sup1))<2) {
        # create output
        output_sens_y <- data.frame(matrix(ncol = 5, nrow = 1))
        colnames(output_sens_y) <- c("popu", "sensibility", "ppv", "npv","specificity")

        output_sens_y$popu="less_50"
        output_sens_y$sensibility=-999
        output_sens_y$ppv=-999
        output_sens_y$npv=-999
        output_sens_y$specificity=-999

    }else {

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
        colnames(output_sens_y) <- c("popu", "sensibility", "ppv", "npv","specificity")

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

    pre_hospit_tmp<-  LocalPatientObservations %>%
      dplyr::filter(concept_type %in% c("DIAG-ICD10","DIAG-ICD9") & days_since_admission <= limit_d_befor) %>%
      dplyr::mutate(previous_hospi = "YES")%>%
      dplyr::select( patient_num, GROUP,previous_hospi) %>%
      dplyr::distinct()


    pre_hospit<-  LocalPatientObservations %>%
        dplyr::filter(concept_type %in% c("DIAG-ICD10","DIAG-ICD9") & days_since_admission <= limit_d_befor ) %>%
        dplyr::mutate(previous_hospi = "YES") %>%
        dplyr::select( patient_num, GROUP,previous_hospi) %>%
        dplyr::distinct()%>%
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
        dplyr::distinct()%>%
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

    ## remove ICD-10 code no specific for elixhauser classification:

    icd_remove_code = c("Z45","Z49","Z50","Z71","Z72","Z94","Z95","Z99")

    loc_pat_obs_elix <- LocalPatientObservations%>%
      filter(! concept_code %in%  icd_remove_code)

    comorb_names_elix <- get_quan_elix_names()
    comorbs_elix <- as.vector(comorb_names_elix$Abbreviation)

    comorb_elix_before <- map_char_elix_codes(
        df = loc_pat_obs_elix,
        comorb_names = comorb_names_elix,
        t1 = -365,
        t2 = limit_d_befor,
        map_type = 'elixhauser',
        truncate = TRUE
    )

    comorb_elix_before <- comorb_elix_before$index_scores %>%
        rename('elixhauser_score' = van_walraven_score)%>%
        mutate(patient_num = as.integer(patient_num))

    output_elix_before <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_before, by = "patient_num")%>%
        dplyr::group_by(GROUP) %>%
        dplyr::summarise(elix_mean_min7=mean(elixhauser_score, na.rm = TRUE),
                         elix_std_min7=sd(elixhauser_score, na.rm = TRUE))%>%
        pivot_longer(c(elix_mean_min7,elix_std_min7),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="patient")%>%
        dplyr::mutate(periode_group="all")%>%
        data.frame()

    output_elix_before_p <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_before, by = "patient_num")%>%
        dplyr::group_by(GROUP,periode_group) %>%
        dplyr::summarise(elix_mean_min7=mean(elixhauser_score, na.rm = TRUE),
                         elix_std_min7=sd(elixhauser_score, na.rm = TRUE))%>%
        pivot_longer(c(elix_mean_min7,elix_std_min7),names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="patient")%>%
        data.frame()

    output_gen=rbind(output_gen,output_elix_before,output_elix_before_p)

    comorb_elix_90 <- map_char_elix_codes(
        df = loc_pat_obs_elix,
        comorb_names = comorb_names_elix,
        t1 = -365,
        t2 = 90,
        map_type = 'elixhauser',
        truncate = TRUE
    )

    comorb_elix_90 <- comorb_elix_90$index_scores %>%
        rename('elixhauser_score' = van_walraven_score)%>%
        mutate(patient_num = as.integer(patient_num))%>%
        data.frame()

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

    output_elix_before_score <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_before, by = "patient_num")%>%
        dplyr::group_by(GROUP) %>%
        dplyr::summarise(across(comorbs_elix, ~sum(.x, na.rm = TRUE)))%>%
        pivot_longer(comorbs_elix,names_to = "variable", values_to = "value")%>%
        dplyr::mutate(concept="elix",
                      periode_group="all",
                      time="before")%>%
        data.frame()

    output_elix_before_score_p <- LocalPatientSummary  %>%
        dplyr::inner_join(comorb_elix_before, by = "patient_num")%>%
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


    out_elix=rbind(output_elix_before_score,output_elix_before_score_p,output_elix_90_score,output_elix_90_score_p)

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
        dplyr::left_join( FourCePhase2.1ards::pheno_ICD[,c("ICD","ICD.String","PheCode","Phenotype")],by =c( "concept_code" = "ICD"))%>%
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
        dplyr::left_join( FourCePhase2.1ards::comp_class,by =c( "concept_code" = "ICD10"))%>%
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
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT")) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="all")%>%
        data.frame()

    out_proc_all_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT")) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="all")%>%
        data.frame()

    # before time
    out_proc_before <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT")) %>%
        dplyr::filter( days_since_admission < 0) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="before")%>%
        data.frame()

    out_proc_before_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT"))  %>%
        dplyr::filter( days_since_admission < 0) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="before")%>%
        data.frame()

    # after time
    out_proc_after <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT"))  %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="after")%>%
        data.frame()

    out_proc_after_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT"))%>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code,concept_type, GROUP,periode_group) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type)%>%
        dplyr::mutate(time="after")%>%
        data.frame()

    # after 90 days after admission time
    out_proc_after90 <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT")) %>%
        dplyr::filter( days_since_admission >= 0) %>%
        dplyr::group_by(concept_code, GROUP, concept_type) %>%
        dplyr::summarise(value=n_distinct(patient_num))%>%
        dplyr::rename(variable=concept_code,
                      concept=concept_type )%>%
        dplyr::mutate(periode_group="all",
                      time="after_90")%>%
        data.frame()

    out_proc_after90_p <- LocalPatientObservations %>%
        dplyr::filter( concept_type %in% c("PROC-ICD10","PROC-ICD9","PROC-CPT"))  %>%
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

    patnum_day <- LocalPatientSummary %>%
        dplyr::select("patient_num","GROUP","periode_group")%>%
        dplyr::distinct()

    patnum_day=patnum_day[rep(1:nrow(patnum_day),times = 3),]
    patnum_day$days_since_admission <- rep(c(6,27,89), each=nrow(patnum_day)/3)

    status7 <- LocalPatientClinicalCourse %>%
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","GROUP","periode_group"))%>%
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
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","GROUP","periode_group"))%>%
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
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","GROUP","periode_group"))%>%
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
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","GROUP","periode_group"))%>%
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
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","GROUP","periode_group"))%>%
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
        dplyr::right_join(patnum_day,by = c("patient_num","days_since_admission","GROUP","periode_group"))%>%
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

    message("Analysis => OK")

    ## ========================================
    ## PART 5 : univariate analysis
    ## ========================================

    # univariate analysis ARDS  versus No ARDS young

    yi_vi_uni <- function(ayg,outcome1,outcome0,popu_eval,popu_ref){

      temp_metha <- ayg %>%
        dplyr::filter (GROUP %in% c(outcome0,outcome1)) %>%
        dplyr::filter(variable %in% c( popu_eval,popu_ref)) %>%
        mutate(
          GROUP2 = case_when(
            GROUP==outcome0 ~ "outcome0",
            GROUP==outcome1 ~ "outcome1"
          ))%>%
        mutate(
          variable2 = case_when(
            variable==popu_eval ~ "popu_eval",
            variable==popu_ref ~ "popu_ref",
          ))%>%
        select(GROUP2, variable2, value, periode_group)%>%
        pivot_wider(names_from = c(GROUP2, variable2), values_from = value)%>%
        replace(is.na(.), 0)

      if (ncol(temp_metha) ==5 ){

      temp_yivi <- metafor::escalc(measure = "RR",
                             ai = outcome1_popu_eval,
                             bi = outcome0_popu_eval,
                             ci = outcome1_popu_ref,
                             di = outcome0_popu_ref,
                             data = temp_metha,
                             append = TRUE)

      output_metha=data.frame("periode_group"=temp_yivi$periode_group,
                             "outcome0" = outcome0,
                             "outcome1" = outcome1,
                             "popu_eval" = popu_eval,
                             "popu_ref" = popu_ref,
                             "yi"=temp_yivi$yi,
                             "vi"=temp_yivi$vi)

      } else  {

        output_metha=data.frame("periode_group"=c("all","P2", "P1"),
                                "outcome0" = outcome0,
                                "outcome1" = outcome1,
                                "popu_eval" = popu_eval,
                                "popu_ref" = popu_ref,
                                "yi"=miss_value,
                                "vi"=miss_value)

      }

      return(output_metha)
    }


    uni_age= yi_vi_uni(ayg=output_gen ,
                         outcome1= "ARDS_18_49",
                         outcome0 ="NO_ARDS_18_49",
                         popu_eval= "26to49",
                         popu_ref= "18to25")

    uni_sex= yi_vi_uni(ayg=output_gen ,
                       outcome1= "ARDS_18_49",
                       outcome0 ="NO_ARDS_18_49",
                       popu_eval= "male",
                       popu_ref= "female")

    uni_prehosp= yi_vi_uni(ayg=output_gen ,
                       outcome1= "ARDS_18_49",
                       outcome0 ="NO_ARDS_18_49",
                       popu_eval= 'previous_hospi',
                       popu_ref= 'no_previous_hospi')

    ## comorbidities

    ## number of patient per groups, period, and befor and after

    num_pat_all <- output_gen %>%
      filter(concept == "patient")%>%
      filter(variable == "number")%>%
      select(GROUP, value, periode_group)%>%
      mutate(time = "all")

    num_pat_before <- output_gen %>%
      filter(concept == "patient")%>%
      filter(variable %in% c('previous_hospi'))%>%
      select(GROUP, value, periode_group)%>%
      mutate(time = "before")

    num_pat_gen =rbind(num_pat_all,num_pat_before)

    ## comorbidities comparison between ards - no ards young

    outcome1="ARDS_18_49"
    outcome0="NO_ARDS_18_49"

    variable_uni= c(unique(elix_score$elix_short),unique(substr(ICD10$code,1,3)))

    comb_concept_G=expand.grid(variable_uni,c(outcome0,outcome1),unique(out_comp$periode_group), unique(num_pat_gen$time))
    colnames(comb_concept_G)=c("variable","GROUP","periode_group","time")

    comb_concept_G<-comb_concept_G %>%
      left_join(num_pat_gen,by = c("GROUP","periode_group", "time"))%>%
      rename(npat = value)

    uni_comor_format <- out_proc_diag_med %>%
      filter(periode_group %in% c("all","before")) %>%
      filter(GROUP %in% c(outcome1, outcome0)) %>%
      filter(variable %in% variable_uni)%>%
      filter(time %in% unique(num_pat_gen$time))%>%
      mutate(popu_eval =value)%>%
      right_join(comb_concept_G,by = c("variable","GROUP","periode_group", "time"))%>%
      mutate(popu_eval =ifelse(is.na(popu_eval),0,popu_eval))%>%
      mutate(popu_ref = npat - popu_eval)%>%
      select(GROUP,variable,periode_group,time,popu_eval, popu_ref)%>%
      mutate( GROUP = case_when(
        GROUP==outcome0 ~ "outcome0",
        GROUP==outcome1 ~ "outcome1"))%>%
      pivot_wider(names_from = c(GROUP), values_from = c(popu_eval, popu_ref))%>%
      data_frame()

    uni_comor <- metafor::escalc(measure = "RR",
                      ai = popu_eval_outcome1,
                      bi = popu_eval_outcome0,
                      ci = popu_ref_outcome1,
                      di = popu_ref_outcome0,
                           data = uni_comor_format,
                           append = TRUE)

    uni_comor <- uni_comor%>%
      select("variable","periode_group","time","yi","vi")


    ## univariate analysis Young  versus Old  ARDS
    #  Age

    outcome1="male"
    outcome0="female"
    popu_eval="ARDS_18_49"
    popu_ref="ARDS_sup_49"

    uni_sex_ards <- output_gen %>%
      dplyr::filter ( variable %in% c(outcome0,outcome1)) %>%
      dplyr::filter( GROUP %in% c( popu_eval,popu_ref)) %>%
      mutate(
        GROUP2 = case_when(
          GROUP==popu_ref ~ "popu_ref",
          GROUP==popu_eval ~ "popu_eval"
        ))%>%
      mutate(
        variable2 = case_when(
          variable==outcome1 ~ "outcome1",
          variable==outcome0 ~ "outcome0",
        ))%>%
      select(GROUP2, variable2, value, periode_group)%>%
      pivot_wider(names_from = c(variable2,GROUP2), values_from = value)%>%
      replace(is.na(.), 0)

    if (ncol(uni_sex_ards) ==5 ){

    uni_sex_ards_yivi <- metafor::escalc(measure = "RR",
                        ai = outcome1_popu_eval,
                        bi = outcome0_popu_eval,
                        ci = outcome1_popu_ref,
                        di = outcome0_popu_ref,
                        data = uni_sex_ards,
                        append = TRUE)

    uni_sex_ards_yivi_out=data.frame("periode_group"=uni_sex_ards_yivi$periode_group,
                            "outcome0" = outcome0,
                            "outcome1" = outcome1,
                            "popu_eval" = popu_eval,
                            "popu_ref" = popu_ref,
                            "yi"=uni_sex_ards_yivi$yi,
                            "vi"=uni_sex_ards_yivi$vi)
    } else  {

      uni_sex_ards_yivi_out=data.frame("periode_group"=c("all","P2", "P1"),
                                       "outcome0" = outcome0,
                                       "outcome1" = outcome1,
                                       "popu_eval" = popu_eval,
                                       "popu_ref" = popu_ref,
                                       "yi"=miss_value,
                                       "vi"=miss_value)

    }


    uni_demo = rbind(uni_age,uni_sex,uni_prehosp,uni_sex_ards_yivi_out )

    #  complication

    popu_eval = "ARDS_18_49"
    popu_ref = "ARDS_sup_49"

    var_uni_comp= c(unique(comp_class$complication_class),unique(substr(ICD10$code,1,3)))

    comb_concept_ARDS=expand.grid(var_uni_comp,c(popu_eval,popu_ref),unique(out_comp$periode_group), "all")
    colnames(comb_concept_ARDS)=c("variable","GROUP","periode_group","time")

    comb_concept_ARDS<-comb_concept_ARDS %>%
      left_join(num_pat_gen,by = c("GROUP","periode_group","time"))%>%
      rename(npat = value)

    uni_comp_format <- out_proc_diag_med %>%
      filter(GROUP %in% c(popu_eval, popu_ref)) %>%
      filter(variable %in% var_uni_comp )%>%
      filter(time=="after")%>%
      mutate(outcome1 =value)%>%
      right_join(comb_concept_ARDS,by = c("variable","GROUP","periode_group"))%>%
      mutate(outcome1 =ifelse(is.na(outcome1),0,outcome1))%>%
      mutate(outcome0 = npat - outcome1)%>%
      select(GROUP,variable,periode_group,outcome1, outcome0)%>%
      mutate( GROUP = case_when(
        GROUP==popu_eval ~ "popu_eval",
        GROUP==popu_ref ~ "popu_ref"))%>%
      pivot_wider(names_from = c(GROUP), values_from = c(outcome1, outcome0))%>%
      data_frame()

    uni_comp <- metafor::escalc(measure = "RR",
                           ai = outcome1_popu_eval,
                           bi = outcome0_popu_eval,
                           ci = outcome1_popu_ref,
                           di = outcome0_popu_ref,
                            data = uni_comp_format,
                            append = TRUE)

    uni_comp <- uni_comp%>%
      select("variable","periode_group","yi","vi")


    message("Univariate analysis  => OK")

    ## ========================================
    ## PART 5 : Multivariate analysis
    ## ========================================

    run_logicregression <-
      function(name, currSiteId,df, depend_var, ind_vars) {

        independ_vars <- paste(ind_vars, collapse = ' + ')
        output <- tryCatch(
          {glm(as.formula(paste(depend_var, '~', independ_vars)),
               family = 'binomial', data = df) %>%
              summary()},
          error = function(cond) {
            message(paste("Error when regressing", depend_var))
            message("Original error message:")
            message(cond)
            message('Skipping for now...')
            return(NULL) # return NA in case of error
          }
        )

        if (!is.null(output)){

        multi=data.frame(output$coefficients)
        multi$siteid=currSiteId
        multi$variable <- rownames(multi)
        multi$name=name
        }
        else { multi = NULL}

      return(multi)
      }

    data_multi <- LocalPatientSummary%>%
      mutate(pre_hospit = ifelse((patient_num %in% pre_hospit_tmp$patient_num),1,0))%>%
      inner_join(comorb_elix_90, by = c("patient_num"))%>%
      filter(GROUP %in% c("ARDS_18_49","NO_ARDS_18_49"))

    LR_ALL <-run_logicregression(
      name = "ALL",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Pulmonary","Renal","DM","Liver","Obesity","Alcohol","Drugs","CHF","HTN"))

    LR_ALL_wout_RENAL <-run_logicregression(
      name = "ALL_wout_RENAL",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Pulmonary","DM","Liver","Obesity","Alcohol","Drugs","CHF","HTN"))

    LR_ALL_wout_Abuses <-run_logicregression(
      name = "ALL_wout_Abuses",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Pulmonary","DM","Liver","Obesity","CHF","HTN"))

    LR_CHF <-run_logicregression(
      name = "CHF",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Obesity","DM","HTN","Alcohol","CHF"))

    LR_OBESITY <-run_logicregression(
      name = "OBESITY",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Obesity","DM","Alcohol"))

    LR_DM <-run_logicregression(
      name = "DM",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Obesity","DM"))

    LR_HT <-run_logicregression(
      name = "HT",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Obesity","DM","Renal","HTN"))

    LR_LIVER <-run_logicregression(
      name = "LIVER",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Obesity","DM","HTN","Liver"))

    LR_RENAL <-run_logicregression(
      name = "RENAL",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Obesity","DM","HTN","Renal"))


    LR_ABUSES <-run_logicregression(
      name = "ABUSES",
      currSiteId = currSiteId,
      df= data_multi,
      depend_var= "ARDS",
      ind_vars= c("sex","Alcohol","Drugs"))

    multiresults <- rbind(LR_ALL,LR_ALL_wout_RENAL,LR_ALL_wout_Abuses,LR_CHF,LR_OBESITY,LR_DM, LR_HT, LR_LIVER, LR_RENAL,LR_ABUSES)

    message("Multivariate analysis  => OK")
    ## ========================================
    ## PART 6 : Saving output
    ## ========================================

    siteid=unique(LocalPatientSummary$siteid)

    output_day$siteid=siteid
    output_gen$siteid=siteid
    output_lab$siteid=siteid
    out_proc_diag_med$siteid=siteid
    output_sens$siteid=siteid
    uni_demo$siteid=siteid
    uni_comor$siteid=siteid
    uni_comp$siteid=siteid

    ## manage std with 1 patient
    output_lab <- output_lab %>%
        dplyr::mutate(std_value = ifelse(npat == 1, 0, std_value),
                      std_log_value = ifelse(npat == 1, 0, std_log_value))

    var_gen_num = c("number","n_severe","npat_nodiag","never_severe","ever_severe","26to49","50to69","70to79","80plus","18to25","female","male","other","no_previous_hospi","previous_hospi","no_rehospit","rehospit","dead","dead_less_28","dead_less_90","out_hospit_28","out_hospit_90","in_hospital", "death", "out_hospital","no_data")

    var_gen_ag =c("mean_freq","std_freq", "mean_los", "std_los")

    ## manage obfuscation

    ## add obfuscation

    obfusc=data.frame("siteid" = siteid, "truefalse" = obfuscation, "value" = obfuscationThreshord)

    if(obfuscation){

        group_less_obf <- output_gen %>%
            dplyr::filter(variable =="number" &  value < obfuscationThreshord)%>%
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
    write.csv(multiresults, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_multiresults",".csv")), row.names = FALSE, na = "")
    write.csv(uni_demo, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_uni_demo",".csv")), row.names = FALSE, na = "")
    write.csv(uni_comor, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_uni_comor",".csv")), row.names = FALSE, na = "")
    write.csv(uni_comp, file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_uni_comp",".csv")), row.names = FALSE, na = "")

    ## save multi variate analysis

    cat(
      "Result is saved in",
      file.path(
        getProjectOutputDirectory()
      ),
      "\nPlease submit the result file by running submitAnalysis()\n"
    )



}

