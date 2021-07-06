library(icd)
library(tidyverse)

# The below was copied from the Comorbidity_Mapping.Rmd file in order to support sourcing it for other analyses.
# df = obs_raw # (4ce LocalPatientObservations.csv)
# comorb_names <- get_charlson_names() or get_elix_names()
# t1 <- earliest time point to consider comorbidities
# t2 <- latest time point to consider comorbidities
# t1 = -365; t2 = -1;
# example above will map all codes up to a year prior but before admission (admission = day 0)
# num_days_prior_admission = -365 indicates that we consider all codes up to a year prior to the first COVID admission as comorbidities
# day_of
# map_type = 'charlson', 'elixhauser' - where charlson will be scored with quan-deyo
# truncate = TRUE # indicates we are using ICD codes truncated to the first 3 characters; set FALSE if you have full ICD codes

map_char_elix_codes <- function(df, comorb_names, t1, t2, map_type, truncate = TRUE) {

  df <- df %>%
    filter(concept_type %in% c("DIAG-ICD10", "DIAG-ICD9"),
           days_since_admission >= t1 & days_since_admission <= t2)

  # Create separate df frames for ICD9 and 10 Codes
  # icd package does not support simultaneous processing of both ICD code types
  # we will recombine after the initial processing
  icd10 <- df %>%
    select(-c(days_since_admission, value)) %>%
    filter(concept_type == "DIAG-ICD10") %>%
    distinct()

  icd9 <- df %>%
    select(-c(days_since_admission, value)) %>%
    filter(concept_type == "DIAG-ICD9") %>%
    distinct()

  # quan prefixes revised charlson & elixhauser mapping
  if (map_type == "charlson") {
    icd10_comorb_map = icd::icd10_map_quan_deyo
    icd9_comorb_map = icd::icd9_map_quan_deyo
  }
  if (map_type == "elixhauser") {
    icd10_comorb_map = icd::icd10_map_quan_elix
    icd9_comorb_map = icd::icd9_map_quan_elix


    ## add regroupment cancer = Lymphoma+ Tumor +Mets

    icd10_comorb_map$Cancer= c(icd10_comorb_map$Lymphoma,  icd10_comorb_map$Tumor, icd10_comorb_map$Mets )

  }

  ## Because the 4CE has truncated ICD codes, we will also truncate the icd package index maps
  # Function to select first 3 characters of the ICD Code in all lists of the index map
  if (truncate == TRUE) {
    icd10_comorb_map <- lapply(icd10_comorb_map, first_3)
    icd9_comorb_map <- lapply(icd9_comorb_map, first_3)

    ### remove specific ICD10 code if trunc

    icd10_comorb_map$Alcohol= c("F10","K70","T51") # remove"E52" ,"G62", ,"G62" "I42", "Z50","Z71","Z72"
    icd10_comorb_map$Arrhythmia = c("I44","I45","I47","I48","I49","R00","T82","Z45") # remove "Z95"
    icd10_comorb_map$Drugs = c("F11","F12","F13","F14","F15","F16","F18","F19") # remove  "Z71" "Z72"
    icd10_comorb_map$Liver = c("B18","I85","K70","K71","K72","K73","K74","K76") # remove   "Z94" ,"I86","I98"
    icd10_comorb_map$PVD = c("I70","I71","I73","I77","I79","K55" ) # remove   "Z95"
    icd10_comorb_map$Renal = c("I12","I13","N18","N19","N25","Z49") # remove   "Z95","Z94","Z99"
    icd10_comorb_map$HTN = c("I10","I11","I12","I13","I15") # ADD code for HTN complicated "I11","I12","I13","I15"

  }


  if (truncate == FALSE) {
    # convert icd code to short format (without decimals) to facilitate mapping
    # where diagnosis code is you non-truncated icd column
    icd10 <- icd10 %>%
      mutate(diagnosis_code = icd::decimal_to_short(diagnosis_code)) %>%
      select(-concept_code) %>%
      rename(concept_code = diagnosis_code)

    icd9 <- icd9 %>%
      mutate(diagnosis_code = icd::decimal_to_short(diagnosis_code)) %>%
      select(-concept_code) %>%
      rename(concept_code = diagnosis_code)

  }

  # perform the mapping
  icd10_map <-
    icd10_comorbid(
      icd10,
      map = icd10_comorb_map,
      return_df = TRUE,
      visit_name = "patient_num",
      return_binary = TRUE,
    )

  icd9_map <-
    icd9_comorbid(
      icd9,
      map = icd9_comorb_map,
      icd_name = "concept_code",
      return_df = TRUE,
      visit_name = "patient_num",
      return_binary = TRUE
    )

  # If multiple rows due to a patient having both ICD 9 and 10 codes, we will take the max of the column
  # This will allow us to capture the 1s indicating that the comorbidity is present
  # the try wrapper is important in cases there is not an instance of a specific comorbidity in the data - try silences errors
  icd_map <- rbind(icd9_map, icd10_map) %>% # results of 9 and 10 mapping
    group_by(patient_num) %>%
    summarise(
      across(everything(), ~ try(max(.x), silent = TRUE)),
      .groups = 'drop'
    )

  ## Calculate Index Scores
  if (map_type == 'charlson') {
    charlson_score <- icd::charlson_from_comorbid(
      icd_map,
      visit_name = "patient_num",
      scoring_system = "charlson",
      hierarchy = TRUE
    ) %>%
      data.frame(charlson_score = .) %>%
      tibble::rownames_to_column("patient_num")
  }

  # need to check If I've done this correctly - it seems to be that their are different versions
  # of the elixhuaser mapping and that in one HTN is combined. I think this is the version needed
  # to run the van_walraven_from_comorb function
  if (map_type == "elixhauser") {
    # combine hypertension into one category
    icd_map <- icd_map %>%
      mutate(HTN = pmax(HTN, HTNcx, na.rm = TRUE)) %>%
      select(-HTNcx)

    icd_map_wout_K <- icd_map%>%
      select(-Cancer)

    van_walraven_score <- icd::van_walraven_from_comorbid(
      icd_map_wout_K,
      visit_name = 'patient_num',
      hierarchy = TRUE
    ) %>%
      data.frame(van_walraven_score = .) %>%
      tibble::rownames_to_column("patient_num")
  }

  if(map_type == 'charlson') {
    index_scores <- icd_map %>%
      full_join(charlson_score, by = "patient_num") %>%
      arrange(desc(charlson_score)) %>%
      select(
        patient_num,
        charlson_score,
        everything()
      )
  } else {
    index_scores <- icd_map %>%
      full_join(van_walraven_score, by = "patient_num") %>%
      arrange(desc(van_walraven_score)) %>%
      select(
        patient_num,
        van_walraven_score,
        everything()
      )
  }


  # Identify the specific codes that mapped
  # unlist the charlson mapping lists
  icd10$concept_code <- as.character(icd10$concept_code)
  icd9$concept_code <- as.character(icd9$concept_code)

  icd10_map <-
    purrr::map_df(icd10_comorb_map, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Abbreviation", "concept_code")) %>%
    mutate(concept_code = as.character(concept_code)) %>%
    distinct() %>%
    # merge the mapping dataframe to the patient level ICD codes
    # this will return all comorbidities that mapped to our patient data
    inner_join(icd10, by = "concept_code")

  icd9_map <-
    purrr::map_df(icd9_comorb_map, ~ as.data.frame(.x), .id = "name") %>%
    `colnames<-`(c("Abbreviation", "concept_code")) %>%
    mutate(concept_code = as.character(concept_code)) %>%
    distinct() %>%
    inner_join(icd9, by = "concept_code")

  # explain_codes will add additional information regarding the code name
  # add if statements in order to handle sites that only have ICD 9 or 10 codes but not both
  if (nrow(icd10_map) > 0) {
    icd10_mapped_table <- unique(icd10_map$concept_code) %>%
      icd::explain_table() %>%
      left_join(icd10_map, ., by = c("concept_code" = "code")) %>%
      select(patient_num, concept_code, Abbreviation, long_desc) %>%
      distinct()
  }

  if (nrow(icd9_map) > 0) {
    icd9_mapped_table <- unique(icd9_map$concept_code) %>%
      icd::explain_table() %>%
      left_join(icd9_map, ., by = c("concept_code" = "code")) %>%
      select(patient_num, concept_code, Abbreviation, long_desc) %>%
      distinct()
  }

  # Bind both ICD 9 and 10 code tables together
  if(nrow(icd9_map) > 0) {
  mapped_codes_table <-
    icd10_mapped_table %>%
    bind_rows(icd9_mapped_table) %>%
    # calculate how many patients had each unique comorbidity/concept_code
    count(concept_code, Abbreviation, long_desc, name = 'n_patients') %>%
    arrange(desc(n_patients))
  } else {
    mapped_codes_table <-
      icd10_mapped_table %>%
      count(concept_code, Abbreviation, long_desc, name = 'n_patients') %>%
      arrange(desc(n_patients))
  }

  map_results <-
    list(
      index_scores = index_scores,
      mapped_codes_table = mapped_codes_table
    )

  map_results

}

# where df takes in the matrix for the initial mapping

get_table1 <- function(
  df, num_pats, comorbidities,
  pat_col = 'patients', ...
)
  {
  npat_col = sym(paste('n', pat_col, sep = '_'))
  proppat_col = sym(paste('prop', pat_col, sep = '_'))
  comorbidities_map = comorbidities$Abbreviation

  df %>%
    select(all_of(comorbidities_map)) %>%
    colSums() %>%
    data.frame(n_patients = .) %>%
    rownames_to_column("Abbreviation") %>%
    blur_it('n_patients', ...) %>%
    mutate(prop_patients = n_patients/num_pats) %>%
    rename(!!npat_col := n_patients,
           !!proppat_col := prop_patients) %>%
    right_join(comorbidities, ., by = "Abbreviation")
}


get_charlson_names <- function(){
  data.frame(
    Comorbidity = do.call(rbind, icd::names_charlson),
    Abbreviation = do.call(rbind, icd::names_charlson_abbrev))
}

get_quan_elix_names <- function(){
  data.frame(
    Comorbidity = do.call(rbind, icd::names_quan_elix),
    Abbreviation = do.call(rbind, icd::names_quan_elix_abbrev))
}

first_3 <- function(x) {
  # retain first 3 characters of the ICD code
  substr(x, 1, 3) %>% unique()
}
