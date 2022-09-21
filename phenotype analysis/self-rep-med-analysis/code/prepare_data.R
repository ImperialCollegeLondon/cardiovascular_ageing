library(tidyverse)
################################################################################
source("code/functions.R")

################################################################################
setwd("data/")
################################################################################
df0 <- read_csv("./raw_df.csv.gz")
################################################################################

#df0 %>%
#  dplyr::rename(alcohol=p1558_i2,
#                smoking=p20116_i2
#                ) %>%
#  mutate(alcohol=if_else(alcohol==-3, NA_real_, alcohol)) %>%
#  mutate(alcohol=as.factor(alcohol)) %>%
#  mutate(smoking=if_else(smoking==-3, NA_real_, smoking)) %>%
#  mutate(smoking=as.factor(smoking)) ->
#  df0
################################################################################
prescp.df <- read_csv("./annotated_primary_care_prescriptions_cmri_34k_cohort.csv.gz") %>%
  mutate(SID=as.integer(SID))
################################################################################
prescp.df %>%
  filter(str_detect(target, fixed("HMG-CoA reductase inhibitor",ignore_case=TRUE)))-> # Statins
  statin.prescp.df
################################################################################
prescp.df %>%
  filter(name=="Metformin")->
  metformin.prescp.df
################################################################################
prescp.df %>%
  filter(name=="Digoxin")->
  digoxin.prescp.df
################################################################################
prescp.df %>%
  filter(str_starts(atc, "C07"))->
  BB.prescp.df
################################################################################
prescp.df %>%
  filter(str_starts(atc, "C09A|C09B"))->
  ACEi.prescp.df
################################################################################
prescp.df %>%
  filter(target=="Type-1 angiotensin II receptor antagonist")->
  ARB.prescp.df
################################################################################
prescp.df %>%
  filter(str_starts(atc, "C08"))->
  CCB.prescp.df
################################################################################
prescp.df %>%
  filter(str_starts(atc, "C03"))->
  diuretics.prescp.df
################################################################################
# Get self-reported medication information

#This is self-reported medication data as determined at the imaging visit!
#  https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=20003

#See also https://biobank.ctsu.ox.ac.uk/crystal/instance.cgi?id=2

# "All participants attended an initial assessment centre.
# A proportion were invited several years later to repeat the assessment.
# For some participants the repeat visit also captured information which was not gathered during their initial visit.
df0 %>%
  select(SID, p20003_i2) %>%
  drop_na(p20003_i2) %>%
  mutate(p20003_i2=stringr::str_replace(p20003_i2, fixed("["), "")) %>%
  mutate(p20003_i2=stringr::str_replace(p20003_i2, fixed("]"), "")) %>%
  tidyr::separate_rows(p20003_i2,sep=",") %>%
  transmute(SID, coding=as.integer(p20003_i2)) %>%
  inner_join(read_tsv("./coding4.tsv"), by="coding") ->
  self.reported.med.df
################################################################################
(statin.prescp.df %>%
   pull(name) %>%
   unique %>%
   stringr::str_to_lower() -> statin_drug_names)

statin_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(SID) %>%
  summarise() ->
  statin.self.reported.df

#statin.self.reported.df %>%
#  write_csv("./Derived/SIDs_statin.self.reported.csv.gz")


get_df.bin.self.reported(statin.self.reported.df,df0) %>%
  mutate(Drug="Statins") ->
  dataframes_drugs.df

statin_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(meaning) %>%
  summarise() %>%
  mutate(Drug="Statins") ->
  drug_names_self_reported.df

  #write_csv("./Derived/dataframe_statin_self_reported.csv.gz")
################################################################################
(metformin.prescp.df %>%
   pull(name) %>%
   unique %>%
   stringr::str_to_lower() -> metformin_drug_names)

metformin_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(SID) %>%
  summarise() ->
  metformin.self.reported.df

#metformin.self.reported.df %>%
#  write_csv("./Derived/SIDs_metformin.self.reported.csv.gz")


get_df.bin.self.reported(metformin.self.reported.df,df0) %>%
  mutate(Drug="Metformin") %>%
  bind_rows(dataframes_drugs.df) ->
  dataframes_drugs.df


metformin_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(meaning) %>%
  summarise() %>%
  mutate(Drug="Metformin") %>%
  bind_rows(drug_names_self_reported.df)->
  drug_names_self_reported.df

################################################################################
(digoxin.prescp.df %>%
   pull(name) %>%
   unique %>%
   stringr::str_to_lower() -> digoxin_drug_names)

digoxin_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(SID) %>%
  summarise() ->
  digoxin.self.reported.df

#metformin.self.reported.df %>%
#  write_csv("./Derived/SIDs_metformin.self.reported.csv.gz")


get_df.bin.self.reported(digoxin.self.reported.df,df0) %>%
  mutate(Drug="Digoxin") %>%
  bind_rows(dataframes_drugs.df) ->
  dataframes_drugs.df
  #write_csv("./Derived/dataframe_metformin_self_reported.csv.gz")


digoxin_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(meaning) %>%
  summarise() %>%
  mutate(Drug="Digoxin") %>%
  bind_rows(drug_names_self_reported.df)->
  drug_names_self_reported.df
################################################################################
(BB.prescp.df %>%
   pull(name) %>%
   unique %>%
   stringr::str_to_lower() -> BB_drug_names)

BB_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(SID) %>%
  summarise() ->
  BB.self.reported.df

#BB.self.reported.df %>%
#  write_csv("./Derived/SIDs_BB.self.reported.csv.gz")


get_df.bin.self.reported(BB.self.reported.df,df0) %>%
  mutate(Drug="BB") %>%
  bind_rows(dataframes_drugs.df) ->
  dataframes_drugs.df

  #write_csv("./Derived/dataframe_BB_self_reported.csv.gz")


BB_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(meaning) %>%
  summarise() %>%
  mutate(Drug="BB") %>%
  bind_rows(drug_names_self_reported.df)->
  drug_names_self_reported.df
################################################################################
(ACEi.prescp.df %>%
   pull(name) %>%
   unique %>%
   stringr::str_to_lower() -> ACEi_drug_names)

ACEi_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(SID) %>%
  summarise() ->
  ACEi.self.reported.df

#ACEi.self.reported.df %>%
#  write_csv("./Derived/SIDs_ACEi.self.reported.csv.gz")


get_df.bin.self.reported(ACEi.self.reported.df,df0) %>%
  mutate(Drug="ACEi") %>%
  bind_rows(dataframes_drugs.df) ->
  dataframes_drugs.df


ACEi_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(meaning) %>%
  summarise() %>%
  mutate(Drug="ACEi") %>%
  bind_rows(drug_names_self_reported.df)->
  drug_names_self_reported.df

  #write_csv("./Derived/dataframe_ACEi_self_reported.csv.gz")
################################################################################
(ARB.prescp.df %>%
   pull(name) %>%
   unique %>%
   stringr::str_to_lower() -> ARB_drug_names)

ARB_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(SID) %>%
  summarise() ->
  ARB.self.reported.df

#ARB.self.reported.df %>%
#  write_csv("./Derived/SIDs_ARB.self.reported.csv.gz")

get_df.bin.self.reported(ARB.self.reported.df,df0) %>%
  mutate(Drug="ARB") %>%
  bind_rows(dataframes_drugs.df) ->
  dataframes_drugs.df


ARB_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(meaning) %>%
  summarise() %>%
  mutate(Drug="ARB") %>%
  bind_rows(drug_names_self_reported.df)->
  drug_names_self_reported.df

  #write_csv("./Derived/dataframe_ARB_self_reported.csv.gz")
################################################################################
(CCB.prescp.df %>%
   pull(name) %>%
   unique %>%
   stringr::str_to_lower() -> CCB_drug_names)

CCB_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(SID) %>%
  summarise() ->
  CCB.self.reported.df

#CCB.self.reported.df %>%
#  write_csv("./Derived/SIDs_CCB.self.reported.csv.gz")

get_df.bin.self.reported(CCB.self.reported.df,df0) %>%
  mutate(Drug="CCB") %>%
  bind_rows(dataframes_drugs.df) ->
  dataframes_drugs.df



CCB_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(meaning) %>%
  summarise() %>%
  mutate(Drug="CCB") %>%
  bind_rows(drug_names_self_reported.df)->
  drug_names_self_reported.df

  #write_csv("./Derived/dataframe_CCB_self_reported.csv.gz")
################################################################################
(diuretics.prescp.df %>%
   pull(name) %>%
   unique %>%
   stringr::str_to_lower() -> diuretics_drug_names)

diuretics_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(SID) %>%
  summarise() ->
  diuretics.self.reported.df

#diuretics.self.reported.df %>%
#  write_csv("./Derived/SIDs_diuretics.self.reported.csv.gz")

get_df.bin.self.reported(diuretics.self.reported.df,df0) %>%
  mutate(Drug="Diuretics") %>%
  bind_rows(dataframes_drugs.df) ->
  dataframes_drugs.df



diuretics_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(meaning) %>%
  summarise() %>%
  mutate(Drug="Diuretics") %>%
  bind_rows(drug_names_self_reported.df)->
  drug_names_self_reported.df
  #write_csv("./Derived/dataframe_diuretics_self_reported.csv.gz")
################################################################################
################################################################################
dataframes_drugs.df %>%
  write_csv("./Derived/dataframes_drugs_self_reported.csv.gz")
################################################################################
drug_names_self_reported.df %>%
  write_csv("./Derived/self_reported_drugs_used.csv")
