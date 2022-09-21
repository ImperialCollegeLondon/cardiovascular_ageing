library(tidyverse)
library(rstanarm)
################################################################################
################################################################################
setwd("~/github/ecg-age-prediction-retrain/self-reported-medication-analyses-MRI/data/")
################################################################################
df0 <- read_csv("./raw_df.csv.gz")
################################################################################
################################################################################
prescp.df <- read_csv("./annotated_primary_care_prescriptions_cmri_34k_cohort.csv.gz") %>%
  mutate(SID=as.integer(SID))
################################################################################
prescp.df %>%
  filter(str_starts(atc, "C07"))->
  BB.prescp.df
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
################################################################################
(BB.prescp.df %>%
   pull(name) %>%
   unique %>%
   stringr::str_to_lower() -> BB_drug_names)

BB_drug_names %>%
  map_dfr(~{ self.reported.med.df %>%
      filter(str_detect(meaning, fixed(.x, ignore_case=TRUE)))}) %>%
  group_by(meaning,SID) %>%
  summarise() %>%
  ungroup()->
  BB.self.reported.df
################################################################################
# only keep drugs or combinations with more than one subjects
BB.self.reported.df %>%
  group_by(meaning) %>%
  summarise(Nsubj=length(unique(SID))) %>%
  filter(Nsubj > 1) %>%
  select(-Nsubj) %>%
  left_join(BB.self.reported.df, by="meaning") ->
  BB.self.reported.df
################################################################################
df %>%
  select(SID) %>%
  left_join(p.df %>% transmute(SID, Drug=meaning, PriorMed=1),by="SID") %>%
  pivot_wider(id_cols=SID, names_from=Drug, values_from=PriorMed, values_fill =0) %>%
  select(-`NA`)->
  tmp

tmp %>%
  pivot_longer(-SID, names_to="Drug", values_to="PriorMed")->
  tmp

# How to treat individuals that have more than one self-reported?
# pick one randomly

set.seed(42)
tmp%>%
  filter(PriorMed==1) %>%
  group_by(SID) %>%
  summarise(N=n()) %>%
  #filter(N>1) %>%
  left_join(tmp %>% filter(PriorMed==1), by="SID") %>%
  sample_frac(1) %>% # shuffle
  group_by(SID) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  select(-N)->
  tmp

tmp %>%
  select(-PriorMed) %>%
  right_join(df %>%
               select(SID,
                     age_delta,
                     Sex,
                     Age,
                     packyrs,
                     alcoholintakegpd,
                     #z_adj_2, # telomere
                     Hypertension,
                     Obesity,
                     `Diabetes Mellitus`,
                     CAD,
                     Hyperchol,
                     `Heart Failure`,
                     bmi,
                     sbp,
                     dbp,
                     pulse_rate
                  ) %>%
               drop_na()
               ,
             by="SID") ->
  df.fit

df.fit %>%
  mutate(Drug=if_else(is.na(Drug), "No BB", Drug))->
  df.fit
################################################################################
M1_stanlmer <- stan_lmer(formula = course ~ 1 + (1 | school),
                         data = GCSE,
                         seed = 349)

post <-
  stan_lmer(
    age_delta ~ 1+ (1|Drug) + Age + Sex +  bmi + sbp + dbp + pulse_rate + packyrs + alcoholintakegpd +Hypertension + Obesity + `Diabetes Mellitus` +
      CAD+Hyperchol + `Heart Failure`,
    QR=TRUE,
    cores=2,
    data = df.fit,
    seed = 12345
  )
post

summary(post)
post %>%
  bayesplot::mcmc_intervals(pars=vars(-sigma, -`(Intercept)`))

post %>%
  bayesplot::mcmc_intervals(pars=vars(starts_with("b[(Intercept)")), prob_outer=.95)


post %>%
  bayesplot::mcmc_areas_ridges(pars=vars(starts_with("b[(Intercept)")))

################################################################################

BB.df <- read_csv("~/github/ecg-age-prediction-retrain/self-reported-medication-analyses-MRI/data/Derived/dataframes_drugs_self_reported.csv.gz") %>%
  filter(Drug=="BB")

post2 <-
  stan_lm(
    age_delta ~ PriorMed + Age + Sex +  bmi + sbp + dbp + pulse_rate + packyrs + alcoholintakegpd +Hypertension + Obesity + `Diabetes Mellitus` +
      CAD+Hyperchol + `Heart Failure`,
    cores=2,
    data = BB.df,
    prior = R2(location = 0.1),
    seed = 12345
  )
post2

post2 %>%
  bayesplot::mcmc_intervals(pars=vars(-sigma, -`(Intercept)`), prob_outer = .95)



post2 %>%
  bayesplot::mcmc_pairs(pars=vars(-sigma, -`(Intercept)`, -R2, -bmi, -sb))
