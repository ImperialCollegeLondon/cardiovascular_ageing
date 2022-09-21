library(UKBRlib)
library(tidyverse)
################################################################################
setwd("data/")
################################################################################
dxproj <- "project-G7b1FVQJ4fX8FYYz4QyPPk4p"
################################################################################
cdnx_mri <- DX$new(project_id=dxproj)
cdnx_mri$download_file(file="file-GB6Qgf0J4fX9QyZvJbG8z0j4")

df <- read_csv("./eid_andcardiacAGE_May2022_primaryphenotype.csv")
df <- df %>%
  dplyr::rename(SID=eid_40616,
                age_delta=catb_delta_with_t1_bc_cole) %>%
  dplyr::select(SID, age_delta) %>%
  dplyr::arrange(age_delta) %>%
  dplyr::group_by(SID) %>%
  dplyr::slice_head(n=1) %>%
  dplyr::ungroup() %>%
  mutate(SID=as.integer(SID))


################################################################################

robj <- UkbRap$new(project_id = dxproj)
sa <- robj$fetch_sex_age(instance = 2, GP_provider = FALSE)

sa %>%
  filter(Instance==2) %>%
  select(SID, Sex, Age) %>%
  drop_na(Sex, Age) %>%
  mutate(SID=as.integer(SID)) %>%
  inner_join(df, by="SID")->
  df


################################################################################
fields.tmp <- robj$fetch_ukb_fields(c(1558, 20116, 6177))
fields.tmp %>%
  as_tibble() %>%
  mutate(SID=as.integer(SID)) %>%
  right_join(df, by="SID")->
  df


################################################################################
quantsR6 <- DxPhenotypesQuant$new(project_id=dxproj)
quant_data <- quantsR6$fetch_data()
quant_data2 <- quant_data %>% dplyr::filter(visit==2)  %>% dplyr::select(SID, bmi, sbp, dbp,  map_adj, pulse_rate, ldl_direct)
quant_data2 %>%
  mutate(SID=as.integer(SID)) %>%
  right_join(df, by="SID")->
  df

################################################################################
## Get visit dates
robj$fetch_visit_time() ->
  df.visit.date

df.visit.date %>%
  filter(visitn==2) ->
  df.visit.date

df.visit.date %>%
  mutate(SID=as.integer(SID)) %>%
  right_join(df, by="SID")->
  df

################################################################################
# self-reported medication

robj$fetch_ukb_fields(c(20003)) %>%
  select(SID, p20003_i2) %>%
  mutate(SID=as.integer(SID)) %>%
  right_join(df, by="SID")->
  df

################################################################################
# get CV risk factors as determined by Mit Shah
cdnx_mri$list_files("/Derived/Phenotypes/cardiac_age_cmri/Groups_MS/") -> files.df

names <- files.df$name
ids<-files.df$id

nfiles <- length(names)
dfs <- list()
for (i in 1:nfiles) {
  cdnx_mri$download_file(file=ids[[i]])
  dfs[[names[[i]]]] <- data.table::fread(names[[i]]) %>% as_tibble()
}

read_csv("./HTN_unselected11049.csv") %>% pull(eid_40616) %>%  as.integer() -> htn_ids
read_csv("./Obese_7089_unselected.csv") %>% pull(eid_40616)  %>%  as.integer()-> obesity_ids
read_csv("./diabetes_2474_unselected.csv") %>%  pull(eid_40616)  %>%  as.integer()-> t2d_ids
read_csv("./cad_2657_unselected.csv") %>% pull(eid_40616)  %>%  as.integer()-> cad_ids
read_csv("./hyperchol_7484_unselected.csv") %>% pull(eid_40616)  %>%  as.integer()-> hyperchol_ids
read_csv("./hf_502_unselected.csv") %>% pull(eid_40616)  %>%  as.integer() -> hf_ids


df %>%
  mutate(Hypertension=SID %in% htn_ids) %>%
  mutate(Obesity=SID %in% obesity_ids) %>%
  mutate(`Diabetes Mellitus` = SID %in% t2d_ids) %>%
  mutate(CAD = SID%in% cad_ids) %>%
  mutate(Hyperchol=SID%in% hyperchol_ids) %>%
  mutate(`Heart Failure`=SID %in% hf_ids) ->
  df
################################################################################
# eid_alcohol - alcoholintakegpd is grams per day consumed.
# smoke_exercise - packyrs_instance2  pack year of smoking
# 500k_telomere - z_adj_2 used in analysis

for(fid in c("file-GFyQpF8J4fXP8k1j4j2yXkg2","file-GFyQpF8J4fX1vXYf4jxk8jf3","file-GFyQpF8J4fXP8y5K4xf19958") ) {
  cdnx_mri$download_file(file=fid)
}
smoking.df <- read_csv("smoke_exercise.csv")
smoking.df %>%
  transmute(SID=as.integer(eid_40616), packyrs_instance2) ->
  smoking.df

alcohol.df <- read_csv("eid_alcohol.csv")
alcohol.df %>%
  transmute(SID=as.integer(eid_40616), alcoholintakegpd) ->
  alcohol.df

telomere.df <- read_csv("500k_telomere.csv") %>%
  transmute(SID=as.integer(eid_40616),
            z_adj_2
  )

df %>%
  left_join(smoking.df, by="SID") %>%
  mutate(packyrs=if_else(is.na(packyrs_instance2), 0, packyrs_instance2)) %>%
  select(-packyrs_instance2) %>%
  left_join(alcohol.df, by="SID") %>%
  mutate(alcoholintakegpd=if_else(is.na(alcoholintakegpd), 0, alcoholintakegpd)) %>%
  left_join(telomere.df, by="SID") ->
  df

################################################################################

df %>%
  write_csv("./raw_df.csv.gz")
