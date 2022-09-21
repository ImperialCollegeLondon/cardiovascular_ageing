library(tidyverse)
library(glmnet)
################################################################################
source("code/functions.R")
################################################################################
setwd("./data/")
################################################################################
dataframes_drugs.df <- read_csv("./Derived/dataframes_drugs_self_reported.csv.gz")
################################################################################
statin.fit.df <- dataframes_drugs.df %>%
  filter(Drug=="Statins") %>%
  select(-Drug)

get_effect.bin.self.reported(statin.fit.df,add_markers = TRUE) %>%
mutate(AdditionalMarkers=TRUE) %>%
bind_rows(get_effect.bin.self.reported(statin.fit.df, add_markers = FALSE) %>%
              mutate(AdditionalMarkers=FALSE)) -> statin.fits.df

statin.fits.df %>%
  mutate(Drug="Statins")->
  fits.df

  #write_csv("./Derived/statin.fits.csv")

################################################################################
fit <- glm(age_delta ~. ,data=statin.fit.df, family = gaussian)
model.matrix(fit) -> X
X[,colnames(X)[colnames(X)!="(Intercept)"]]  -> X
N <- dim(X)[1]
p <- dim(X)[2]

X <- scale(X)
Y <- statin.fit.df$age_delta %>% scale
fit.glmnet.cv <- cv.glmnet(x=X, y=Y, family="gaussian",
                              intercept = TRUE, alpha=1)

fits.glmnet <- list("Statins"=fit.glmnet.cv)

#saveRDS(fit.glmnet.cv,"./Derived/statin.fit.glm.Rds")

################################################################################
metformin.fit.df <- dataframes_drugs.df %>%
  filter(Drug=="Metformin") %>%
  select(-Drug)

get_effect.bin.self.reported(metformin.fit.df,add_markers = TRUE) %>%
  mutate(AdditionalMarkers=TRUE) %>%
  bind_rows(get_effect.bin.self.reported(metformin.fit.df, add_markers = FALSE) %>%
              mutate(AdditionalMarkers=FALSE)) -> metformin.fits.df

metformin.fits.df %>%
  mutate(Drug="Metformin") %>%
  bind_rows(fits.df)->
  fits.df


  #write_csv("./Derived/metformin.fits.csv")

################################################################################
fit <- glm(age_delta ~. ,data=metformin.fit.df, family = gaussian)
model.matrix(fit) -> X
X[,colnames(X)[colnames(X)!="(Intercept)"]]  -> X
N <- dim(X)[1]
p <- dim(X)[2]

X <- scale(X)
Y <- metformin.fit.df$age_delta %>% scale
fit.glmnet.cv <- cv.glmnet(x=X, y=Y, family="gaussian",
                           intercept = TRUE, alpha=1)
fits.glmnet[["Metformin"]] <- fit.glmnet.cv

#saveRDS(fit.glmnet.cv,"./Derived/metformin.fit.glm.Rds")
################################################################################
digoxin.fit.df <- dataframes_drugs.df %>%
  filter(Drug=="Digoxin") %>%
  select(-Drug)

get_effect.bin.self.reported(digoxin.fit.df,add_markers = TRUE) %>%
  mutate(AdditionalMarkers=TRUE) %>%
  bind_rows(get_effect.bin.self.reported(digoxin.fit.df, add_markers = FALSE) %>%
              mutate(AdditionalMarkers=FALSE)) -> digoxin.fits.df

digoxin.fits.df %>%
  mutate(Drug="Digoxin") %>%
  bind_rows(fits.df)->
  fits.df


#write_csv("./Derived/metformin.fits.csv")

################################################################################
fit <- glm(age_delta ~. ,data=digoxin.fit.df, family = gaussian)
model.matrix(fit) -> X
X[,colnames(X)[colnames(X)!="(Intercept)"]]  -> X
N <- dim(X)[1]
p <- dim(X)[2]

X <- scale(X)
Y <- digoxin.fit.df$age_delta %>% scale
fit.glmnet.cv <- cv.glmnet(x=X, y=Y, family="gaussian",
                           intercept = TRUE, alpha=1)
fits.glmnet[["Digoxin"]] <- fit.glmnet.cv

#saveRDS(fit.glmnet.cv,"./Derived/metformin.fit.glm.Rds")

################################################################################
BB.fit.df <- dataframes_drugs.df %>%
  filter(Drug=="BB") %>%
  select(-Drug)

get_effect.bin.self.reported(BB.fit.df,add_markers = TRUE) %>%
  mutate(AdditionalMarkers=TRUE) %>%
  bind_rows(get_effect.bin.self.reported(BB.fit.df, add_markers = FALSE) %>%
              mutate(AdditionalMarkers=FALSE)) -> BB.fits.df

BB.fits.df %>%
  mutate(Drug="BB") %>%
  bind_rows(fits.df)->
  fits.df

  #write_csv("./Derived/BB.fits.csv")
################################################################################
fit <- glm(age_delta ~. ,data=BB.fit.df, family = gaussian)
model.matrix(fit) -> X
X[,colnames(X)[colnames(X)!="(Intercept)"]]  -> X
N <- dim(X)[1]
p <- dim(X)[2]

X <- scale(X)
Y <- BB.fit.df$age_delta %>% scale
fit.glmnet.cv <- cv.glmnet(x=X, y=Y, family="gaussian",
                           intercept = TRUE, alpha=1)
fits.glmnet[["BB"]] <- fit.glmnet.cv
#saveRDS(fit.glmnet.cv,"./Derived/BB.fit.glm.Rds")
################################################################################
ACEi.fit.df <- dataframes_drugs.df %>%
  filter(Drug=="ACEi") %>%
  select(-Drug)

get_effect.bin.self.reported(ACEi.fit.df,add_markers = TRUE) %>%
  mutate(AdditionalMarkers=TRUE) %>%
  bind_rows(get_effect.bin.self.reported(ACEi.fit.df, add_markers = FALSE) %>%
              mutate(AdditionalMarkers=FALSE)) -> ACEi.fits.df

ACEi.fits.df %>%
  mutate(Drug="ACEi") %>%
  bind_rows(fits.df)->
  fits.df
  #write_csv("./Derived/ACEi.fits.csv")
################################################################################
fit <- glm(age_delta ~. ,data=ACEi.fit.df, family = gaussian)
model.matrix(fit) -> X
X[,colnames(X)[colnames(X)!="(Intercept)"]]  -> X
N <- dim(X)[1]
p <- dim(X)[2]

X <- scale(X)
Y <- ACEi.fit.df$age_delta %>% scale
fit.glmnet.cv <- cv.glmnet(x=X, y=Y, family="gaussian",
                           intercept = TRUE, alpha=1)
fits.glmnet[["ACEi"]] <- fit.glmnet.cv
#saveRDS(fit.glmnet.cv,"./Derived/ACEi.fit.glm.Rds")
################################################################################
ARB.fit.df <- dataframes_drugs.df %>%
  filter(Drug=="ARB") %>%
  select(-Drug)

get_effect.bin.self.reported(ARB.fit.df,add_markers = TRUE) %>%
  mutate(AdditionalMarkers=TRUE) %>%
  bind_rows(get_effect.bin.self.reported(ARB.fit.df, add_markers = FALSE) %>%
              mutate(AdditionalMarkers=FALSE)) -> ARB.fits.df

ARB.fits.df %>%
  mutate(Drug="ARB") %>%
  bind_rows(fits.df)->
  fits.df
  #write_csv("./Derived/ARB.fits.csv")
################################################################################
fit <- glm(age_delta ~. ,data=ARB.fit.df, family = gaussian)
model.matrix(fit) -> X
X[,colnames(X)[colnames(X)!="(Intercept)"]]  -> X
N <- dim(X)[1]
p <- dim(X)[2]

X <- scale(X)
Y <- ARB.fit.df$age_delta %>% scale
fit.glmnet.cv <- cv.glmnet(x=X, y=Y, family="gaussian",
                           intercept = TRUE, alpha=1)
fits.glmnet[["ARB"]] <- fit.glmnet.cv
#saveRDS(fit.glmnet.cv,"./Derived/ARB.fit.glm.Rds")
################################################################################
CCB.fit.df <- dataframes_drugs.df %>%
  filter(Drug=="CCB") %>%
  select(-Drug)

get_effect.bin.self.reported(CCB.fit.df,add_markers = TRUE) %>%
  mutate(AdditionalMarkers=TRUE) %>%
  bind_rows(get_effect.bin.self.reported(CCB.fit.df, add_markers = FALSE) %>%
              mutate(AdditionalMarkers=FALSE)) -> CCB.fits.df

CCB.fits.df %>%
  mutate(Drug="CCB") %>%
  bind_rows(fits.df)->
  fits.df
  #write_csv("./Derived/CCB.fits.csv")
################################################################################
fit <- glm(age_delta ~. ,data=CCB.fit.df, family = gaussian)
model.matrix(fit) -> X
X[,colnames(X)[colnames(X)!="(Intercept)"]]  -> X
N <- dim(X)[1]
p <- dim(X)[2]

X <- scale(X)
Y <- CCB.fit.df$age_delta %>% scale
fit.glmnet.cv <- cv.glmnet(x=X, y=Y, family="gaussian",
                           intercept = TRUE, alpha=1)
fits.glmnet[["CCB"]] <- fit.glmnet.cv
#saveRDS(fit.glmnet.cv,"./Derived/CCB.fit.glm.Rds")
################################################################################
diuretics.fit.df <- dataframes_drugs.df %>%
  filter(Drug=="Diuretics") %>%
  select(-Drug)

get_effect.bin.self.reported(diuretics.fit.df,add_markers = TRUE) %>%
  mutate(AdditionalMarkers=TRUE) %>%
  bind_rows(get_effect.bin.self.reported(diuretics.fit.df, add_markers = FALSE) %>%
              mutate(AdditionalMarkers=FALSE)) -> diuretics.fits.df

diuretics.fits.df %>%
  mutate(Drug="Diuretics") %>%
  bind_rows(fits.df)->
  fits.df

#write_csv("./Derived/diuretics.fits.csv")
################################################################################
fit <- glm(age_delta ~. ,data=diuretics.fit.df, family = gaussian)
model.matrix(fit) -> X
X[,colnames(X)[colnames(X)!="(Intercept)"]]  -> X
N <- dim(X)[1]
p <- dim(X)[2]

X <- scale(X)
Y <- diuretics.fit.df$age_delta %>% scale
fit.glmnet.cv <- cv.glmnet(x=X, y=Y, family="gaussian",
                           intercept = TRUE, alpha=1)
fits.glmnet[["Diuretics"]] <- fit.glmnet.cv
#saveRDS(fit.glmnet.cv,"./Derived/diuretics.fit.glm.Rds")
################################################################################
################################################################################
fits.df %>%
  write_csv("./Derived/fits.csv")

saveRDS(fits.glmnet, "./Derived/fits.glm.Rds")
