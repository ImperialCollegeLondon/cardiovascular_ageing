library(tidyverse)
library(kableExtra)
################################################################################
################################################################################
get_df.bin.self.reported <- function(p.df,df) {
  df %>%
    left_join(p.df %>% mutate(PriorMed=1),by="SID") %>%
    mutate(PriorMed=if_else(is.na(PriorMed), 0, PriorMed)) %>%
    select(SID,
           age_delta,
           Sex,
           Age,
           packyrs,
           alcoholintakegpd,
           #z_adj_2, # telomere
           PriorMed,
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
}
################################################################################
################################################################################
get_effect.bin.self.reported <- function(tmp, add_markers=TRUE) {

  #tmp <- get_df.bin.self.reported(p.df)



  fracPriorMed <- tmp %>%
    pull(PriorMed) %>%
    mean

  if(add_markers) {
    lm(age_delta~Sex+poly(Age, 2)+packyrs+alcoholintakegpd+PriorMed+bmi+Obesity+CAD+
         Hypertension+ `Diabetes Mellitus`+Hyperchol+`Heart Failure`+sbp+dbp+pulse_rate ,
       data=tmp) -> fit
  } else {
    lm(age_delta~Sex+poly(Age, 2)+packyrs+alcoholintakegpd+PriorMed+bmi+Obesity+CAD+
         Hypertension+ `Diabetes Mellitus`+Hyperchol+`Heart Failure`,
       data=tmp) -> fit
  }


  fit %>%
    broom::glance() %>%
    pivot_longer(everything(), names_to="Estimator", values_to="Estimate") ->
    fit.glance

  fit%>%
    broom::tidy(conf.int=TRUE) %>%
    mutate(fracPriorMed=fracPriorMed) %>%
    mutate(r.squared=fit.glance %>% filter(Estimator=="r.squared") %>% pull(Estimate),
           adj.r.squared=fit.glance %>% filter(Estimator=="adj.r.squared") %>% pull(Estimate),
           nobs=fit.glance %>% filter(Estimator=="nobs") %>% pull(Estimate)
    )
}
################################################################################
################################################################################
generate_table.self.reported <- function(tmp, drugname) {


  tmp %>%
    drop_na(age_delta,
            Sex,
            Age,
            packyrs,
            alcoholintakegpd,
            #z_adj_2, # telomere
            PriorMed,
            Hypertension,
            Obesity,
            `Diabetes Mellitus`,
            CAD,
            Hyperchol,
            `Heart Failure`,
            bmi,
            dbp,
            sbp,
            pulse_rate
    ) %>%
    select(age_delta,
           Sex,
           Age,
           packyrs,
           alcoholintakegpd,
           #z_adj_2, # telomere
           PriorMed,
           Hypertension,
           Obesity,
           `Diabetes Mellitus`,
           CAD,
           Hyperchol,
           `Heart Failure`,
           dbp,
           sbp,
           pulse_rate,
           bmi) %>%
    mutate(PriorMed=if_else(PriorMed==1, drugname,paste0("No ",drugname))) %>%
    gtsummary::tbl_summary(
      by=PriorMed,
      label  = list(
        age_delta   ~ "Cardiac Age gap",
        CAD ~ "Coronary Artery Disease",
        Hyperchol~"Hypercholesterolemia",
        packyrs~ "Packs year of smoking",
        alcoholintakegpd~ "Alcohol is g/day consumed"
      ),
      missing_text = "Missing",
      digits = list(age_delta ~ c(1, 1))

    ) %>%
    gtsummary::add_p()
}
################################################################################
################################################################################
get_gp_sids <- function(prescp.df, df0) {
  df0 %>%
    mutate(SID=as.integer(SID)) %>%
    select(dt_visit,SID) %>%
    left_join(prescp.df,by="SID") ->
    prescp.df

  prescp.df %>%
    drop_na(issue_date) %>%
    mutate(delta=lubridate::time_length(dt_visit-issue_date, unit="month")) %>%
    filter(delta>0) ->
    tmp


  tmp%>%
    group_by(SID) %>%
    filter(delta<=3) %>%
    ungroup() %>%
    select(SID) %>%
    left_join(tmp , by="SID") %>%
    filter(delta<=6) %>%
    group_by(SID) %>%
    summarise(N=n()) %>%
    filter(N>=2) %>%
    select(SID)


}
################################################################################
################################################################################
get_pretty_fit_table <- function(drug.fits.df, drug_name="Statins") {
  drug.fits.df %>%
    mutate(Adjustment=if_else(AdditionalMarkers==TRUE, "With SBP/DBP/HR", "Without SBP/DBP/HR")) %>%
    select(-fracPriorMed,-r.squared,-adj.r.squared, -Drug,-AdditionalMarkers,-nobs) %>%
    mutate(across(estimate:conf.high	, ~round(.,digits=3))) %>%
    arrange(term, Adjustment) %>%
    kbl(caption = sprintf("%s: Coefficients from linear models",drug_name)) %>%
    kable_classic(full_width = F, html_font = "Cambria")
}
################################################################################
################################################################################
get_pretty_fit_forest_plot <- function(drug.fits.df) {
  drug.fits.df %>%
    filter(term != "(Intercept)") %>%
    mutate(Adjustment=if_else(AdditionalMarkers==TRUE, "With SBP/DBP/HR", "Without SBP/DBP/HR")) %>%
    mutate(Adjustment=factor(Adjustment, levels=c("Without SBP/DBP/HR","With SBP/DBP/HR"), ordered=TRUE)) %>%
    ggplot(aes(x=Adjustment, y=estimate, ymin=conf.low, ymax=conf.high))+
    geom_point()+
    geom_linerange()+
    coord_flip()+
    facet_wrap(~term, scales="free_x",nrow=5)+
    geom_hline(yintercept=0, color="red", alpha=.5, linetype="dashed")
}
################################################################################
################################################################################
get_stratified_age_gap_histogram <- function(drug.df, drug_name) {
  drug.df %>%
    mutate(PriorMed=if_else(PriorMed==0, sprintf("No %s",drug_name), drug_name)) %>%
    group_by(PriorMed) %>%
    mutate(Mean=mean(age_delta),
           Median=median(age_delta)
    ) %>%
    ungroup() %>%
    ggplot(aes(x=age_delta))+
    geom_histogram(bins=30)+
    geom_vline(aes(xintercept=Median), color="green", alpha=.5)+
    facet_wrap(~PriorMed,nrow=2,scales="free_y")
}
################################################################################
################################################################################
get_stratified_age_gap_quantiles <- function(drug.df, drug_name) {
  drug.df %>%
    mutate(PriorMed=if_else(PriorMed==0, sprintf("No %s",drug_name), drug_name)) %>%
    group_by(PriorMed) %>%
    summarise(`1%`=quantile(age_delta,prob=.01) %>% round(2),
              `2.5%`=quantile(age_delta, prob=0.025)%>% round(2),
              `25%`=quantile(age_delta, prob=0.25)%>% round(2),
              `50%`=quantile(age_delta, prob=0.5)%>% round(2),
              `75%`=quantile(age_delta, prob=0.75)%>% round(2),
              `97.5%`=quantile(age_delta, prob=0.975)%>% round(2),
              `99%`=quantile(age_delta, prob=0.99)%>% round(2),
              N=n()
    )
}
################################################################################
################################################################################
get_fit_glance <- function(drug.fits.df) {
  drug.fits.df %>%
    group_by(AdditionalMarkers) %>%
    summarise(N=dplyr::first(nobs),
              fracPriorMed=dplyr::first(fracPriorMed),
              adj.r.squared=dplyr::first(adj.r.squared),
              r.squared=dplyr::first(r.squared)
    ) %>%
    mutate(across(fracPriorMed:r.squared,~round(.,3)))
}
