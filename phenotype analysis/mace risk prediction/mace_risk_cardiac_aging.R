##############################################################################
## 1) Generation of cumulative incidence curves

## assumption: data is a dataset with the following columns
## - age_delta_quantile: quartile of observation (1,2,3 or 4)
## - FUT: time to MACE event or censoring
## - Event: indicator: 1 MACE event, 0 no MACE event observed
## - DI: indicator: 1 all-cause death observed, 0 no all-cause death observed
## - TTD: time to death event or censoring 


## take 99% quartile as cutoff for result display
cutoff <- quantile(data$FUT, 0.975, na.rm=T)

## generate competing event objects
data$StatusComp <- data$Event
data$MaceComp <- data$FUT
for (i in 1:nrow(data)) {
  if (data$DI[i]==1 & data$Event[i]!=1) {
    data$StatusComp[i] <- 2
    data$MaceComp[i] <- data$TTD[i]
  }
}

## calculate curves
ci_fit_mri <- cuminc(data$MaceComp, data$StatusComp, data$age_delta_quantile)


## generate figure
ci_fit_mri_plot <- 
  ci_fit_mri %>% 
  list_modify("Tests" = NULL) %>% 
  map_df(`[`, c("time", "est"), .id = "id") %>% 
  filter(id %in% c("1 1", "2 1","3 1","4 1")) %>% 
  mutate(Group = recode(
    id, 
    "1 1" = "Quartile 1", 
    "2 1" = "Quartile 2",
    "3 1" = "Quartile 3",
    "4 1" = "Quartile 4")
  )
ci_fit_mri_plot$time <- ci_fit_mri_plot$time/365
y_cutoff <- 0.1

pdf(file="cumulative_incidence_curves.pdf", width=8, height=5)
plot_mri <- 
  ggplot(ci_fit_mri_plot, aes(x = time, y = est, color = Group)) +
  geom_step(lwd = 1.2)  +
  ylim(c(0, y_cutoff)) +
  theme_classic() +
  theme(plot.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "bottom") +
  labs(x = "Time to event since imaging visit (years)", 
       y = "Cumulative incidence",
       title = "") 
print(plot_mri)
dev.off()



##############################################################################
## 2) Model-based assessment

## this assumes a data frame with the following columns
## - age_delta: cardiovascular age delta as described in the paper
## - Sex
## - Age
## - Alcohol intake (alcoholintakegpd)
## - Smoking
## - Obesity
## - coronary artery disease (CAD)
## - Hypercholesterolaemia (Hyperchol)
## - Diabetes
## - prior: MACE event prior to baseline? (binary)


## age delta as a cont variable
SurvObj <- survival::Surv(data$FUT, data$Event)
full_MRI <- coxph(SurvObj~age_delta+Sex+Age+smoking+alcoholintakegpd+Hypertension+Obesity+CAD+Hyperchol+diabetes+strata(prior), data=data)
min_MRI <- coxph(SurvObj~age_delta+Sex+Age+strata(prior), data=data)

result_plot_mri <- data.frame(Lower=NA, Upper=NA, Mean=NA, Model=c("Fully adjusted", "Minimally adjusted"))
result_plot_mri[1,1:2] <- confint(full_MRI)[1,]
result_plot_mri[2,1:2] <- confint(min_MRI)[1,]
result_plot_mri[,3] <- rowMeans(result_plot_mri[,1:2])

pcont = ggplot(data=result_plot_mri,
               aes(x = Model,y = Mean, ymin = Lower, ymax = Upper ))+
  geom_pointrange(aes(col=Model))+
  geom_hline(yintercept =1, linetype=2)+
  scale_color_manual(values=c("grey59", "grey3"))+
  xlab('')+ ylab("Hazard ratio per 1 year increase")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper,col=Model),width=1,cex=1)+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y.left = element_text(hjust=0,vjust = 1,angle=0,face="bold"))+
  coord_flip()


## age delta as quartiles (4 vs 1)
result_plot_mri_quad <- data.frame(Lower=NA, Upper=NA, Mean=NA, Model=c("Fully adjusted", "Minimally adjusted"))

data$quartiles_mri_low_high <- data$age_delta_quantile
data$quartiles_mri_low_high[which(data$quartiles_mri_low_high==2 | data$quartiles_mri_low_high==3)] <- NA # ignore middle quartiles
full_MRI_quad <- coxph(SurvObj~quartiles_mri_low_high+Sex+Age+packyrs_instance2+alcoholintakegpd+Hypertension+Obesity+CAD+Hyperchol+diabetes+strata(prior), data=data)
min_MRI_quad <- coxph(SurvObj~quartiles_mri_low_high+Sex+Age+strata(prior), data=data)


result_plot_mri_quad[1,1:2] <- confint(full_MRI_quad)[1,]
result_plot_mri_quad[2,1:2] <- confint(min_MRI_quad)[1,]
result_plot_mri_quad[,3] <- rowMeans(result_plot_mri_quad[,1:2])

pquad = ggplot(data=result_plot_quad,
               aes(x = Model,y = Mean, ymin = Lower, ymax = Upper ))+
  geom_pointrange(aes(col=Model))+
  geom_hline(yintercept =1, linetype=2)+
  scale_color_manual(values=c("grey40", "grey3", "grey70"))+
  xlab('')+ ylab("Hazard ratio (Quartile 1 vs. Quartile 4)")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper,col=Model),width=1,cex=1)+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y.left = element_text(hjust=0,vjust = 1,angle=0,face="bold"))+
  coord_flip()

## export fig
pdf(file="quartile_comparison.pdf", width=8, height=3)
p <- ggarrange(pquad, pcont, 
               align='h', labels=c('A', 'B'),
               common.legend = T)
p
dev.off()

