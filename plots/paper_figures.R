### Load functions data
source("functions.R")
source("plots/paper_figure_aesthetics.R")

library(ggplot2)
library(cowplot)
library(dplyr)
library(magrittr)
library(ggpubr)

###############
### FIGURES ###
###############

### Load in corresponding data in order to render plots
# note: all data are for baseline 'low community incidence' scenario

# load in prevalence data for lot 3
load("analysis_CTC/data_analysis_CTC/lot3/df_prevalence_summarized_3ltcfs.Rdata")
df_prevalenceLongSummarizedInfected_3ltcfs$LTCF = factor(df_prevalenceLongSummarizedInfected_3ltcfs$LTCF, labels = vec_3ltcfs_labels_control)

# load in cumulative incidence data for lot 3 with slight adjustment
load("analysis_CTC/data_analysis_CTC/lot3/df_incidence_nosocomial_3ltcfs.Rdata")
df_incidence_nosocomial_3ltcfs$LTCF = factor(df_incidence_nosocomial_3ltcfs$LTCF, labels = vec_3ltcfs_labels_control)

# organize and summarize
dfInc_cumul_onset = df_incidence_nosocomial_3ltcfs%>%
  mutate(incCommAllTotal = incCommCasesTotal + incIndexTotal)%>%
  dplyr::select(-c(immunity, IPC, Network, NetworkIPC, incNosocomialPA, incNosocomialPE, 
                   incIndexPA, incIndexPE, incCommCasesPA, incCommCasesPE, incCommCasesTotal, incIndexTotal))%>%
  pivot_longer(-c(LTCF, SIM), names_to = 'measure', values_to = 'incidence')%>%
  mutate(`infection onset` = ifelse(grepl("Nosocomial", measure), 'hospital', 'community'))
dfInc_cumul_onset_summarized = dfInc_cumul_onset%>%
  group_by(LTCF, `infection onset`)%>%
  summarise(incidence_mean = mean(incidence))

dfInc_cumul_cat = df_incidence_nosocomial_3ltcfs%>%
  mutate(incTotalPA = incNosocomialPA + incCommCasesPA + incIndexPA,
         incTotalPE = incNosocomialPE + incCommCasesPE + incIndexPE)%>%
  dplyr::select(LTCF, SIM, incTotalPA, incTotalPE)%>%
  pivot_longer(-c(LTCF, SIM), names_to = 'measure', values_to = 'incidence')%>%
  mutate(`type of individual` = ifelse(grepl("PA", measure), 'patient', 'staff'))
dfInc_cumul_cat_summarized = dfInc_cumul_cat%>%
  group_by(LTCF, `type of individual`)%>%
  summarise(incidence_mean = mean(incidence))


# load in incidence_t data for lot 3 with slight adjustment
load("analysis_CTC/data_analysis_CTC/lot3/df_incidence_t_3ltcfs.Rdata")
load("analysis_CTC/data_analysis_CTC/lot3/df_incidence_t_3ltcfs_summarized.Rdata")

dfInc_3ltcfs$LTCF = factor(dfInc_3ltcfs$LTCF, labels = vec_3ltcfs_labels_control)
dfInc_3ltcfs_summarized$LTCF = factor(dfInc_3ltcfs_summarized$LTCF, labels = vec_3ltcfs_labels_control)


# incidence by type of individual (combined community- and hospital-onset)
dfInc_3ltcfs_cat = filter(dfInc_3ltcfs, cat %in% c('patients', 'staff'))%>%
  group_by(Date, immunity, IPC, Network, SIM, cat, LTCF, day)%>%
  summarise(Nb= sum(Nb))
dfInc_3ltcfs_cat_summarized = dfInc_3ltcfs_cat%>%
  group_by(Date, immunity, IPC, Network, cat, LTCF, day)%>%
  summarise(Nb_mean = mean(Nb),
            Nb_median = median(Nb),
            Nb_975 = quantile(Nb, 0.975),
            Nb_025 = quantile(Nb, 0.025))

# load in superspreading data 
load("analysis_CTC/data_analysis_CTC/lot3/df_incidence_super_neverspreading_3ltcfs_DENSITY_infectors.Rdata")
load("analysis_CTC/data_analysis_CTC/lot3/df_incidence_super_neverspreading_3ltcfs_DENSITY_infectees.Rdata")

df_incidence_super_neverspreading_frequencies_3ltcfs_DENSITY_infectees$LTCF = factor(df_incidence_super_neverspreading_frequencies_3ltcfs_DENSITY_infectees$LTCF, labels = vec_3ltcfs_labels_control)
df_incidence_super_neverspreading_frequencies_3ltcfs_DENSITY_infectors$LTCF = factor(df_incidence_super_neverspreading_frequencies_3ltcfs_DENSITY_infectors$LTCF, labels = vec_3ltcfs_labels_control)


##################
### PREVALENCE ###
##################


pA_prevalence = ggplot(df_prevalenceLongSummarizedInfected_3ltcfs, aes(x = day, y = mean, fill = status))+
  geom_bar(stat = 'identity', colour = 'black', size = 0.1)+
  ylab('mean infection prevalence')+
  xlab('time (days)')+
  theme_bw()+
  scale_fill_manual("infection status", values= cols_infection_status[c(1,2,3,5)], labels = c('symptomatic', 'asymptomatic', 'pre-symptomatic', 'exposed'))+
  facet_wrap(facets = vars(LTCF))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'bottom')+
  xlim(0.5,15.5)
pA_prevalence


#####################################
### INCIDENCE_T (NOSOCOMIAL ONLY) ###
#####################################

pA_incidence_t_noso = ggplot(dfInc_3ltcfs%>%filter(cat %in% c('patients & staff'), Measure %in% 'hospital'),
                              aes(x = day, y = Nb))+
  geom_point(alpha = 0.1, size = 1, shape = 21, colour = 'white', fill = 'black')+
  #geom_ribbon(alpha = 0.2, colour = NA)+
  geom_line(dfInc_3ltcfs_summarized%>%filter(cat %in% c('patients & staff'), Measure %in% 'hospital'), 
            mapping = aes(y = Nb_mean), lwd = 0.8, alpha = 1, colour = 'red')+
  facet_wrap(facets = vars(LTCF), ncol = 3)+
  theme_bw()+
  ylab('daily infection incidence')+xlab('time (days)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(shape = guide_legend(override.aes = list(alpha=1)))+
  xlim(0.5,15.5)
pA_incidence_t_noso

###########################
### INCIDENCE_T (ONSET) ###
###########################

pA_incidence_t_onset = ggplot(dfInc_3ltcfs%>%filter(cat %in% c('patients & staff')),
                            aes(x = day, y = Nb, group = Measure, colour = Measure, fill = Measure, shape = Measure))+
  geom_point(alpha = 0.1, size = 1.5)+
  #geom_ribbon(alpha = 0.2, colour = NA)+
  geom_line(dfInc_3ltcfs_summarized%>%filter(cat %in% c('patients & staff')), mapping = aes(y = Nb_mean), lwd = 0.8, alpha = 1)+
  facet_wrap(facets = vars(LTCF), ncol = 3)+
  theme_bw()+
  ylab('daily infection incidence')+xlab('time (days)')+
  scale_shape_manual("infection onset", values = c(22,23), labels = c('community', 'nosocomial'))+
  scale_colour_manual("infection onset", values = cols_onset, labels = c('community', 'nosocomial'))+
  scale_fill_manual("infection onset", values = cols_onset, labels = c('community', 'nosocomial'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(shape = guide_legend(override.aes = list(alpha=1)))
pA_incidence_t_onset

#########################
### INCIDENCE_T (CAT) ###
#########################

pA_incidence_t_cat = ggplot(dfInc_3ltcfs_cat,
                          aes(x = day, y = Nb, group = cat, colour = cat, fill = cat, shape = cat))+
  geom_point(alpha = 0.1, size = 1.5)+
  geom_line(dfInc_3ltcfs_cat_summarized%>%filter(!cat %in% c('patients & staff')), mapping = aes(y = Nb_mean), lwd = 0.8, alpha = 1)+
  #geom_ribbon(alpha = 0.2, colour = NA)+
  #geom_line(stat = 'identity', size = 1, alpha = 1)+
  facet_wrap(facets = vars(LTCF), ncol = 3)+
  theme_bw()+
  ylab('daily infection incidence')+xlab('time (days)')+
  scale_shape_manual("type of individual", values = c(24,21), labels = c('staff', 'patient'))+
  scale_colour_manual("type of individual", values = cols_classes[c(2:3)], labels = c('staff', 'patient'))+
  scale_fill_manual("type of individual", values = cols_classes[c(2:3)], labels = c('staff', 'patient'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(shape = guide_legend(override.aes = list(alpha=1)))
pA_incidence_t_cat

##########################################
### INCIDENCE_T (RAINBOW FOR EACH SIM) ###
##########################################
incidence_label = df_incidence_nosocomial_3ltcfs%>%group_by(LTCF)%>%
  summarise(I = mean(incNosocomialTotal),
            Imin = min(incNosocomialTotal),
            Imax = max(incNosocomialTotal))

pA_incidence_t_noso_rainbow = ggplot(dfInc_3ltcfs%>%filter(cat %in% c('patients & staff'), Measure %in% 'hospital'),
                             aes(x = day, y = Nb, group = SIM, colour = factor(SIM)))+
  geom_line(alpha = 0.6, size = 0.1)+
  #geom_ribbon(alpha = 0.2, colour = NA)+
  geom_line(dfInc_3ltcfs_summarized%>%filter(cat %in% c('patients & staff'), Measure %in% 'hospital'), 
            mapping = aes(y = Nb_mean, group = NA), lwd = 1, alpha = 1, colour = 'black')+
  facet_wrap(facets = vars(LTCF), ncol = 3)+
  theme_bw()+
  ylab('infection incidence')+xlab('time (days)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(shape = guide_legend(override.aes = list(alpha=1)), colour = "none")+
  xlim(0.5,15.5)+
  geom_text(incidence_label%>%mutate(SIM = 100), x = 12, y = 16, 
            mapping = aes(label = paste0("I = ",round(I,1), " (", Imin, " - ", Imax,")")), colour = 'black')
pA_incidence_t_noso_rainbow

################
### Figure 1 ###
################
# NB: table describing interventions in place for each LTCF (Figure 1A) made in PPT

pA = plot_grid(pA_prevalence, pA_incidence_t_noso_rainbow, nrow = 2, align = 'v', axis = 'lr', labels = c('B', 'C'))
pA

#################
### Figure S2 ###
#################
### Cumulative incidence (onset)
pA_incidence_cumul_onset = ggplot(dfInc_cumul_onset, aes(x = incidence, fill = `infection onset`))+
  geom_histogram(stat = 'bin', binwidth = 2, position = 'identity', alpha = 0.6)+
  #scale_alpha_manual(values = c(0.8,0.2))+
  theme_bw()+
  ylab('number of simulations')+xlab('cumulative infection incidence (at 2 weeks)')+
  facet_wrap(facets = vars(LTCF), scales = 'free_y', ncol = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(dfInc_cumul_onset_summarized, mapping = aes(xintercept = incidence_mean, colour = `infection onset`), linetype = 2)+
  scale_fill_manual(values = cols_onset, labels = c('community', 'nosocomial'))+
  scale_colour_manual(values = cols_onset, labels = c('community', 'nosocomial'))+
  scale_x_continuous(breaks = seq(0,80,10))
pA_incidence_cumul_onset

### Cumulative incidence (type of individual)
pA_incidence_cumul_cat = ggplot(dfInc_cumul_cat, aes(x = incidence, fill = `type of individual`))+
  geom_histogram(stat = 'bin', binwidth = 2, position = 'identity', alpha = 0.6)+
  #scale_alpha_manual(values = c(0.8,0.2))+
  theme_bw()+
  ylab('number of simulations')+xlab('cumulative infection incidence (at 2 weeks)')+
  facet_wrap(facets = vars(LTCF), scales = 'free_y', ncol = 3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values = cols_classes[c(3,2)], guide = guide_legend(reverse = TRUE) )+
  geom_vline(dfInc_cumul_cat_summarized, mapping = aes(xintercept = incidence_mean, colour = `type of individual`), linetype = 2)+
  scale_colour_manual(values = cols_classes[c(3,2)], guide = guide_legend(reverse = TRUE) )+
  scale_x_continuous(breaks = seq(0,60,10))
pA_incidence_cumul_cat


#################
### Figure S3 ###
#################
### Super-spreading

### infectors
pA_superspreading_infectors = ggplot(df_incidence_super_neverspreading_frequencies_3ltcfs_DENSITY_infectors, 
                                     aes(x = x,
                                         ymin = ymin,
                                         ymax = ymax,
                                         group = magnitude,
                                         fill = magnitude))+
  geom_ribbon(colour = 'black', alpha = 1)+
  facet_grid(cols = vars(LTCF))+
  # Enclosing lines
  #geom_line(aes(y = ymin))+
  #geom_line(aes(y = ymax))+
  theme_bw()+
  xlab("proportion of infected individuals")+
  ylab("")+
  scale_fill_manual(name = "contribution to nosocomial transmission", labels = vec_superspreading_labels_long, 
                    values =  c("orange", "violet", "darkgrey"),
                    guide = guide_legend(reverse = T))+
  scale_y_continuous(breaks = 1:3, labels = c("low spreader", "super spreader", "non spreader"))+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0", "0.25", "0.5", "0.75", "1"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pA_superspreading_infectors

### infectees
pA_superspreading_infectees = ggplot(df_incidence_super_neverspreading_frequencies_3ltcfs_DENSITY_infectees, 
                                     aes(x = x,
                                         ymin = ymin,
                                         ymax = ymax,
                                         group = magnitude,
                                         fill = magnitude))+
  geom_ribbon(colour = 'black', alpha = 1)+
  facet_grid(cols = vars(LTCF))+
  # Enclosing lines
  #geom_line(aes(y = ymin))+
  #geom_line(aes(y = ymax))+
  theme_bw()+
  xlab("proportion of nosocomial transmission events")+
  ylab("")+
  scale_fill_manual(name = "source of nosocomial acquisition", labels = vec_superspreading_labels_long, 
                    values =  c("orange", "violet"),
                    guide = guide_legend(reverse = T))+
  scale_y_continuous(breaks = 1:2, labels = c("low spreader", "super spreader"))+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0", "0.25", "0.5", "0.75", "1"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pA_superspreading_infectees

### combined
p_superspreading = plot_grid(pA_superspreading_infectors, pA_superspreading_infectees, nrow = 2, labels = c('A', 'B'),  align = 'v', axis = 'lr')
p_superspreading


#################
### Figure S4 ###
#################
### Test sensitivity curves

df_sens_adjusted = read.csv("kucirka_adjusted.csv")

df_sens_adjusted$Test = factor(df_sens_adjusted$Test,
                               levels = c('sens', 'sens_adj_rel', 'sens_adj_const'),
                               labels = c('RT-PCR', 'Ag-RDT (A)', 'Ag-RDT (B)'))

pSensCurves = ggplot(df_sens_adjusted, aes(x = day, y = sensitivity*100, ymin = 0, ymax = sensitivity, colour = Test, shape = Test))+
  geom_point()+
  geom_line()+
  theme_classic()+
  ylab("probability of positive test result (%)")+
  xlab("days since SARS-CoV-2 exposure")+
  scale_colour_manual("type of test", values = cols_typescreening[c(1,3,2)])+
  scale_shape_manual("type of test", values = c(15,19,19))
pSensCurves


########################################################
### PCR DYNAMICS (LOW COMMUNITY INCIDENCE SCENARIO)  ###
########################################################

load("analysis_CTC/data_analysis_CTC/lot3/df_sensitivityMeans_timevarying.Rdata")
load("Surveillance/output/Outcomes_lot3/dataIncidenceBeforeAfterBaseline.Rdata")
load("Surveillance/output/Outcomes_lot3/dataIncidenceBeforeAfterBaseline_Means.Rdata")

df_sensitivityMeans_timevarying$LTCF = factor(df_sensitivityMeans_timevarying$LTCF, labels = vec_3ltcfs_labels_control)
dataIncidenceBeforeAfterBaseline$LTCF = factor(dataIncidenceBeforeAfterBaseline$LTCF, labels = vec_3ltcfs_labels_control)
dataIncidenceBeforeAfterBaseline_Means$LTCF = factor(dataIncidenceBeforeAfterBaseline_Means$LTCF, labels = vec_3ltcfs_labels_control)

### sensitivity, time-varying
df_sensitivityMeans_timevarying = df_sensitivityMeans_timevarying%>%
  mutate(time = abs(as.numeric(as.Date("2009-07-05") - as.Date(date))))%>%
  filter(onset != 'none')

df_sensitivityMeans_timevarying_Means = df_sensitivityMeans_timevarying%>%
  group_by(time, LTCF, onset)%>%
  dplyr::summarise(sensitivity_mean_means = mean(sensitivity_mean, na.rm = T))

#################
### Figure S5 ###
#################
### Sensitivity dynamics 

pB_sensitivity = ggplot(df_sensitivityMeans_timevarying, 
                        aes(x = time, y = as.numeric(sensitivity_mean)*100, 
                            group = interaction(simulationCTC, onset), colour = onset))+
  geom_line(size = 0.08, alpha = 0.4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(facets = vars(LTCF), ncol = 3)+
  ylab("mean probability of\npositive RT-PCR result (%)")+xlab("time (days)")+
  scale_colour_manual("infection onset", values = cols_onset, labels = c('community', 'nosocomial'))+
  geom_line(df_sensitivityMeans_timevarying_Means, mapping = aes(y = sensitivity_mean_means*100, group = onset), size = 1)+
  geom_line(df_sensitivityMeans_timevarying_Means,
            mapping = aes(y = sensitivity_mean_means*100, group = onset), size = 1, colour = 'black', linetype = 3)
pB_sensitivity


#################
### Figure S6 ###
#################
### Cumul incidence with and without routine testing

pB_incidenceBeforeAfterBaseline = ggplot(dataIncidenceBeforeAfterBaseline, aes(x = Incidence, fill = Testing))+
  geom_histogram(stat = 'bin', binwidth = 3, position = 'identity', alpha = 0.7)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(facets = vars(LTCF), ncol = 3, scales = 'free_y')+
  ylab("number of simulations")+xlab("cumulative nosocomial incidence (at 2 weeks)")+
  geom_vline(dataIncidenceBeforeAfterBaseline_Means, mapping = aes(xintercept = Incidence_mean, colour = Testing), linetype = 2)+
  scale_fill_manual("testing strategy", values = cols_basePCR, labels = c("routine RT-PCR testing", "no testing"))+
  scale_colour_manual("testing strategy", values = cols_basePCR, labels = c("routine RT-PCR testing", "no testing"))
pB_incidenceBeforeAfterBaseline



#########################################
### UNRAVELLING PROP AVERTED OUTCOME  ###
#########################################
# NB: not in manuscript but an example of raw data, 
# illustrative of underlying binomial infection prevention outcome,
# relatively uninformative nature of its quantiles,
# and calculation of bootstrap confidence intervals

load("Surveillance/output/Outcomes_lot3/dataIncidenceBeforeAfterBaseline.Rdata")
load("Surveillance/output/Outcomes_lot3/dataIncidenceBeforeAfterBaseline_Means.Rdata")
load("Surveillance/output/Outcomes_lot3/dataBootstrapDemo_fullD.Rdata")
load("Surveillance/output/Outcomes_lot3/dataBootstrapDemo_CIs.Rdata")

# fix LTCF labels for all data
dataIncidenceBeforeAfterBaseline$LTCF = factor(dataIncidenceBeforeAfterBaseline$LTCF, labels = vec_3ltcfs_labels_control)
dataIncidenceBeforeAfterBaseline_Means$LTCF = factor(dataIncidenceBeforeAfterBaseline_Means$LTCF, labels = vec_3ltcfs_labels_control)
d_all = d_all%>%mutate(LTCF = ifelse(LTCF == 1, vec_3ltcfs_labels_control[1],
                             ifelse(LTCF == 2, vec_3ltcfs_labels_control[2], vec_3ltcfs_labels_control[3])))
ci_all = ci_all%>%mutate(LTCF = ifelse(LTCF == 1, vec_3ltcfs_labels_control[1],
                                     ifelse(LTCF == 2, vec_3ltcfs_labels_control[2], vec_3ltcfs_labels_control[3])))

### raw efficacy of PCR
dataIncidenceBeforeAfterBaseline_propAverted = dataIncidenceBeforeAfterBaseline%>%
  pivot_wider(names_from = "Testing", values_from = "Incidence")%>%
  mutate(propAverted = 1-(`Baseline PCR testing`/`No testing`))

dataIncidenceBeforeAfterBaseline_propAverted_Summarized = dataIncidenceBeforeAfterBaseline_propAverted%>%
  group_by(simulationCTC, LTCF)%>%
  dplyr::summarise(propAverted_mean = mean(propAverted))

dataIncidenceBeforeAfterBaseline_propAverted_Summarized_Summarized = dataIncidenceBeforeAfterBaseline_propAverted_Summarized%>%
  group_by(LTCF)%>%
  dplyr::summarise(propAverted_mean_mean = mean(propAverted_mean, na.rm = T))

dataPropAvertedRawPCR = left_join(dataIncidenceBeforeAfterBaseline_propAverted, dataIncidenceBeforeAfterBaseline_propAverted_Summarized)%>%
  dplyr::select(-c(`No testing`, `Baseline PCR testing`))%>%
  pivot_longer(-c(simulationCTC, testingRound, LTCF), values_to = 'propAverted', names_to = 'measure')


# factor SIM ID
labels_rawDataPlot = c("individual estimates", "means per outbreak")

pC_rawPropAvertedPCR = ggplot(dataPropAvertedRawPCR, aes(x = propAverted, y = simulationCTC, shape = measure, alpha = measure, size = measure,
                                                         colour = measure, fill = factor(simulationCTC)))+
  geom_point()+
  scale_alpha_manual(values = c(0.1,1), labels = labels_rawDataPlot)+
  scale_colour_manual(values = c('white', 'black'), labels = labels_rawDataPlot)+
  scale_shape_manual(values = c(21, 23), labels = labels_rawDataPlot)+
  scale_size_manual(values = c(1.5,3), labels = labels_rawDataPlot)+
  scale_linetype_manual(values = c(2,2,2))+
  facet_grid(cols = vars(LTCF))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(fill = F, colour = F, alpha = F)+
  scale_y_reverse(breaks = c(1, 25, 50, 75, 100))+
  xlab("proportion of nosocomial infections averted by routine RT-PCR testing")+
  ylab("outbreak ID")+
  geom_vline(dataIncidenceBeforeAfterBaseline_propAverted_Summarized_Summarized,
             mapping = aes(xintercept = propAverted_mean_mean), linetype = 2)
pC_rawPropAvertedPCR

### showing distribution of outcomes and including boxplots
pC_distributionPropAvertedPCR = ggplot(dataIncidenceBeforeAfterBaseline_propAverted%>%mutate(count = -250), 
       aes(x = propAverted))+
  geom_histogram(bins = 50, fill = 'darkgrey', colour = 'darkgrey')+
  facet_wrap(facets = vars(LTCF), ncol = 3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_boxplot(mapping = aes(y = count), width = 200, fill = 'lightgrey')+
  geom_hline(yintercept = 0)+
  ylab('number of simulations')+
  xlab("proportion of nosocomial infections averted by routine RT-PCR testing")+
  geom_vline(dataIncidenceBeforeAfterBaseline_propAverted_Summarized_Summarized,
             mapping = aes(xintercept = propAverted_mean_mean), linetype = 2)
pC_distributionPropAvertedPCR


### showing distributions of bootstrap means and calculated CIs
vec_ylims = c(-75, -150, -225)

pC_boots = ggplot(d_all, 
       aes(x = proportion, fill = factor(nboots), colour = factor(nboots), alpha = factor(nboots)))+
  geom_histogram(bins = 50, position = 'identity', colour = NA)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab('density')+xlab('mean proportion of nosocomial infections averted by routine RT-PCR testing')+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  geom_point(data = ci_all, 
             mapping = aes(x = proportion, colour = factor(nboots)), 
             y = rep(vec_ylims,3),
             inherit.aes = F)+
  facet_wrap(facets = vars(LTCF), ncol = 3, scales = 'free_x')+
  geom_hline(yintercept = 0)+
  ylim(-225, 900)+
  geom_rect(data = ci_all, mapping = aes(xmin = min, xmax = max), ymin = rep(vec_ylims, 3), ymax = rep(vec_ylims, 3))+
  scale_fill_manual("# of bootstrap\nreplicates", values = cols_bootstrap, labels = expression(10^2,10^3,10^4))+
  scale_colour_manual("# of bootstrap\nreplicates", values = cols_bootstrap, labels = expression(10^2,10^3,10^4))+
  scale_alpha_manual("# of bootstrap\nreplicates", values = rep(0.4,3), labels = expression(10^2,10^3,10^4))
pC_boots


################
### Figure 2 ###
################

### FIGURE 2A
# Figure 2 is an illustrative figure demonstrating how our surveillance algorithm works
# for simplicity, two additional community-onset cases introduced over simulation time (PE-201-LOR, PE-706-PAH) are ommitted from the transmission chain in Figure 2A
# thus evaluating only the transmission chain resulting from the original 4 index cases
# Note: the transmission chain was illustrated by hand, but data informing it are provided below

### FOR LOW INCIDENCE, IPC 1, NETWORK 1, SIMULATION 22
fileroot = f_fileroot("lot3", "0.20.2",1,22,1)
dataStatus = f_read_csv(filepath_CTC, fileroot, "statusByDay")
dataSecondaryCases = f_read_csv(filepath_CTC, fileroot, "secondaryCases")
dataIndex = f_read_csv(filepath_CTC, fileroot, "FirstIndex")
dataCommCases = f_read_csv(filepath_CTC, fileroot, 'CommCases')

### PCR sensitivty (for Figure 2A)
dataSensitivityPCR = f_dataSensitivity(dataStatus, "PCR")%>%as.data.frame()%>%
  filter(id %in% c('PE-212-BOM', 'PA-637-VIJ', 'PE-015-DAP'))
# infections with symptoms: 
# PE 212 on July 10: 0.613
# PA 637 on July 14: 0.752
# PE-015 on July 6: 0.799

### Ag-RDT sensitivity (for Figure 2A)
dataSensitivityRDT = f_dataSensitivity(dataStatus, "RDT2")%>%as.data.frame()%>%
  filter(id %in% c('PE-212-BOM', 'PA-616-TAJ', 'PE-015-DAP'))
# screening on day 2, July 7:
# PE-212: 0
# PA-616: 0.659
# PE-015: 0.611



### FIGURE 2B
### load data for all 100 testing runs for a single CTC outbreak, compare infections averted by different strategies

load("Surveillance/output/Output_prepared_lot3/m_surveillance_prepared_0202_1_1.Rdata")

data_all = dataSurveillanceOutput0202_1_1%>%
  filter(stratSurveillance %in% c('adm_sym', 'd02', 'adm_sym_d02'), 
         typeDailyScreening == "RDT2", 
         targetDailyScreening == "All",
         sensitivity_function == "time-varying")


# ommission of introduction cases (see above) required a small modification for Figure 2B: there was one nosocomial infection on the final simulation day 14 from an introduction: PE-201-LOR transmitted to PA-451-BRS
# PE-201-LOR was infected too late to be detected by screening, and was pre-symptomatic at the end of simulation, so not detected by symptomatic testing either
# however, on 2 surveillance runs PE-201-LOR had non-covid but covid-like symptoms, received PCR, tested positive and was isolated, preventing infection of PA-451-BRS
# on both of these runs, all 22 nosocomial infections were prevented; so this outcome had to be modified to 21 in order to exclude PA-451-BRS's infection from Figure 2B, as it was excluded from Figure 2A

df_exploreSurveillance = filter(data_all, simulationCTC == 22)%>%
  mutate(infAverted = ifelse(infAverted == 13, 12,
                             ifelse(infAverted == 22, 21, infAverted)))

df_exploreSurveillance_both = df_exploreSurveillance%>%filter(stratSurveillance == 'adm_sym_d02')

df_exploreSurveillance_order = df_exploreSurveillance_both[order(df_exploreSurveillance_both$infAverted),]$testingRound

df_exploreSurveillance$testingRound = factor(df_exploreSurveillance$testingRound, levels = df_exploreSurveillance_order)
df_exploreSurveillance$stratSurveillance = factor(df_exploreSurveillance$stratSurveillance,
                                levels = c('adm_sym', 'd02', 'adm_sym_d02'),
                                labels = c('routine RT-PCR testing', 'Ag-RDT screening (day 2)', 'routine RT-PCR testing +\nAg-RDT screening (day 2)'))

p_exploreSurveillance = ggplot(df_exploreSurveillance, aes(x = testingRound, y = infAverted, colour = stratSurveillance))+
  geom_point(alpha = 0.75)+
  theme_classic()+
  ylab('infections averted\ndue to testing and isolation')+
  xlab('stochastic surveillance run')+
  scale_x_discrete(breaks = df_exploreSurveillance_order[c(1,25,50,75,100)], labels = c(1,25,50,75,100))+
  scale_colour_manual('surveillance strategy', values = c('#00B0F0', '#FF847B', '#D25FFF'))
p_exploreSurveillance




##########################
### SCREENING EFFICACY ### 
##########################

### PLOT: OVERALL EFFICACY ACROSS ALL STRATEGIES
# Not provided in manscript but left here for completeness

# load in prevalence data for lot 3
load("Surveillance/output/Outcomes_lot3/dataPropAverted_3ltcfs_compareCIs.Rdata")

dataPropAverted_3ltcfs_compareCIs_final = dataPropAverted_3ltcfs_compareCIs%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2', sensitivity_function == 'time-varying', methodCI == 'bootstrap')

dataPropAverted_3ltcfs_compareCIs_final$stratSurveillance = factor(dataPropAverted_3ltcfs_compareCIs_final$stratSurveillance, 
                                                             levels = all_strats_3ltcfs_labels[c(3,2,1,4:29)],
                                                             labels = all_strats_3ltcfs_labels_UPDATED[c(3,2,1,4:29)])

dataPropAverted_3ltcfs_compareCIs_final$LTCF = factor(dataPropAverted_3ltcfs_compareCIs_final$LTCF, labels = vec_3ltcfs_labels_control)

pD_propAvertedAll_total = ggplot(dataPropAverted_3ltcfs_compareCIs_final,
                                 aes(y = stratSurveillance, x = central, xmin = lower, xmax = upper, 
                                     colour = targetDailyScreening, shape = targetDailyScreening))+
  geom_errorbar(width = 0.5)+
  geom_point(size = 1)+
  theme_bw()+
  facet_grid(cols = vars(LTCF))+
  xlab('proportion of nosocomial infections averted')+ylab('')+
  scale_colour_manual("screening target", values = cols_classes)+
  scale_shape_manual("screening target", values = c(8, 17, 19, 15))+
  scale_y_discrete(limits = rev)+
  geom_hline(yintercept = c(8.5, 17.5, 26.5), linetype = 2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlim(0,0.8)
pD_propAvertedAll_total



### PLOT: REFINED EFFICACY ACROSS DIFFERENT SCREENING MODALITIES
# Overlay different timings, distinguishing between "1-round Ag-RDT", "routine RT-PCR + 1-round Ag-RDT", "routine RT-PCR + 2-round Ag-RDT"

dataPropAverted_3ltcfs_compareCIs_final_refined = dataPropAverted_3ltcfs_compareCIs_final%>%
  filter(!stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[1:3])%>%
  mutate(modality = case_when(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[4:12] ~ "1-round Ag-RDT",
                              stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[13:21] ~ "routine RT-PCR + 1-round Ag-RDT",
                              stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[22:29] ~ "routine RT-PCR + 2-round Ag-RDT"),
         timing = case_when(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(4,13)] ~ 1,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(5,14,22)] ~ 2,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(6,15,23)] ~ 3,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(7,16,24)] ~ 4,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(8,17,25)] ~ 5,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(9,18,26)] ~ 6,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(10,19,27)] ~ 7,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(11,20,28)] ~ 8,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(12,21,29)] ~ 9
                            ))

dataPropAverted_3ltcfs_compareCIs_final_refined_baselinePCR = dataPropAverted_3ltcfs_compareCIs_final%>%
  filter(stratSurveillance == all_strats_3ltcfs_labels_UPDATED[3])%>%
  slice(rep(1:n(), each=9))%>%
  mutate(timing = row_number()%%9,
         timing = ifelse(timing == 0, 9, timing))

################
### FIGURE 3 ###
################

pD_propAvertedAll_modality = ggplot(dataPropAverted_3ltcfs_compareCIs_final_refined%>%filter(targetDailyScreening == 'patients & staff'),
       aes(x = timing, y = central*100, ymin = lower*100, ymax = upper*100, 
           colour = modality))+#, shape = targetDailyScreening))+
  geom_errorbar(width = 0.5)+
  geom_point(size = 1)+
  theme_bw()+
  facet_grid(cols = vars(LTCF))+
  xlab('screening date (days since outbreak detection)')+
  ylab('reduction in nosocomial incidence (%)')+
  scale_x_continuous(breaks = c(1:9))+
  scale_colour_manual("screening implementation", values = cols_implementation)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylim(15,80)+
  geom_hline(data = dataPropAverted_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(yintercept = central*100), alpha = 0.5)+
  geom_hline(data = dataPropAverted_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(yintercept = lower*100), alpha = 0.3, linetype = 2)+
  geom_hline(data = dataPropAverted_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(yintercept = upper*100), alpha = 0.3, linetype = 2)+
  geom_text(data = dataPropAverted_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(y = lower*100-1.2, x = 7.5, label = "routine RT-PCR"), 
            colour = 'black', size = 2.5)
pD_propAvertedAll_modality


#################
### Figure S8 ###
#################

# load in efficacy data (propAverted) for lot 3
load("Surveillance/output/Outcomes_lot3/dataPropAverted_3ltcfs_compareCIs.Rdata")
dataPropAverted_3ltcfs_compareCIs$LTCF = factor(dataPropAverted_3ltcfs_compareCIs$LTCF, labels = vec_3ltcfs_labels_control)

### pE: different targets (patients, staff, both)
dat_pE_targets = dataPropAverted_3ltcfs_compareCIs%>%filter(typeDailyScreening %in% 'RDT2',
                                                                    methodCI == "bootstrap",
                                                                    sensitivity_function == "time-varying",
                                                                    stratSurveillance %in% all_strats_3ltcfs_labels[c(13,22:29)])
dat_pE_targets$stratSurveillance = factor(dat_pE_targets$stratSurveillance, levels = all_strats_3ltcfs_labels[c(13,22:29)])


pE_targets = ggplot(dat_pE_targets, 
                    aes(x = stratSurveillance, y = central*100, ymin = lower*100, ymax = upper*100, 
                        colour = targetDailyScreening, fill = targetDailyScreening, group = targetDailyScreening, shape = targetDailyScreening))+
  geom_errorbar(width = 0.5)+
  geom_point(size = 1.5)+
  #geom_line()+
  #geom_ribbon(alpha = 0.2, colour = NA)+
  theme_bw()+
  facet_wrap(facets = vars(LTCF))+
  xlab('screening date (days since outbreak detection)')+
  ylab('reduction in nosocomial incidence (%)')+
  scale_colour_manual("screening target", values = cols_classes[c(2,3,1)])+
  scale_fill_manual("screening target", values = cols_classes[c(2,3,1)])+
  scale_shape_manual("screening target", values = c(17, 19, 8))+
  scale_x_discrete(labels = c('1', paste0('1\n&\n', 2:9)))+
  theme(axis.title.x = element_text(size = 10),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
pE_targets


### pE: different types of screening (PCR vs 2 RDTs)
dat_pE_typescreening = dataPropAverted_3ltcfs_compareCIs%>%filter(targetDailyScreening == "patients & staff",
                                                            methodCI == "bootstrap",
                                                            sensitivity_function == "time-varying",
                                                            stratSurveillance %in% all_strats_3ltcfs_labels[c(13,22:29)])
dat_pE_typescreening$stratSurveillance = factor(dat_pE_typescreening$stratSurveillance, levels = all_strats_3ltcfs_labels[c(13,22:29)])


pE_typescreening = ggplot(dat_pE_typescreening, 
                          aes(x = stratSurveillance, y = central*100, ymin = lower*100, ymax = upper*100, 
                              colour = typeDailyScreening, fill = typeDailyScreening, group = typeDailyScreening, shape = typeDailyScreening))+
  geom_errorbar(width = 0.5)+
  geom_point(size = 1.5)+
  #geom_line()+
  #geom_ribbon(alpha = 0.2, colour = NA)+
  theme_bw()+
  facet_wrap(facets = vars(LTCF))+
  xlab('screening date (days since outbreak detection)')+
  ylab('reduction in nosocomial incidence (%)')+
  scale_colour_manual("type of test\nused for screening", values = cols_typescreening, labels = c("RT-PCR", "Ag-RDT (B)", "Ag-RDT (A)"))+
  scale_fill_manual("type of test\nused for screening", values = cols_typescreening, labels = c("RT-PCR", "Ag-RDT (B)", "Ag-RDT (A)"))+
  scale_shape_manual("type of test\nused for screening", values = c(15, 19, 8), labels = c("RT-PCR", "Ag-RDT (B)", "Ag-RDT (A)"))+
  theme(axis.title.x = element_text(size = 10),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  guides(colour = guide_legend(reverse = TRUE),
         fill = guide_legend(reverse = TRUE),
         shape = guide_legend(reverse = TRUE))+
  scale_x_discrete(labels = c('1', paste0('1\n&\n', 2:9)))
pE_typescreening


### pE: different sensitivity functions (time-varying, uniform, perfect)
dat_pE_sensfunc = dataPropAverted_3ltcfs_compareCIs%>%filter(targetDailyScreening == "patients & staff",
                                                             typeDailyScreening %in% 'RDT2',
                                                                  methodCI == "bootstrap",
                                                                  stratSurveillance %in% all_strats_3ltcfs_labels[c(13,22:29)])
dat_pE_sensfunc$stratSurveillance = factor(dat_pE_sensfunc$stratSurveillance, levels = all_strats_3ltcfs_labels[c(13,22:29)])

pE_sensfunc = ggplot(dat_pE_sensfunc, 
                     aes(x = stratSurveillance, y = central*100, ymin = lower*100, ymax = upper*100, 
                         colour = sensitivity_function, fill = sensitivity_function, group = sensitivity_function, shape = sensitivity_function))+
  geom_errorbar(width = 0.5)+
  geom_point(size = 1.5)+
  #geom_line()+
  #geom_ribbon(alpha = 0.2, colour = NA)+
  theme_bw()+
  facet_wrap(facets = vars(LTCF))+
  xlab('screening date (days since outbreak detection)')+
  ylab('reduction in nosocomial incidence (%)')+
  scale_colour_manual("sensitivity function", values = cols_sensfuncs)+
  scale_fill_manual("sensitivity function", values = cols_sensfuncs)+
  scale_shape_manual("sensitivity function", values = c(3, 8, 4))+
  theme(axis.title.x = element_text(size = 10),
        #axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_x_discrete(labels = c('1', paste0('1\n&\n', 2:9)))
pE_sensfunc





#########################################
### SCREENING EFFICACY OTHER OUTCOMES ###
#########################################

### SCREENING EFFICACY 2: number of infections averted and incidence rate ratio

# load in prevalence data for lot 3
load("Surveillance/output/Outcomes_lot3/dataInfAverted_3ltcfs_compareCIs.Rdata")
load("Surveillance/output/Outcomes_lot3/dataIRR_3ltcfs_compareCIs.Rdata")

dataInfAverted_3ltcfs_compareCIs$LTCF = factor(dataInfAverted_3ltcfs_compareCIs$LTCF, labels = vec_3ltcfs_labels_control)
dataIRR_3ltcfs_compareCIs$LTCF = factor(dataIRR_3ltcfs_compareCIs$LTCF, labels = vec_3ltcfs_labels_control)

### freshen up data as needed for efficacy plots
dataInfAverted_3ltcfs_compareCIs_final = dataInfAverted_3ltcfs_compareCIs%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2', sensitivity_function == 'time-varying', methodCI == 'bootstrap')

dataInfAverted_3ltcfs_compareCIs_final$stratSurveillance = factor(dataInfAverted_3ltcfs_compareCIs_final$stratSurveillance, 
                                                                   levels = all_strats_3ltcfs_labels[c(3,2,1,4:29)],
                                                                   labels = all_strats_3ltcfs_labels_UPDATED[c(3,2,1,4:29)])

dataIRR_3ltcfs_compareCIs_final = dataIRR_3ltcfs_compareCIs%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2', sensitivity_function == 'time-varying', methodCI == 'bootstrap')

dataIRR_3ltcfs_compareCIs_final$stratSurveillance = factor(dataIRR_3ltcfs_compareCIs_final$stratSurveillance, 
                                                                  levels = all_strats_3ltcfs_labels[c(3,2,1,4:29)],
                                                                  labels = all_strats_3ltcfs_labels_UPDATED[c(3,2,1,4:29)])

### PLOT: REFINED EFFICACY ACROSS DIFFERENT SCREENING MODALITIES
# Overlay different timings, distinguishing between "1-round Ag-RDT", "routine RT-PCR + 1-round Ag-RDT", "routine RT-PCR + 2-round Ag-RDT"

# Infaverted:
dataInfAverted_3ltcfs_compareCIs_final_refined = dataInfAverted_3ltcfs_compareCIs_final%>%
  filter(!stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[1:3])%>%
  mutate(modality = case_when(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[4:12] ~ "1-round Ag-RDT",
                              stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[13:21] ~ "routine RT-PCR + 1-round Ag-RDT",
                              stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[22:29] ~ "routine RT-PCR + 2-round Ag-RDT"),
         timing = case_when(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(4,13)] ~ 1,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(5,14,22)] ~ 2,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(6,15,23)] ~ 3,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(7,16,24)] ~ 4,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(8,17,25)] ~ 5,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(9,18,26)] ~ 6,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(10,19,27)] ~ 7,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(11,20,28)] ~ 8,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(12,21,29)] ~ 9
         ))

dataInfAverted_3ltcfs_compareCIs_final_refined_baselinePCR = dataInfAverted_3ltcfs_compareCIs_final%>%
  filter(stratSurveillance == all_strats_3ltcfs_labels_UPDATED[3])%>%
  slice(rep(1:n(), each=9))%>%
  mutate(timing = row_number()%%9,
         timing = ifelse(timing == 0, 9, timing))

# and IRRs:
dataIRR_3ltcfs_compareCIs_final_refined = dataIRR_3ltcfs_compareCIs_final%>%
  filter(!stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[1:3])%>%
  mutate(modality = case_when(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[4:12] ~ "1-round Ag-RDT",
                              stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[13:21] ~ "routine RT-PCR + 1-round Ag-RDT",
                              stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[22:29] ~ "routine RT-PCR + 2-round Ag-RDT"),
         timing = case_when(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(4,13)] ~ 1,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(5,14,22)] ~ 2,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(6,15,23)] ~ 3,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(7,16,24)] ~ 4,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(8,17,25)] ~ 5,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(9,18,26)] ~ 6,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(10,19,27)] ~ 7,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(11,20,28)] ~ 8,
                            stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(12,21,29)] ~ 9
         ))

dataIRR_3ltcfs_compareCIs_final_refined_baselinePCR = dataIRR_3ltcfs_compareCIs_final%>%
  filter(stratSurveillance == all_strats_3ltcfs_labels_UPDATED[3])%>%
  slice(rep(1:n(), each=9))%>%
  mutate(timing = row_number()%%9,
         timing = ifelse(timing == 0, 9, timing))

#################
### FIGURE S7 ###
#################

pF_infAvertedAll_modality = ggplot(dataInfAverted_3ltcfs_compareCIs_final_refined%>%filter(targetDailyScreening == 'patients & staff'),
                                    aes(x = timing, y = central, ymin = lower, ymax = upper, 
                                        colour = modality))+#, shape = targetDailyScreening))+
  geom_errorbar(width = 0.5)+
  geom_point(size = 1)+
  theme_bw()+
  facet_grid(cols = vars(LTCF))+
  xlab('screening date (days since outbreak detection)')+
  ylab('number of nosocomial infections averted')+
  scale_x_continuous(breaks = c(1:9))+
  scale_colour_manual("screening implementation", values = cols_implementation)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(data = dataInfAverted_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(yintercept = central), alpha = 0.5)+
  geom_hline(data = dataInfAverted_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(yintercept = lower), alpha = 0.3, linetype = 2)+
  geom_hline(data = dataInfAverted_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(yintercept = upper), alpha = 0.3, linetype = 2)
pF_infAvertedAll_modality



pF_IRR_modality = ggplot(dataIRR_3ltcfs_compareCIs_final_refined%>%filter(targetDailyScreening == 'patients & staff'),
                                   aes(x = timing, y = central, ymin = lower, ymax = upper, 
                                       colour = modality))+#, shape = targetDailyScreening))+
  geom_errorbar(width = 0.5)+
  geom_point(size = 1)+
  theme_bw()+
  facet_grid(cols = vars(LTCF))+
  xlab('screening date (days since outbreak detection)')+
  ylab('incidence rate ratio\n(relative to no testing or screening)')+
  scale_x_continuous(breaks = c(1:9))+
  scale_colour_manual("screening implementation", values = cols_implementation)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(data = dataIRR_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(yintercept = central), alpha = 0.5)+
  geom_hline(data = dataIRR_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(yintercept = lower), alpha = 0.3, linetype = 2)+
  geom_hline(data = dataIRR_3ltcfs_compareCIs_final_refined_baselinePCR, 
             mapping = aes(yintercept = upper), alpha = 0.3, linetype = 2)
pF_IRR_modality



##########################
### TESTING EVALUATION ### 
##########################

load("Surveillance/output/Outcomes_lot3/dataTestPositivity_3ltcfs_compareCIs_plottable.Rdata")
load("Surveillance/output/Outcomes_lot3/dataTestSensitivity_3ltcfs_compareCIs_plottable.Rdata")
load("Surveillance/output/Outcomes_lot3/dataTestSpecificity_3ltcfs_compareCIs_plottable.Rdata")
load("Surveillance/output/Outcomes_lot3/dataTestNPV_3ltcfs_compareCIs_plottable.Rdata")
load("Surveillance/output/Outcomes_lot3/dataTestPPV_3ltcfs_compareCIs_plottable.Rdata")
load("Surveillance/output/Outcomes_lot3/dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable.Rdata")
load("Surveillance/output/Outcomes_lot3/dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable.Rdata")

dataTestPositivity_3ltcfs_compareCIs_plottable$LTCF = factor(dataTestPositivity_3ltcfs_compareCIs_plottable$LTCF, labels = vec_3ltcfs_labels_control)
dataTestSensitivity_3ltcfs_compareCIs_plottable$LTCF = factor(dataTestSensitivity_3ltcfs_compareCIs_plottable$LTCF, labels = vec_3ltcfs_labels_control)
dataTestSpecificity_3ltcfs_compareCIs_plottable$LTCF = factor(dataTestSpecificity_3ltcfs_compareCIs_plottable$LTCF, labels = vec_3ltcfs_labels_control)
dataTestNPV_3ltcfs_compareCIs_plottable$LTCF = factor(dataTestNPV_3ltcfs_compareCIs_plottable$LTCF, labels = vec_3ltcfs_labels_control)
dataTestPPV_3ltcfs_compareCIs_plottable$LTCF = factor(dataTestPPV_3ltcfs_compareCIs_plottable$LTCF, labels = vec_3ltcfs_labels_control)
dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable$LTCF = factor(dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable$LTCF, labels = vec_3ltcfs_labels_control)
dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable$LTCF = factor(dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable$LTCF, labels = vec_3ltcfs_labels_control)

### Positivity
dataTestPositivity_3ltcfs_compareCIs_plottable_final = dataTestPositivity_3ltcfs_compareCIs_plottable%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2',
         sensitivity_function == 'time-varying', methodCI == 'bootstrap',
         !stratSurveillance %in% c('Symptoms (routine PCR)', 'Admissions (routine PCR)'))
dataTestPositivity_3ltcfs_compareCIs_plottable_final$stratSurveillance = factor(dataTestPositivity_3ltcfs_compareCIs_plottable_final$stratSurveillance,
                                                                                 levels = rev(all_strats_3ltcfs_labels[c(3:29)]),
                                                                                 labels = rev(all_strats_3ltcfs_labels_UPDATED[c(3:29)]))

### Sensitivity
dataTestSensitivity_3ltcfs_compareCIs_plottable_final = dataTestSensitivity_3ltcfs_compareCIs_plottable%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2',
         sensitivity_function == 'time-varying', methodCI == 'bootstrap',
         !stratSurveillance %in% c('Symptoms (routine PCR)', 'Admissions (routine PCR)'))
dataTestSensitivity_3ltcfs_compareCIs_plottable_final$stratSurveillance = factor(dataTestSensitivity_3ltcfs_compareCIs_plottable_final$stratSurveillance,
                                                                                 levels = rev(all_strats_3ltcfs_labels[c(3:29)]),
                                                                                 labels = rev(all_strats_3ltcfs_labels_UPDATED[c(3:29)]))

### Specificity
dataTestSpecificity_3ltcfs_compareCIs_plottable_final = dataTestSpecificity_3ltcfs_compareCIs_plottable%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2',
         sensitivity_function == 'time-varying', methodCI == 'bootstrap',
         !stratSurveillance %in% c('Symptoms (routine PCR)', 'Admissions (routine PCR)'))
dataTestSpecificity_3ltcfs_compareCIs_plottable_final$stratSurveillance = factor(dataTestSpecificity_3ltcfs_compareCIs_plottable_final$stratSurveillance,
                                                                                 levels = rev(all_strats_3ltcfs_labels[c(3:29)]),
                                                                                 labels = rev(all_strats_3ltcfs_labels_UPDATED[c(3:29)]))

### NPV
dataTestNPV_3ltcfs_compareCIs_plottable_final = dataTestNPV_3ltcfs_compareCIs_plottable%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2',
         sensitivity_function == 'time-varying', methodCI == 'bootstrap',
         !stratSurveillance %in% c('Symptoms (routine PCR)', 'Admissions (routine PCR)'))
dataTestNPV_3ltcfs_compareCIs_plottable_final$stratSurveillance = factor(dataTestNPV_3ltcfs_compareCIs_plottable_final$stratSurveillance,
                                                                                 levels = rev(all_strats_3ltcfs_labels[c(3:29)]),
                                                                                 labels = rev(all_strats_3ltcfs_labels_UPDATED[c(3:29)]))

### PPV
dataTestPPV_3ltcfs_compareCIs_plottable_final = dataTestPPV_3ltcfs_compareCIs_plottable%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2',
         sensitivity_function == 'time-varying', methodCI == 'bootstrap',
         !stratSurveillance %in% c('Symptoms (routine PCR)', 'Admissions (routine PCR)'))
dataTestPPV_3ltcfs_compareCIs_plottable_final$stratSurveillance = factor(dataTestPPV_3ltcfs_compareCIs_plottable_final$stratSurveillance,
                                                                         levels = rev(all_strats_3ltcfs_labels[c(3:29)]),
                                                                         labels = rev(all_strats_3ltcfs_labels_UPDATED[c(3:29)]))

### cases detected per test
dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable_final = dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2',
         sensitivity_function == 'time-varying', methodCI == 'bootstrap',
         !stratSurveillance %in% c('Symptoms (routine PCR)', 'Admissions (routine PCR)'))
dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable_final$stratSurveillance = factor(dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable_final$stratSurveillance,
                                                                         levels = rev(all_strats_3ltcfs_labels[c(3:29)]),
                                                                         labels = rev(all_strats_3ltcfs_labels_UPDATED[c(3:29)]))

### cases averted per test
dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable_final = dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable%>%
  filter(targetDailyScreening == 'no screening' | typeDailyScreening == 'RDT2',
         sensitivity_function == 'time-varying', methodCI == 'bootstrap',
         !stratSurveillance %in% c('Symptoms (routine PCR)', 'Admissions (routine PCR)'))
dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable_final$stratSurveillance = factor(dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable_final$stratSurveillance,
                                                                                      levels = rev(all_strats_3ltcfs_labels[c(3:29)]),
                                                                                      labels = rev(all_strats_3ltcfs_labels_UPDATED[c(3:29)]))
#################
### Figure S9 ###
#################

### 1. PCR plot: TPV, NPV, PPV, NPV, infections detected/1000; infections averted/1000 {LTCFs as colours, outcomes as panels}
### 2. Ag-RDT plot: TPV, NPV, PPV, NPV, infections detected/1000; infections averted/1000
### 3. Fig 4: for two-round screening, compare infections detected and averted per 1,000 Ag-RDT tests (panel per outcome), stratified by cat
### 4. Fig 5: Cost ratio
pH_pointsize = 1.5

### 1. PCR performance
dataEfficiency_PCR = rbind(dataTestSensitivity_3ltcfs_compareCIs_plottable_final%>%mutate(central = central * 100, lower = lower * 100, upper = upper * 100,
                                                                                          outcome = 'true positive rate (%)'),
                           dataTestSpecificity_3ltcfs_compareCIs_plottable_final%>%mutate(central = central * 100, lower = lower * 100, upper = upper * 100,
                                                                                          outcome = 'true negative rate (%)'),
                           dataTestPPV_3ltcfs_compareCIs_plottable_final%>%mutate(central = central * 100, lower = lower * 100, upper = upper * 100,
                                                                                  outcome = 'positive predictive value (%)'),
                           dataTestNPV_3ltcfs_compareCIs_plottable_final%>%mutate(central = central * 100, lower = lower * 100, upper = upper * 100,
                                                                                  outcome = 'negative predictive value (%)'),
                           dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable_final%>%mutate(central = central * 1000, lower = lower * 1000, upper = upper * 1000,
                                                                                               outcome = 'cases detected/1,000 RT-PCR tests'),
                           dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable_final%>%mutate(central = central * 1000, lower = lower * 1000, upper = upper * 1000,
                                                                                              outcome = 'cases averted/1,000 RT-PCR tests'))%>%
  filter(stratSurveillance == "routine RT-PCR (admissions + symptoms)")

dataEfficiency_PCR$LTCF = factor(dataEfficiency_PCR$LTCF, levels = rev(vec_3ltcfs_labels_control))
dataEfficiency_PCR$outcome = factor(dataEfficiency_PCR$outcome,
                                    levels = c('true positive rate (%)',
                                               'true negative rate (%)',
                                               'positive predictive value (%)',
                                               'negative predictive value (%)',
                                               'cases detected/1,000 RT-PCR tests',
                                               'cases averted/1,000 RT-PCR tests'))

pPerformancePCR = ggplot(dataEfficiency_PCR, aes(x = central, xmin = lower, xmax = upper, y = LTCF, colour = LTCF))+
  geom_point(size = pH_pointsize)+
  geom_errorbar(width = 0.3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(facets = vars(outcome), scales = 'free_x', ncol = 2, nrow = 3)+
  ylab('')+
  xlab('surveillance outcome (routine RT-PCR testing)')+
  scale_y_discrete(labels = c("LTCF 3", "LTCF 2", "LTCF 1"))+
  scale_colour_manual("", values = rev(cols_ltcf), guide = guide_legend(reverse = T))
pPerformancePCR

##################
### Figure S#0 ###
##################

### 2. Ag-RDT performance (sens, spec, NPV, PPV)
# 2a sensitivity
pPerformanceAgRDTa = ggplot(dataTestSensitivity_3ltcfs_compareCIs_plottable_final%>%
         filter(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(13,22:29)]),
       aes(y = stratSurveillance, x = central*100, xmin = lower*100, xmax = upper*100, 
           colour = targetDailyScreening, shape = targetDailyScreening))+#, shape = targetDailyScreening))+
  geom_point(size = pH_pointsize)+
  geom_errorbar(width = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(cols = vars(LTCF))+
  ylab('screening date')+
  xlab('true positive rate (%)')+
  scale_y_discrete(labels = c(paste0('1 & ',9:2,''),  '1'))+
  scale_colour_manual("screening target", values = cols_classes[2:4])+
  scale_shape_manual("screening target", values = shapes_cat)

# 2b specificity
pPerformanceAgRDTb = ggplot(dataTestSpecificity_3ltcfs_compareCIs_plottable_final%>%
         filter(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(13,22:29)]),
       aes(y = stratSurveillance, x = central*100, xmin = lower*100, xmax = upper*100, 
           colour = targetDailyScreening, shape = targetDailyScreening))+#, shape = targetDailyScreening))+
  geom_point(size = pH_pointsize)+
  geom_errorbar(width = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(cols = vars(LTCF))+
  ylab('screening date')+
  xlab('true negative rate (%)')+
  scale_y_discrete(labels = c(paste0('1 & ',9:2,''),  '1'))+
  scale_x_continuous(breaks = c(99.69, 99.7, 99.71))+
  scale_colour_manual("screening target", values = cols_classes[2:4])+
  scale_shape_manual("screening target", values = shapes_cat)

# 2c PPV
pPerformanceAgRDTc = ggplot(dataTestPPV_3ltcfs_compareCIs_plottable_final%>%
         filter(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(13,22:29)]),
       aes(y = stratSurveillance, x = central*100, xmin = lower*100, xmax = upper*100, 
           colour = targetDailyScreening, shape = targetDailyScreening))+#, shape = targetDailyScreening))+
  geom_point(size = pH_pointsize)+
  geom_errorbar(width = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(cols = vars(LTCF))+
  ylab('screening date')+
  xlab('positive predictive value (%)')+
  scale_y_discrete(labels = c(paste0('1 & ',9:2,''),  '1'))+
  scale_colour_manual("screening target", values = cols_classes[2:4])+
  scale_shape_manual("screening target", values = shapes_cat)

# 2d NPV
pPerformanceAgRDTd = ggplot(dataTestNPV_3ltcfs_compareCIs_plottable_final%>%
         filter(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(13,22:29)]),
       aes(y = stratSurveillance, x = central*100, xmin = lower*100, xmax = upper*100, 
           colour = targetDailyScreening, shape = targetDailyScreening))+#, shape = targetDailyScreening))+
  geom_point(size = pH_pointsize)+
  geom_errorbar(width = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(cols = vars(LTCF))+
  ylab('screening date')+
  xlab('negative predictive value (%)')+
  scale_y_discrete(labels = c(paste0('1 & ',9:2,''),  '1'))+
  scale_colour_manual("screening target", values = cols_classes[2:4])+
  scale_shape_manual("screening target", values = shapes_cat)

# 2e cases detected/1,000
pPerformanceAgRDTe = ggplot(dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable_final%>%
                              filter(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(13,22:29)]),
                            aes(y = stratSurveillance, x = central*1000, xmin = lower*1000, xmax = upper*1000, 
                                colour = targetDailyScreening, shape = targetDailyScreening))+#, shape = targetDailyScreening))+
  geom_point(size = pH_pointsize)+
  geom_errorbar(width = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(cols = vars(LTCF))+
  ylab('screening date')+
  xlab('cases detected/1,000 Ag-RDT tests')+
  scale_y_discrete(labels = c(paste0('1 & ',9:2,''),  '1'))+
  scale_colour_manual("screening target", values = cols_classes[2:4])+
  scale_shape_manual("screening target", values = shapes_cat)

pPerformanceAgRDTf = ggplot(dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable_final%>%
                              filter(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[c(13,22:29)]),
                            aes(y = stratSurveillance, x = central*1000, xmin = lower*1000, xmax = upper*1000, 
                                colour = targetDailyScreening, shape = targetDailyScreening))+#, shape = targetDailyScreening))+
  geom_point(size = pH_pointsize)+
  geom_errorbar(width = 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(cols = vars(LTCF))+
  ylab('screening date')+
  xlab('cases averted/1,000 Ag-RDT tests')+
  scale_y_discrete(labels = c(paste0('1 & ',9:2,''),  '1'))+
  scale_colour_manual("screening target", values = cols_classes[2:4])+
  scale_shape_manual("screening target", values = shapes_cat)


pPerformanceAgRDT_all = ggarrange(pPerformanceAgRDTa, pPerformanceAgRDTb,pPerformanceAgRDTc, pPerformanceAgRDTd, pPerformanceAgRDTe, pPerformanceAgRDTf,
          ncol= 1, common.legend = T, legend = "right", labels = c('A', 'B', 'C', 'D', 'E', 'F'))
pPerformanceAgRDT_all


################
### Figure 4 ###
################

pos = position_dodge(0.5)
p2roundAgRDTdetected = ggplot(dataCasesDetectedPerTest_3ltcfs_compareCIs_plottable_final%>%
                                filter(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[25]),
                              aes(x = LTCF, y = central*1000, ymin = lower*1000, ymax = upper*1000, 
                                  fill = targetDailyScreening, shape = targetDailyScreening))+#, shape = targetDailyScreening))+
  geom_histogram(stat = 'identity', position = pos, colour = 'black', width = 0.5)+
  geom_errorbar(position = pos, width = 0.1)+
  theme_minimal()+
  geom_hline(yintercept = 0)+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab('')+
  ylab('apparent efficiency\n(cases detected / 1,000 Ag-RDT tests)')+
  scale_x_discrete(labels = c('LTCF 1\n(low control)', 'LTCF 2\n(moderate control)', 'LTCF 3\n(high control)'))+
  scale_fill_manual("screening target", values = cols_classes[c(2:4)])+
  scale_shape_manual("screening target", values = shapes_cat)
p2roundAgRDTdetected

p2roundAgRDTaverted = ggplot(dataCasesAvertedPerTest_3ltcfs_compareCIs_plottable_final%>%
                                filter(stratSurveillance %in% all_strats_3ltcfs_labels_UPDATED[25]),
                              aes(x = LTCF, y = central*1000, ymin = lower*1000, ymax = upper*1000, 
                                  fill = targetDailyScreening, shape = targetDailyScreening))+#, shape = targetDailyScreening))+
  geom_histogram(stat = 'identity', position = pos, colour = 'black', width = 0.5)+
  geom_errorbar(position = pos, width = 0.1)+
  theme_minimal()+
  geom_hline(yintercept = 0)+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab('')+
  ylab('marginal real efficiency\n(cases averted / 1,000 Ag-RDT tests)')+
  scale_x_discrete(labels = c('LTCF 1\n(low control)', 'LTCF 2\n(moderate control)', 'LTCF 3\n(high control)'))+
  scale_fill_manual("screening target", values = cols_classes[c(2:4)])+
  scale_shape_manual("screening target", values = shapes_cat)
p2roundAgRDTaverted

p2roundAgRDTdetectedAverted = ggarrange(p2roundAgRDTdetected+ylim(0,20.5), p2roundAgRDTaverted+ylim(0,20.5), 
          ncol = 2, common.legend = T, legend = 'right', labels = c('A','B'))
p2roundAgRDTdetectedAverted


################
### Figure 5 ###
################

load("Surveillance/output/Outcomes_lot3/dataCostsPerCaseAverted_3ltcfs_plottable.Rdata")

# vary PCR cost
p4econ_pcr = ggplot(dataCostPerCaseAverted_plottable%>%filter(cost_rdt == 5),
                 aes(x = factor(cost_pcr), y = cost_per_caseAverted, ymin = cost_per_caseAverted_lower, ymax = cost_per_caseAverted_upper, 
                     colour = stratSurveillance, fill = stratSurveillance))+
  facet_wrap(facets = vars(LTCF), nrow = 1)+
  geom_histogram(stat = 'identity', position = "dodge", colour = 'black', size = 0.2, width = 0.7)+
  geom_errorbar(position = "dodge", colour = 'black', width = 0.7)+
  theme_bw()+
  geom_hline(yintercept = 0)+
  ylab('cost-effectiveness (€ / case averted)')+
  xlab('unit cost per RT-PCR test (€)')+
  theme(legend.direction = 'horizontal', legend.position = 'bottom')+
  scale_fill_manual("surveillance strategy", values = cols_strats3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 8))+ 
  scale_y_continuous(labels = function(x) gsub("000$", "k", x)) +
  guides(fill = guide_legend(nrow = 2))
p4econ_pcr


# vary Ag-RDT cost
p4econ_rdt = ggplot(dataCostPerCaseAverted_plottable%>%filter(cost_pcr == 50),
                    aes(x = factor(cost_rdt), y = cost_per_caseAverted, ymin = cost_per_caseAverted_lower, ymax = cost_per_caseAverted_upper, 
                        colour = stratSurveillance, fill = stratSurveillance))+
  facet_wrap(facets = vars(LTCF), nrow = 1)+
  geom_histogram(stat = 'identity', position = "dodge", colour = 'black', size = 0.2, width = 0.7)+
  geom_errorbar(position = "dodge", colour = 'black', width = 0.7)+
  theme_bw()+
  geom_hline(yintercept = 0)+
  ylab('cost-effectiveness (€ / case averted)')+
  xlab('unit cost per Ag-RDT test (€)')+
  theme(legend.direction = 'horizontal', legend.position = 'bottom')+
  scale_fill_manual("surveillance strategy", values = cols_strats3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 8))+ 
  scale_y_continuous(labels = function(x) gsub("000$", "k", x)) +
  guides(fill = guide_legend(nrow = 2))
p4econ_rdt


p4econ_rdt_pcr = ggarrange(p4econ_rdt, p4econ_pcr, nrow = 2, common.legend = T, legend = 'bottom', labels = c('A', 'B'))
p4econ_rdt_pcr
