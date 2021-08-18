library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# load raw data from Kucirka et al. Ann Intern Med 2020
df_kucirka_raw = read.csv2("kucirka_raw.csv")

# take median FNR, remove percent sign, add as new column, divide by 100 and calculate sensitivity as 1-fnr
fnr_raw = str_split(as.character(df_kucirka_raw$False.negative.rate..95..CI.),'%', simplify=T)[,1]
fnr_numeric = as.numeric(fnr_raw)/100
df_kucirka_raw$fnr = fnr_numeric

# clean dataframe for kucirka data
df_kucirka_fnr = df_kucirka_raw%>%
  mutate(sens = 1-fnr)%>%
  dplyr::select(day, fnr, sens)

# extrapolate past 21 days using exponential decay model
model = lm(sens ~ exp(-day*0.055), data = filter(df_kucirka_fnr, day > 10))
sens_predicted = predict(model, data.frame(day = 22:43))%>%as.data.frame(col.names = c('sens'))
colnames(sens_predicted) = 'sens'

# extrapolated data
df_kucirka_extrapolated = sens_predicted%>%
  mutate(day = row_number()+21)%>%
  mutate(sens = ifelse(sens<0,0,sens))%>%
  mutate(fnr = 1-sens)

# combine original and extrapolated data
df_kucirka = rbind(df_kucirka_fnr%>%mutate(Source = 'Kucirka et al.'), df_kucirka_extrapolated%>%mutate(Source = 'Extrapolated')%>%rbind(.,data.frame(sens = 0.37, day = 21, fnr = 0.63, Source = "Extrapolated")))


###############################################
### Antigen rapid diagnostic tests (Ag-RDT) ###
###############################################

### Calculate sensitivity relative to RT-PCR
# data entered manually from meta-analysis from Brummer et al. medrxiv 2020

# Ag-RDT: constant reduction over full infectious period (C)
C_median = 0.738;

# Ag-RDT: less reduction <7 days from symptom onset (Early), compared to >=7 days (Late) (E vs L)
E_median = 0.875; L_median = 0.641; 

# incubation period: assume 5 days
incubation_period = 5

# distance function to minize to find shape parameters
func_optim = function(par){
  
  gamma = par[1]
  beta = par[2]
  
  dat = df_kucirka%>%
    mutate(days_from_peak = day-8)%>%
    mutate(sens_adj_rel = sens*(1-gamma)*exp(-abs(days_from_peak^2)*beta))
  
  dat_early = filter(dat, day < incubation_period + 7)
  dat_late = filter(dat, day >= incubation_period + 7)
  
  # what are the distances?
  distance_overall = C_median - sum(dat$sens_adj_rel)/sum(dat$sens)
  
  distance_early = E_median - sum(dat_early$sens_adj_rel)/sum(dat_early$sens)
  
  distance_late = L_median - sum(dat_late$sens_adj_rel)/sum(dat_late$sens)
  
  # final distance
  distances_sum_sq = abs(distance_overall) + abs(distance_early) + abs(distance_late)
  
  return(distances_sum_sq)
}

# find shape parameters by minimizing distance function func_optim
result_optim = optim(par = c(0,0.1), func_optim)
gamma = result_optim$par[1]
beta = result_optim$par[2]

# calculate Ag-RDT sensitivity using ovrall estimate from Brummer et al. (absolute reduction) and shape parameters (relative reduction)
df_sens = df_kucirka%>%
  mutate(days_from_peak = day-8)%>%
  mutate(sens_adj_const = sens*C_median)%>%
  mutate(sens_adj_rel = sens*(1-gamma)*exp(-days_from_peak^2*beta))

### calculate area under sensitivity curve before and after cut-offs defined in Brummer
dat_early = filter(df_sens, day <= incubation_period+6)
dat_late = filter(df_sens, day >= incubation_period+7)
sens_rel_early = sum(dat_early$sens_adj_rel)/sum(dat_early$sens)
sens_rel_late = sum(dat_late$sens_adj_rel)/sum(dat_late$sens)

# put data in long format
df_sens_long = df_sens%>%
  dplyr::select(-c(fnr, Source))%>%
  pivot_longer(-c(day, days_from_peak), names_to = 'Test', values_to = 'sensitivity')


###############################
### Plot sensitivity curves ###
###############################

# adjust factors for plotting
df_sens_long_plot = df_sens_long

df_sens_long_plot$Test = factor(df_sens_long_plot$Test, levels = c('sens', 'sens_adj_const', 'sens_adj_rel'),
                                labels = c('RT-PCR (Kucirka et al.)', 'Ag-RDT (constant reduction)', 'Ag-RDT (relative reduction)'))


p_sens_curves = ggplot(df_sens_long_plot, aes(x = day, y = sensitivity*100, ymin = 0, ymax = sensitivity, colour = Test, fill = Test))+
  geom_point()+
  geom_line()+
  theme_bw()+
  ylab("diagnostic sensitivity (%)")+
  xlab("days since SARS-CoV-2 exposure")+
  geom_vline(xintercept = 11.5, size = 2, alpha = 0.5)+
  geom_text(x = 7, y = 5, label = paste0("cumulative sensitivity\n<7 days before\nsymptom onset: \n ", round(sens_rel_early, 4)*100, "%"), 
            colour = 'black', size = 3)+
  geom_text(x = 16, y = 5, label = paste0("cumulative sensitivity\n>=7 days before\nsymptom onset: \n ",round(sens_rel_late,4)*100, "%"), 
            colour = 'black', size = 3)
p_sens_curves



###################################
### Save final sensitivity data ###
###################################

# ### Main analysis: 
# # RT-PCR: extrapolated sensitivity from Kurcika et al.
# # Ag-RDT A: uniform reduction in sensitivity relative to RT-PCR, using estimate from Brummer et al.
# # Ag-RDT B: time-varying reduction in sensitivity relative to RT-PCR, using exponential function with shape parameters defined above

# write.csv(df_sens_long, paste0(filepath_sensitivity_curves, "kucirka_adjusted.csv"))

### Sensitivty analysis 1
# a uniform sensitivity function
# for a sensitivity analysis to see to what degree time-varying sensitivity drives stochasticity in our results
sens_pcr_uniform = 0.7
sens_rdt_uniform = sens_pcr_uniform*0.738
df_sens_long_uniform = df_sens_long%>%
  mutate(sensitivity = case_when(Test == "sens" ~ sens_pcr_uniform,
                                 Test != "sens" ~ sens_rdt_uniform))
#write.csv(df_sens_long_uniform, paste0(filepath_sensitivity_curves, "kucirka_uniform.csv"))

### Sensitivity analsis 2
# a perfect sensitivity function
sens_pcr_perfect = 1
sens_rdt_perfect = 1
df_sens_long_perfect = df_sens_long%>%
  mutate(sensitivity = case_when(Test == "sens" ~ sens_pcr_perfect,
                                 Test != "sens" ~ sens_rdt_perfect))
#write.csv(df_sens_long_perfect, paste0(filepath_sensitivity_curves, "kucirka_perfect.csv"))
