source("functions.R")


###########################
### CALCULATE INCIDENCE ###
###########################
lot = "lot3"

### initialize loop
vec_ltcf = 1:3
vec_SIM = 1:100
m_incidence_nosocomial = matrix(ncol = 13, nrow = length(vec_ltcf)*length(vec_SIM))

### Run loop
qounter = 0
for(LTCF_x in vec_ltcf){
  if(LTCF_x == 1){immunity_x = "0.20.2"; IPC_x = 1; Network_x = 1}
  if(LTCF_x == 2){immunity_x = "0.20.2"; IPC_x = 1; Network_x = 2}
  if(LTCF_x == 3){immunity_x = "0.50.5"; IPC_x = 2; Network_x = 2}

  print(paste0("LTCF ", LTCF_x))
  for(SIM in vec_SIM){

    qounter = qounter + 1
    fileroot = f_fileroot(lot, immunity_x,IPC_x,SIM,Network_x)
    dataSecondaryCases = f_read_csv(filepath_CTC, fileroot, "secondaryCases")
    dataIndex = f_read_csv(filepath_CTC, fileroot, "FirstIndex")
    dataCommCases = f_read_csv(filepath_CTC, fileroot, "CommCases")

    incNosocomial = f_incNosocomial(dataSecondaryCases)
    incIndex = f_incIndex(dataIndex)
    incCommCases = f_incCommCases(dataCommCases)

    incNosocomialPA = sum(filter(incNosocomial, cat == 'Patient')$Nb)
    incNosocomialPE = sum(filter(incNosocomial, cat == 'Staff')$Nb)
    incNosocomialNb = sum(incNosocomial$Nb)
    if(incNosocomialNb != incNosocomialPA + incNosocomialPE){warning("incNosocomial: patients and staff don't add up")}

    incIndexPA = sum(filter(incIndex, cat == 'Patient')$Nb)
    incIndexPE = sum(filter(incIndex, cat == 'Staff')$Nb)
    incIndexNb = sum(incIndex$Nb)
    if(incIndexNb != incIndexPA + incIndexPE){warning("incIndex: patients and staff don't add up")}

    incCommCasesPA = sum(filter(incCommCases, cat == 'Patient')$Nb)
    incCommCasesPE = sum(filter(incCommCases, cat == 'Staff')$Nb)
    incCommCasesNb = sum(incCommCases$Nb)
    if(incCommCasesNb != incCommCasesPA + incCommCasesPE){warning("incCommCases: patients and staff don't add up")}

    incVector = c(incNosocomialPA, incNosocomialPE, incNosocomialNb,
                  incIndexPA, incIndexPE, incIndexNb,
                  incCommCasesPA, incCommCasesPE, incCommCasesNb)

    m_incidence_nosocomial[qounter,] = c(immunity_x, IPC_x, SIM, Network_x, incVector)
  }
}

# re-orient data
df_incidence_nosocomial = m_incidence_nosocomial%>%
  as.data.frame()
colnames(df_incidence_nosocomial) = c('immunity', 'IPC', 'SIM', 'Network',
                                      'incNosocomialPA', 'incNosocomialPE', 'incNosocomialTotal',
                                      'incIndexPA', 'incIndexPE', 'incIndexTotal',
                                      'incCommCasesPA', 'incCommCasesPE', 'incCommCasesTotal')

# make incidence columns numeric
df_incidence_nosocomial[,5:13] = sapply(sapply(df_incidence_nosocomial[,5:13], as.character),as.numeric)

# add IPC_Network column to string them together
df_incidence_nosocomial = df_incidence_nosocomial%>%
  mutate(NetworkIPC = paste0("Network",Network,"_IPC",IPC))%>%
  mutate(LTCF = paste0(immunity, "_", NetworkIPC))

# Label LTCFs (NB: "risk" labels later replaced by "control" labels)
vec_3ltcfs_labels = c('LTCF 1 (high risk)', 'LTCF 2 (moderate risk)', 'LTCF 3 (low risk)')

df_incidence_nosocomial_3ltcfs = df_incidence_nosocomial%>%
  mutate(LTCF = ifelse(LTCF == '0.20.2_Network1_IPC1', vec_3ltcfs_labels[1],
                            ifelse(LTCF == '0.20.2_Network2_IPC1', vec_3ltcfs_labels[2], vec_3ltcfs_labels[3])))


######################
### SUPERSPREADING ###
######################

vec_SIM = 1:100
vec_ltcf = 1:3

vec_thresholdSuperspreading = 3:5

m_incidence_superspreading = matrix(ncol = 17, nrow = length(vec_ltcf)*length(vec_SIM)*length(vec_thresholdSuperspreading))

### Run loop

qounter = 0
for(LTCF_x in vec_ltcf){
  if(LTCF_x == 1){immunity = "0.20.2"; IPC = 1; Network = 1}
  if(LTCF_x == 2){immunity = "0.20.2"; IPC = 1; Network = 2}
  if(LTCF_x == 3){immunity = "0.50.5"; IPC = 2; Network = 2}
  for(SIM in vec_SIM){
    for(thresholdSuperspreading in vec_thresholdSuperspreading){
      
      qounter = qounter + 1
      fileroot = f_fileroot(lot, immunity,IPC,SIM,Network)
      dataSecondaryCases = f_read_csv(filepath_CTC, fileroot, "secondaryCases")
      dataIndex = f_read_csv(filepath_CTC, fileroot, "FirstIndex")
      dataCommCases = f_read_csv(filepath_CTC, fileroot, "CommCases")
      
      # ids of nosocomial cases
      idsNosocomial = f_idsNosocomial(dataSecondaryCases)
      
      if(nrow(idsNosocomial) == 0){m_incidence_superspreading[qounter,] = c(immunity, IPC, SIM, Network, thresholdSuperspreading, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA); next()}
      
      # neverspreading
      dataNeverspreaders = f_neverSpreaders(dataIndex, dataCommCases, dataSecondaryCases)

      # superspreading
      dataSuperspreaders = f_superSpreaders(idsNosocomial)
      shareOfSuperspreading = f_shareOfSuperspreading(dataNeverspreaders, dataSuperspreaders, threshold = thresholdSuperspreading)

      vec_Neverspreading = filter(shareOfSuperspreading,Superspreader == 0)%>%dplyr::select(-Superspreader)
      vec_Spreading = filter(shareOfSuperspreading,Superspreader == 1)%>%dplyr::select(-Superspreader)
      vec_Superspreading = filter(shareOfSuperspreading,Superspreader == 2)%>%dplyr::select(-Superspreader)
      
      if(nrow(vec_Neverspreading) == 0){vec_Neverspreading = c(CountSpreaders_Never = 0, CountInfectees_Never = 0, FreqSpreaders_Never = 0, FreqInfectees_Never = 0)}
      if(nrow(vec_Spreading) == 0){vec_Spreading = c(CountSpreaders_Spreader = 0, CountInfectees_Spreader = 0, FreqSpreaders_Spreader = 0, FreqInfectees_Spreader = 0)}
      if(nrow(vec_Superspreading) == 0){vec_Superspreading = c(CountSpreaders_Super = 0, CountInfectees_Super = 0, FreqSpreaders_Super = 0, FreqInfectees_Super = 0)}
      
      m_incidence_superspreading[qounter,] = c(immunity, IPC, SIM, Network, thresholdSuperspreading, as.character(c(vec_Neverspreading, vec_Spreading, vec_Superspreading)))
      
    }
  }
}

# re-orient data
df_incidence_superspreading = m_incidence_superspreading%>%
  as.data.frame()

colnames(df_incidence_superspreading) = c('immunity', 'IPC', 'SIM', 'Network', 'thresholdSuperspreading', 
                                          'countSpreadersNever', 'countInfecteesNever', 'freqSpreadersNever', 'freqInfecteesNever',
                                          'countSpreadersSpreader', 'countInfecteesSpreader', 'freqSpreadersSpreader', 'freqInfecteesSpreader',
                                          'countSpreadersSuper', 'countInfecteesSuper', 'freqSpreadersSuper', 'freqInfecteesSuper')

df_incidence_superspreading = filter(df_incidence_superspreading, countSpreadersSuper != "<NA>")

# make incidence columns numeric
df_incidence_superspreading[,6:17] = sapply(sapply(df_incidence_superspreading[,6:17], as.character),as.numeric)

df_incidence_superspreading = df_incidence_superspreading%>%
  mutate(NetworkIPC = paste0("Network",Network,"_IPC",IPC))%>%
  mutate(LTCF = paste0(immunity, "_", NetworkIPC))%>%
  pivot_longer(-c(immunity, IPC, SIM, Network, thresholdSuperspreading, NetworkIPC, LTCF), names_to = "measure", values_to = "outcome")%>%
  mutate(direction = ifelse(grepl("Infectee", measure), "infectee", "infector"),
         magnitude = ifelse(grepl("Never", measure), "never", ifelse(grepl("Super", measure), "super", "spreader")),
         measure = ifelse(grepl("count", measure), "count", "frequency"))


df_incidence_superspreading_3ltcfs = df_incidence_superspreading%>%
  mutate(LTCF = ifelse(LTCF == '0.20.2_Network1_IPC1', vec_3ltcfs_labels[1],
                            ifelse(LTCF == '0.20.2_Network2_IPC1', vec_3ltcfs_labels[2], vec_3ltcfs_labels[3])))

######################
### SUPERSPREADING ###
######################
vec_superspreading_labels = c('low spreader', 'super spreader', 'non spreader')

### plot with both frequencies side-by-side
df_incidence_superspreading_frequencies = df_incidence_superspreading%>%
  filter(measure == 'frequency')
### and specific to 3 LTCFs
df_incidence_superspreading_frequencies_3ltcfs = df_incidence_superspreading_3ltcfs%>%
  filter(measure == 'frequency')


### Fix up factor labels, remove infectee/never combo (always zero) etc.

df_incidence_super_neverspreading_frequencies_3ltcfs = df_incidence_superspreading_frequencies_3ltcfs%>%
  filter(!c(direction == 'infectee' & magnitude == 'never'))%>%
  dplyr::select(-c(immunity, IPC, Network, thresholdSuperspreading, NetworkIPC))

df_incidence_super_neverspreading_frequencies_3ltcfs$magnitude = factor(df_incidence_super_neverspreading_frequencies_3ltcfs$magnitude,
                                                                        levels = c('spreader', 'super', 'never'),
                                                                        labels = vec_superspreading_labels)


### DENSITY: MAKE VIOLINS THE SAME SIZE

# following instructions here: https://stackoverflow.com/questions/47174825/same-area-for-all-violins-independent-of-facets-in-ggplot2
# determine density
# limit density (using from/to) to min and max of frequency

df_incidence_super_neverspreading_frequencies_3ltcfs_DENSITY_infectors = df_incidence_super_neverspreading_frequencies_3ltcfs%>%
  filter(measure == 'frequency', direction == 'infector')%>%
  group_by(magnitude, LTCF) %>%
  do({
    dens <- density(.$outcome, from = min(.$outcome), to = max(.$outcome))
    tibble(x = c(head(dens$x, 1), dens$x, tail(dens$x, 1)), #Add 0s at end to close lines
           y = c(0, dens$y, 0))
  }) %>% 
  ungroup()%>% 
  mutate(y_adj = ifelse(LTCF == vec_3ltcfs_labels[3] & magnitude == vec_superspreading_labels[2], y*0.25, y))%>%
  mutate(ymin = as.numeric(magnitude) - y_adj/max(y_adj), # Add offset for factor levels
         ymax = as.numeric(magnitude) + y_adj/max(y_adj))
                                                              
df_incidence_super_neverspreading_frequencies_3ltcfs_DENSITY_infectees = df_incidence_super_neverspreading_frequencies_3ltcfs%>%
  filter(measure == 'frequency', direction == 'infectee')%>%
  group_by(magnitude, LTCF) %>%
  do({
    dens <- density(.$outcome, from = min(.$outcome), to = max(.$outcome))
    tibble(x = c(head(dens$x, 1), dens$x, tail(dens$x, 1)), #Add 0s at end to close lines
           y = c(0, dens$y, 0))
  }) %>% 
  ungroup()%>% 
  mutate(ymin = as.numeric(magnitude) - y/max(y), # Add offset for factor levels
         ymax = as.numeric(magnitude) + y/max(y))

###########################
### INCIDENCE OVER TIME ###
###########################

vec_ltcf = 1:3
vec_SIM = 1:100

# empty df to save data
incAll_total = data.frame()

for(LTCF_x in vec_ltcf){
  if(LTCF_x == 1){immunity = "0.20.2"; IPC = 1; Network = 1}
  if(LTCF_x == 2){immunity = "0.20.2"; IPC = 1; Network = 2}
  if(LTCF_x == 3){immunity = "0.50.5"; IPC = 2; Network = 2}
  for(SIM in vec_SIM){
    
    fileroot = f_fileroot(lot, immunity,IPC,SIM,Network)
    dataSecondaryCases = f_read_csv(filepath_CTC, fileroot, "secondaryCases")
    dataCommCases = f_read_csv(filepath_CTC, fileroot, "CommCases")
    dataIndex = f_read_csv(filepath_CTC, fileroot, "FirstIndex")
    
    incNosocomial = f_incNosocomial(dataSecondaryCases)
    incCommCases = f_incCommCases(dataCommCases)
    incIndex = f_incIndex(dataIndex)%>%
      mutate(Measure = 'community-onset cases')
    
    incAll = rbind(incNosocomial,
                   incCommCases,
                   incIndex)%>%
      mutate(immunity = immunity,
             IPC = IPC,
             Network = Network,
             SIM = SIM)
    
    incAll_Total = incAll%>%
      group_by(Date, immunity, IPC, Network, SIM, Measure)%>%
      dplyr::summarise(Nb = sum(Nb))%>%
      mutate(cat = 'Total')%>%
      as.data.frame()
    
    incAll_total = rbind(incAll_total, rbind(incAll, incAll_Total))
  }
}


### Complete missing data
# i.e. put in zeroes for days where there are no cases
dfIncAll = incAll_total%>%
  complete(Date, nesting(cat, immunity, IPC, Network, SIM, Measure), fill = list(Nb = 0))


### Combine with data for zero incidence on day zero
vec_Measure = c('community-onset cases', 'secondary cases')

dfIncAll_total = dfIncAll

### Combine 3 LTCFs
dfIncLTCF1 = dfIncAll_total%>%filter(IPC == 1, Network == 1, immunity == "0.20.2")%>%
  mutate(LTCF = vec_3ltcfs_labels[1])

dfIncLTCF2 = dfIncAll_total%>%filter(IPC == 1, Network == 2, immunity == "0.20.2")%>%
  mutate(LTCF = vec_3ltcfs_labels[2])

dfIncLTCF3 = dfIncAll_total%>%filter(IPC == 2, Network == 2, immunity == "0.50.5")%>%
  mutate(LTCF = vec_3ltcfs_labels[3])

dfInc_3ltcfs = rbind(dfIncLTCF1, dfIncLTCF2, dfIncLTCF3)%>%
  mutate(day = as.numeric(abs(as.Date('2009-07-05') - Date)))

dfInc_3ltcfs_summarized = dfInc_3ltcfs%>%
  group_by(Date, cat, LTCF, Measure)%>%
  summarise(Nb_median = median(Nb),
            Nb_mean = mean(Nb),
            Nb_025 = quantile(Nb, 0.025),
            Nb_250 = quantile(Nb, 0.25),
            Nb_750 = quantile(Nb, 0.75),
            Nb_975 = quantile(Nb, 0.975))%>%
  mutate(day = as.numeric(abs(as.Date('2009-07-05') - Date)))


### factor patient/staff category
dfInc_3ltcfs$cat = factor(dfInc_3ltcfs$cat, levels = c('Staff', 'Patient', 'Total'),
                                     labels = c('staff', 'patients', 'patients & staff'))
dfInc_3ltcfs$Measure = factor(dfInc_3ltcfs$Measure, levels = c('community-onset cases', 'secondary cases'),
                              labels = c('community', 'hospital'))

dfInc_3ltcfs_summarized$cat = factor(dfInc_3ltcfs_summarized$cat, levels = c('Staff', 'Patient', 'Total'),
                                     labels = c('staff', 'patients', 'patients & staff'))
dfInc_3ltcfs_summarized$Measure = factor(dfInc_3ltcfs_summarized$Measure, levels = c('community-onset cases', 'secondary cases'),
                                     labels = c('community', 'hospital'))

### incidence by type of individual (combined community- and hospital-onset)
dfInc_3ltcfs_cat = filter(dfInc_3ltcfs, cat %in% c('patients', 'staff'))%>%
  group_by(Date, immunity, IPC, Network, SIM, cat, LTCF, day)%>%
  summarise(Nb= sum(Nb))
dfInc_3ltcfs_cat_summarized = dfInc_3ltcfs_cat%>%
  group_by(Date, immunity, IPC, Network, cat, LTCF, day)%>%
  summarise(Nb_mean = mean(Nb),
            Nb_median = median(Nb),
            Nb_975 = quantile(Nb, 0.975),
            Nb_025 = quantile(Nb, 0.025))
