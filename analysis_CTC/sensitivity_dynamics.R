source("functions.R")


######################################################
### OVERALL TIME-VARYING NATURE OF PCR SENSITIVITY ###
######################################################

### For each CTC simulation in each scenario, determine mean time-varying PCR sensitivity among infected individuals

df_sensitivityMeans_timevarying = data.frame()
df_sensitivityMeans_uniform = data.frame()

for(LTCF_i in 1:3){
  for(simulationCTC_i in 1:100){
    print(paste0('simulation ', simulationCTC_i, ' of 100'))

    ### Set relevant fileroot
    if(LTCF_i == 1){fileroot_sens = f_fileroot(lot, "0.20.2",1,simulationCTC_i,1)}
    if(LTCF_i == 2){fileroot_sens = f_fileroot(lot, "0.20.2",1,simulationCTC_i,2)}
    if(LTCF_i == 3){fileroot_sens = f_fileroot(lot, "0.50.5",2,simulationCTC_i,2)}

    ### Infection status (and vector of who not infected on any given day)
    dataStatus_sens = f_read_csv(filepath_CTC, fileroot_sens, "statusByDay")
    whichNotInfected_sens = which(!as.matrix(dataStatus_sens[,4:ncol(dataStatus_sens)]) %in% c('1', '2I', '2II', '3I', '3II', '3III'))

    ### Introductions
    dataCommCases_sens = f_read_csv(filepath_CTC, fileroot_sens, 'CommCases')

    ### Index cases
    dataIndex_sens = f_read_csv(filepath_CTC, fileroot_sens, "FirstIndex")

    ### Hospital-onset cases
    dataSecondaryCases_sens = f_read_csv(filepath_CTC, fileroot_sens, "secondaryCases")

    ### Determine which individuals are hospital-onset cases vs. community-onset (for the latter, combine index cases and introductions)
    casesNosocomial_sens = dataSecondaryCases_sens$Cases
    casesCommunity_sens = c(as.character(dataCommCases_sens$Individual), as.character(dataIndex_sens$Individual))


    ### Sensitivty curves for each (time-varying and uniform)
    dataSensitivityPCRtimevarying_sens = f_dataSensitivity(dataStatus_sens, "PCR", sensitivity_function = "time-varying")
    dataSensitivityPCRtimevarying_sens[,4:ncol(dataSensitivityPCRtimevarying_sens)][whichNotInfected_sens] = NA

    dataSensitivityPCRuniform_sens = f_dataSensitivity(dataStatus_sens, "PCR", sensitivity_function = "uniform")
    dataSensitivityPCRuniform_sens[,4:ncol(dataSensitivityPCRuniform_sens)][whichNotInfected_sens] = NA

    ### Build time-varying curves
    df_SensitivityPCRtimevarying_sens = dataSensitivityPCRtimevarying_sens%>%
      as.data.frame()%>%
      mutate(onset = ifelse(id %in% casesNosocomial_sens, 'hospital-onset',
                            ifelse(id %in% casesCommunity_sens, 'community-onset',
                                   'none')))%>%
      pivot_longer(-c(id, cat, ward, onset), names_to = "date", values_to = "sensitivity")%>%
      mutate(sensitivity = as.numeric(as.character(sensitivity)))%>%
      group_by(onset, date)%>%
      summarise(sensitivity_mean = mean(sensitivity, na.rm = T))%>%
      mutate(LTCF = vec_3ltcfs_labels[LTCF_i],
             simulationCTC = simulationCTC_i)%>%
      as.data.frame()

    df_sensitivityMeans_timevarying = rbind(df_sensitivityMeans_timevarying, df_SensitivityPCRtimevarying_sens)

    ### Build uniform curves
    df_SensitivityPCRuniform_sens = dataSensitivityPCRuniform_sens%>%
      as.data.frame()%>%
      mutate(onset = ifelse(id %in% casesNosocomial_sens, 'hospital-onset',
                            ifelse(id %in% casesCommunity_sens, 'community-onset',
                                   'none')))%>%
      pivot_longer(-c(id, cat, ward, onset), names_to = "date", values_to = "sensitivity")%>%
      mutate(sensitivity = as.numeric(as.character(sensitivity)))%>%
      group_by(onset, date)%>%
      summarise(sensitivity_mean = mean(sensitivity, na.rm = T))%>%
      mutate(LTCF = vec_3ltcfs_labels[LTCF_i],
             simulationCTC = simulationCTC_i)%>%
      as.data.frame()

    df_sensitivityMeans_uniform = rbind(df_sensitivityMeans_uniform, df_SensitivityPCRuniform_sens)
  }
}

#### ORAGNIZE DATA

# time-varying
df_sensitivityMeans_timevarying = df_sensitivityMeans_timevarying%>%
  mutate(time = abs(as.numeric(as.Date("2009-07-05") - as.Date(date))))%>%
  filter(onset != 'none')

df_sensitivityMeans_timevarying_Means = df_sensitivityMeans_timevarying%>%
  group_by(time, LTCF, onset)%>%
  dplyr::summarise(sensitivity_mean_means = mean(sensitivity_mean, na.rm = T))

# uniform
df_sensitivityMeans_uniform = df_sensitivityMeans_uniform%>%
  mutate(time = abs(as.numeric(as.Date("2009-07-05") - as.Date(date))))%>%
  filter(onset != 'none')

df_sensitivityMeans_uniform_Means = df_sensitivityMeans_uniform%>%
  group_by(time, LTCF, onset)%>%
  dplyr::summarise(sensitivity_mean_means = mean(sensitivity_mean, na.rm = T))
