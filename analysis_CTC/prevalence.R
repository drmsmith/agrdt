source("functions.R")

##################
### PREVALENCE ###
##################
lot = "lot3"

### define vectors for loop as column names for prevalence output matrix (counts of different infection statuses)
vec_ltcf = 1:3
vec_SIM = 1:100
colnames_prev = c('countAbs', 'count0', 'count1', 'count2I', 'count2II', 'count3I', 'count3II', 'count3III', 'count4', 'day',
                  'immunity', 'IPC', 'network','simCTC')
n_days = 15

### empty matrix
m_prevalence = matrix(ncol = length(colnames_prev), nrow = length(vec_ltcf)*length(vec_SIM)*n_days)
colnames(m_prevalence) = colnames_prev


### LOOP
# Use a loop to go through all files and calculate prevalence using:
# Daily status data (prevalence)

qounter = 1
for(LTCF_x in vec_ltcf){
  if(LTCF_x == 1){immunity = "0.20.2"; IPC = 1; Network = 1}
  if(LTCF_x == 2){immunity = "0.20.2"; IPC = 1; Network = 2}
  if(LTCF_x == 3){immunity = "0.50.5"; IPC = 2; Network = 2}
  print(paste0("Immunity ", immunity, ", IPC ", IPC, " and Network", Network, "."))
  
  for(SIM in vec_SIM){
    
    dataPrevSummarized_i = as.matrix(f_prevalence_summarize(lot, immunity, IPC, Network, SIM))
    
    m_prevalence[qounter:(qounter+n_days-1),] = dataPrevSummarized_i
    
    qounter = qounter + n_days
  }
}



### make matrix into data.frame with numeric variables (with exception of immunity)
df_prevalence = as.data.frame(m_prevalence)
df_prevalence[,c(1:10, 12:14)] <- sapply(df_prevalence[,c(1:10, 12:14)],as.character)
df_prevalence[,c(1:10, 12:14)] <- sapply(df_prevalence[,c(1:10, 12:14)],as.numeric)


### simplifications: pivot data for plotitng, combine the two presymptomatics when pivoting the data, and also combine mild and severe symptoms
df_prevalenceLong = df_prevalence%>%
  mutate(count2 = count2I + count2II)%>%
  mutate(count3IV = count3II + count3III)%>%
  dplyr::select(-c(count2I, count2II, count3II, count3III))%>%
  pivot_longer(-c(day, immunity, IPC, network, simCTC), names_to = 'status', values_to = 'count')%>%
  mutate(networkIPC = paste0("Network", network, "_IPC", IPC))%>%
  mutate(LTCF = paste0(immunity, "_", networkIPC))

### calculate mean and median incidence across all 100 CTC simulations
df_prevalenceLongSummarized = df_prevalenceLong%>%group_by(day, immunity, networkIPC, status)%>%
  dplyr::summarize(median = median(count),
                   mean = mean(count))

### data for 3 selected LTCFs
vec_3ltcfs_labels = c('LTCF 1 (high risk)', 'LTCF 2 (moderate risk)', 'LTCF 3 (low risk)')

df_prevalenceLong_3ltcfs = df_prevalenceLong%>%
  mutate(LTCF = ifelse(LTCF == '0.20.2_Network1_IPC1', vec_3ltcfs_labels[1],
                       ifelse(LTCF == '0.20.2_Network2_IPC1', vec_3ltcfs_labels[2], vec_3ltcfs_labels[3])))

df_prevalenceLongSummarized_3ltcfs = df_prevalenceLong_3ltcfs%>%group_by(day, LTCF, status)%>%
  dplyr::summarize(median = median(count),
                   mean = mean(count))


### select and factor status
vec_statuses = c('count1', 'count2', 'count3I', 'count3IV')
vec_status_labels = c('exposed', 'pre-symptomatic', 'asymptomatic', 'symptomatic (mild or severe)')

# Selected LTCFs
df_prevalenceLongSummarizedInfected_3ltcfs = df_prevalenceLongSummarized_3ltcfs%>%
  filter(status %in% vec_statuses)

df_prevalenceLongSummarizedInfected_3ltcfs$status = factor(df_prevalenceLongSummarizedInfected_3ltcfs$status, levels = rev(vec_statuses), labels = rev(vec_status_labels))
