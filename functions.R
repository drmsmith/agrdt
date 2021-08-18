library(tidyr)
library(magrittr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(simpleboot)

# filepaths
filepath = paste0(getwd(),"/")
filepath_noel = filepath
filepath_CTC= paste0(filepath, "Model_herdI_incidenceLow/output/")
filepath_surveillance = paste0(filepath, "Surveillance/output/")

colnames_CTC = c("id", "cat", "ward", paste0("2009-07-0",6:9), paste0("2009-07-",10:20), "Last")

### FILEROOTS
### Function to generate stable fileroot usable on all file types (statusByDay, secondaryCases, etc)
f_fileroot = function(lot, immunity, IPC, SIM, Network){
  if(!immunity %in% c("0.20.2", "0.50.5", "0.80.2")){warning("Wrong immunity");break()}
  if(!IPC %in% 1:2){warning("Wrong IPC");break()}
  if(!SIM %in% 1:100){warning("Wrong SIM");break()}
  if(!Network %in% 1:2){warning("Wrong Network");break()}
  if(!lot %in% c("lot1", "lot2", "lot3")){warning("Wrong lot");break()}
  
  ### LOTS 1 AND 2 "high incidence" scenarios
  if(lot %in% c("lot1", "lot2")){
    fileroot = paste0("_",immunity,"Christmasnetwork", Network,"IPC",IPC,
                      "_SIM_",SIM,"_Bacteria_Coronaviridae_Alphaletovirus_Covid19___.csv")
  }
  
  
  ### LOT 3: baseline "low incidence" scenario; these filepaths are saved with "low_" before "SIM
  if(lot %in% c("lot3")){
    fileroot = paste0("_",immunity,"Christmasnetwork", Network,"IPC",IPC,
                      "low_SIM_",SIM,"_Bacteria_Coronaviridae_Alphaletovirus_Covid19___.csv")
  }
  return(fileroot)
}


### READ CSV
### Functions to read corresponding CSVs depending on filetype
### with some minor processing: 
### (1) fix colnames for CTC data 
### (2) split cases from secondaryCases data connected by "_"

### NB. MANUALLY UPDATED "low" TO REFLECT LOW INCIDENCE SIMULATIONS FOR LOT3

f_read_csv = function(filepath, fileroot, filetype){
  if(!filetype %in% c("FirstIndex","statusByDay", "secondaryCases", "CommCases")){warning("Wrong data type");break()}
  
  # NB: typo in CTC data so filetypes inconsistent; so adjust accordingly:
  if(filetype == "secondaryCases"){filetype1 = "secondaryCases"; filetype2 = "secondayCases"}else{
    filetype1 = filetype; filetype2 = filetype
  }
  
  # CommCases data does not have underscore before fileroot, unlike others, so must remove first character from fileroot string
  if(filetype == "CommCases"){fileroot = str_sub(fileroot, 2)}
  
  # read data
  dat = read.csv2(paste0(filepath, filetype1,
                         #"/Bacteria_Coronaviridae_Alphaletovirus_Covid19___/christmas/",
                         "/christmas/",
                         filetype2,fileroot))
  
  # Fix statusByDay column names
  if(filetype == "statusByDay"){colnames(dat) = colnames_CTC; dat = dat%>%dplyr::select(-Last)}
  
  # Fix secondaryCases: make full line list (instead of cases strung by "_")
  if(filetype == "secondaryCases"){dat = separate_rows(dat, Cases, length, Pathway, sep = "_")%>%
    dplyr::select(-c(Nb,length))}
  
  return(dat)
}

# ## TEST
# fileroot = f_fileroot("lot3", "0.20.2",1,1,1)
# dataStatus = f_read_csv(filepath_CTC, fileroot, "statusByDay")
# dataSecondaryCases = f_read_csv(filepath_CTC, fileroot, "secondaryCases")
# dataIndex = f_read_csv(filepath_CTC, fileroot, "FirstIndex")
# dataCommCases = f_read_csv(filepath_CTC, fileroot, 'CommCases')



#################
### INCIDENCE ###
#################
# Functions that calculate incidence from CTC output data

### INCIDENCE NOSOCOMIAL
### Function to return total number of nosocomial infections across different scenarios
f_incNosocomial = function(dataSecondaryCases){
  incNosocomial = dataSecondaryCases%>%
    mutate(cat = ifelse(str_detect(Cases, "PA-"), 'Patient', 'Staff'))%>%
    group_by(Date, cat)%>%summarise(Nb = n())%>%
    mutate(Measure = "secondary cases")%>%
    as.data.frame()%>%
    mutate(Date = as.Date(Date))
  
  return(incNosocomial)
}
# ## TEST
# f_incNosocomial(dataSecondaryCases)

### INCIDENCE INDEX
### Function to return index cases from day 1
f_incIndex = function(dataIndex){
  incIndex = dataIndex%>%
    mutate(cat = ifelse(str_detect(Individual, "PA-"), 'Patient', 'Staff'))%>%
    group_by(Date, cat)%>%summarise(Nb = n())%>%
    mutate(Measure = "index cases")%>%
    as.data.frame()%>%
    mutate(Date = as.Date(Date))
  return(incIndex)
}
# ## TEST
# f_incIndex(dataIndex)

### INCIDENCE COMM-ONSET
### Function to return community-onset cases other than the initial 4 index cases
f_incCommCases = function(dataCommCases){
  incCommCases = dataCommCases%>%
    mutate(cat = ifelse(str_detect(Individual, "PA-"), 'Patient', 'Staff'))%>%
    group_by(Date, cat)%>%summarise(Nb = n())%>%
    mutate(Measure = "community-onset cases")%>%
    as.data.frame()%>%
    mutate(Date = as.Date(Date))
  return(incCommCases)
}
# ## TEST
# f_incCommCases(dataCommCases)

#####################
### IDS FOR TREES ###
#####################
### f_ids: determining who is infected and when, for pruning of transmission chains

f_idsNosocomial = function(dataSecondaryCases){
  idsNosocomial = dataSecondaryCases%>%
    mutate(cat = ifelse(str_detect(Cases, "PA-"), 'Patient', 'Staff'))%>%
    as.data.frame()%>%
    mutate(Date = as.Date(Date))%>%
    mutate(Infector = Individual)%>%
    mutate(Infectee = Cases)%>%
    dplyr::select(-Cases, -Individual)
  
  return(idsNosocomial)
}
# ## TEST
# idsNosocomial = f_idsNosocomial(dataSecondaryCases)


######################
### SUPERSPREADERS ###
######################

# NEVERSPREADERS - function to identify who is infected but never spreads
f_neverSpreaders = function(dataIndex, dataCommCases, dataSecondaryCases){
  # vector of IDs of all infected individuals: index cases, comm-onset and hosp-onset
  idsInfectees = c(as.character(dataIndex$Individual), as.character(dataCommCases$Individual), as.character(dataSecondaryCases$Cases))
  
  # vector of IDs of all infector indivdiuals
  idsInfectors = unique(as.character(dataSecondaryCases$Individual))
  
  # which of idsInfectees are in infectors
  idsNeverspreaders = unique(idsInfectees[which(!idsInfectees %in% idsInfectors)])
  
  ### NB: Because of transient carriage (HCW vectors), IDs don't always add up: can have more infectors than infectees, and this is by design
  
  dataNeverspreaders = data.frame(numInfectees = length(idsInfectees),
                                  numInfectors = length(idsInfectors),
                                  numNeverspreaders = length(idsNeverspreaders))%>%
    mutate(shareOfNeverspreading = numNeverspreaders/(numInfectors + numNeverspreaders))
                                  
  return(dataNeverspreaders)
}
# ## TEST
# dataNeverspreaders = f_neverSpreaders(dataIndex, dataCommCases, dataSecondaryCases)


# Function to calculate the frequency of the number of infectees per infection
f_superSpreaders = function(idsNosocomial){
  
  dataSuperspreaders = idsNosocomial%>%
    group_by(Infector)%>%
    summarise(NbInfectees = n())%>%
    group_by(NbInfectees)%>%
    summarise(Count = n())
  
  return(dataSuperspreaders)
}
# ## TEST
# dataSuperspreaders = f_superSpreaders(idsNosocomial)

### Function to calculate the share of infectors that are superspreaders, and the share of transmissions they account for
# In the article, "3" is taken as the superspreading threshold 
f_shareOfSuperspreading = function(dataNeverspreaders, dataSuperspreaders, threshold = 3){
  if(!threshold %in% 2:5){warning('Superspreader threshold is not 2, 3, 4 or 5')}
  
  dataSuperspreaders_with_Neverspreaders = data.frame(NbInfectees = 0, Count = dataNeverspreaders$numNeverspreaders)%>%
    rbind(., dataSuperspreaders)
  
  dataShareOfSuperspreading = dataSuperspreaders_with_Neverspreaders%>%
    mutate(Superspreader = ifelse(NbInfectees > threshold, 2, ifelse(NbInfectees > 0, 1, 0)))%>%
    group_by(Superspreader)%>%summarise(CountSpreaders = sum(Count),
                                        CountInfectees = sum(NbInfectees*Count))%>%
    mutate(FreqSpreaders = CountSpreaders/sum(CountSpreaders),
           FreqInfectees = CountInfectees/sum(CountInfectees))
  
  return(dataShareOfSuperspreading)
}

# ## TEST
# f_shareOfSuperspreading(dataNeverspreaders, dataSuperspreaders, 3)


##################
### PREVALENCE ###
##################

### Function to summarize prevalence data (dataStatus) 
f_prevalence_summarize = function(lot, immunity, IPC, Network, SIM){
  
  fileroot = f_fileroot(lot, immunity,IPC,SIM,Network)
  dataStatus = f_read_csv(filepath_CTC, fileroot, "statusByDay")
  
  dataPrev = dataStatus[,4:ncol(dataStatus)]%>%sapply(.,summary)
  
  # loop through each day and count prevalence of each class of infection statuses
  df_dataPrevSummarized = data.frame()
  for(date in 1:length(dataPrev)){
    
    dataPrev_i = do.call(rbind, dataPrev[date])%>%t()
    
    if('-' %in% rownames(dataPrev_i)){countAbs = dataPrev_i['-',]}else{countAbs = 0} #absent from LTCF
    if('0' %in% rownames(dataPrev_i)){count0 = dataPrev_i['0',]}else{count0 = 0} # susceptible
    if('1' %in% rownames(dataPrev_i)){count1 = dataPrev_i['1',]}else{count1 = 0} # exposed
    if('2I' %in% rownames(dataPrev_i)){count2I = dataPrev_i['2I',]}else{count2I = 0} #pre-asymptomatic
    if('2II' %in% rownames(dataPrev_i)){count2II = dataPrev_i['2II',]}else{count2II = 0} #pre-symptomatic
    if('3I' %in% rownames(dataPrev_i)){count3I = dataPrev_i['3I',]}else{count3I = 0} # asymptomatic
    if('3II' %in% rownames(dataPrev_i)){count3II = dataPrev_i['3II',]}else{count3II = 0} # mild symptomatic
    if('3III' %in% rownames(dataPrev_i)){count3III = dataPrev_i['3III',]}else{count3III = 0} # severe symptomatic
    if('4' %in% rownames(dataPrev_i)){count4 = dataPrev_i['4',]}else{count4 = 0} # recovered / immunized
    
    df_dataPrev_i = data.frame(countAbs = countAbs, count0 = count0, count1 = count1,count2I = count2I, count2II = count2II, 
                               count3I = count3I, count3II = count3II, count3III = count3III, count4 = count4)%>%data.frame()%>%
      mutate(day = date)%>%
      mutate(immunity = immunity)%>%
      mutate(IPC = IPC)%>%
      mutate(Network = Network)%>%
      mutate(SimCTC = SIM)
    
    rownames(df_dataPrev_i) = date
    
    df_dataPrevSummarized = rbind(df_dataPrevSummarized, df_dataPrev_i)
  }
  
  return(df_dataPrevSummarized)
}

# # Test
# f_prevalence_summarize("lot3", "0.20.2", 1, 1, 1)



#############################
### SENSITIVITY OVER TIME ###
#############################

##############################
### f_sensitivity_by_time ###
##############################

# load data frame with PCR sensitivity since date of infection
df_kucirka_adjusted =  read.csv(paste0(filepath, "kucirka_adjusted.csv"))%>%
  mutate(Test = ifelse(Test == 'sens', 'PCR', ifelse(Test == 'sens_adj_const', 'RDT1', 'RDT2')))

df_kucirka_uniform =  read.csv(paste0(filepath, "kucirka_uniform.csv"))%>%
  mutate(Test = ifelse(Test == 'sens', 'PCR', ifelse(Test == 'sens_adj_const', 'RDT1', 'RDT2')))

df_kucirka_perfect =  read.csv(paste0(filepath, "kucirka_perfect.csv"))%>%
  mutate(Test = ifelse(Test == 'sens', 'PCR', ifelse(Test == 'sens_adj_const', 'RDT1', 'RDT2')))

#### NB: 
# PCR = RT-PCR as per Kucirka estimates
# RDT1 = Ag-RDT (A): rapid antigen testing with constant reduction in sensitivity relative to Kucirka
# RDT2 = Ag-RDT (B): rapid antigen testing with relative reduction in sensitivity depending on time from symptom onset

# function to build a PCR sensitivity matrix for each CTC status matrix
# accounting for the individual's infection status on the first day of infection (i.e. accounting for admission)
f_dataSensitivity = function(dataStatus, testType, sensitivity_function = 'time-varying'){
  
  if(!testType %in% c('PCR', 'RDT1', 'RDT2')){warning('Wrong testType'); break()}
  
  if(sensitivity_function == 'time-varying'){df_sensitivity_by_time = dplyr::filter(df_kucirka_adjusted, Test == testType)}
  if(sensitivity_function == 'uniform'){df_sensitivity_by_time = dplyr::filter(df_kucirka_uniform, Test == testType)}
  if(sensitivity_function == 'perfect'){df_sensitivity_by_time = dplyr::filter(df_kucirka_perfect, Test == testType)}
  
  dataSensitivity = as.matrix(dataStatus);
  dataSensitivity[,4:ncol(dataStatus)] = 0
  
  # for each row in the CTC output, fill a matrix with that individual's test sensitivity
  for(statusrow in 1:nrow(dataStatus)){
    
    # first select only individuals that become infected
    suivi = as.character(varhandle::unfactor(as.matrix(dataStatus[statusrow,4:ncol(dataStatus)])))
    if(T %in% (c('1','2I', '2II','3I','3II','3III') %in% suivi)){
      
      # find the days they are infected (can still include '4' because, even though people enter the population immune, cases won't)
      which_d_inf = which(suivi %in% c('1','2I', '2II','3I','3II','3III', '4'))
      
      # NB: replacement staff have same ID as individual they replace; thus in rare cases the same ID is associated with multiple distinct infections
      # if not all infected days are consecutive (so a distinct infection situation), take only the consecutive days:
      # start by setting which_d_inf1 and which_d_inf2 to 0
      which_d_inf1 = 0; which_d_inf2 = 0
      first_d_inf1 = 0; first_d_inf2 = 0
      first_d_status1 = 0; first_d_status2 = 0
      
      # then if there is a jump in infected days, split up the two infected periods
      if(T %in% c(diff(which_d_inf) %in% 2:30)){
        
        # if there are 3+ infected periods, return warning (there should not be)
        if(length(which(diff(which_d_inf)>1))>1){warning('There are 3 infections for the same patient ID')}
        
        break_infection_vec = which(diff(which_d_inf) != 1)
        which_d_inf1 = which_d_inf[1:break_infection_vec]
        which_d_inf2 = which_d_inf[(break_infection_vec+1):length(which_d_inf)]
        
        if(length(c(which_d_inf1, which_d_inf2)) != length(which_d_inf)){warning('infection dates not splitting properly'); print(statusrow)}
      }
      
      # and then also pull out the first day of infection and its status
      first_d_inf = which_d_inf[1]
      first_d_status = suivi[first_d_inf]
      if(first_d_inf != which(suivi %in% c('1','2I', '2II','3I','3II','3III'))[1]){warning('Not correctly finding first day of infection'); print(statusrow)}
      
      # or the two 'first days of infection', if multiple separate infections
      if(sum(which_d_inf2) != 0){
        first_d_inf1 = which_d_inf1[1]
        first_d_inf2 = which_d_inf2[1]
        
        first_d_status1 = suivi[first_d_inf1]
        first_d_status2 = suivi[first_d_inf2]
      }
      
      ## create a vector of test sensitivities over time
      # in CTC modeler, first day of infection status corresponds to infection the previous day,
      # so this corresponds well to first value in df_sensitivity_by_time (which is 'day 1' after infection)
      
      ## first vector if there is only one distinct infection
      # it is possible that first infection status is not '1' (e.g. for staff cases coming in from community), so must account for this
      if(sum(which_d_inf2) == 0){
        sensitivity_vector = c(rep(0,first_d_inf-1),df_sensitivity_by_time$sensitivity[1:length(which_d_inf)])
        
        if(first_d_status %in% c('2I', '2II')){
          duration_exposed = sample(2:5,1);
          dur_init = duration_exposed+1
          sensitivity_vector = c(rep(0,first_d_inf-1),df_sensitivity_by_time$sensitivity[dur_init:(dur_init+length(which_d_inf)-1)])
          
        }
        
        if(first_d_status %in% c('3I', '3II', '3III')){
          duration_exposed = sample(2:5,1);
          duration_presymptomatic = sample(1:3,1)
          dur_init = duration_exposed+duration_exposed+1
          sensitivity_vector = c(rep(0,first_d_inf-1),df_sensitivity_by_time$sensitivity[dur_init:(dur_init+length(which_d_inf)-1)])
          
        }
        
        if(first_d_status %in% c('-','0','4')){warning('first day of infection has bad status')}
        
        if(length(sensitivity_vector)<length(suivi)){
          sensitivity_vector = append(sensitivity_vector, rep(0, length(suivi) - length(sensitivity_vector)))
        }
      }else{ # and if there are two distinct infections
        
        # first half is calculated same as above, but include initial 0 sensitivities prior to infection, and add on the 'gap between infections' at the end
        gap_between_infections = first(which_d_inf2) - last(which_d_inf1) - 1
        sensitivity_vector1 = c(rep(0,first_d_inf1-1),
                                df_sensitivity_by_time$sensitivity[1:length(which_d_inf1)],
                                rep(0,gap_between_infections))
        
        # but first check if first status is not '1', and update if necessary
        if(first_d_status1 %in% c('2I','2II')){duration_exposed = sample(2:5,1); dur_init = duration_exposed+1;
        sensitivity_vector1 = c(rep(0,first_d_inf1-1),
                                df_sensitivity_by_time$sensitivity[dur_init:(dur_init+length(which_d_inf1)-1)],
                                rep(0,gap_between_infections))
        }
        if(first_d_status1 %in% c('3I','3II','3III')){duration_exposed = sample(2:5,1); duration_presymptomatic = sample(1:3,1); dur_init = duration_exposed+duration_presymptomatic+1;
        sensitivity_vector1 = c(rep(0,first_d_inf1-1),
                                df_sensitivity_by_time$sensitivity[dur_init:(dur_init+length(which_d_inf1)-1)],
                                rep(0,gap_between_infections))
        }
        
        # second half is calculated in the simplest way using the duration of infection and data from df_sensitivity_by_time
        sensitivity_vector2 = df_sensitivity_by_time$sensitivity[1:length(which_d_inf2)]
        # but check if first status is not '1'
        if(first_d_status2 %in% c('2I','2II')){duration_exposed = sample(2:5,1); 
        dur_init = duration_exposed+1;
        sensitivity_vector2 = df_sensitivity_by_time$sensitivity[dur_init:(dur_init+length(which_d_inf2)-1)]
        }
        if(first_d_status2 %in% c('3I', '3II', '3III')){duration_exposed = sample(2:5,1); duration_presymptomatic = sample(1:3,1); 
        dur_init = duration_exposed+duration_presymptomatic+1;
        sensitivity_vector2 = df_sensitivity_by_time$sensitivity[dur_init:(dur_init+length(which_d_inf2)-1)]
        }
        
        sensitivity_vector = c(sensitivity_vector1, sensitivity_vector2)
        
        # and if there is remaining lag at the end (which happens if patient is discharged),
        # then finish off with zeroes
        
        if(length(sensitivity_vector)<length(suivi)){
          sensitivity_vector = append(sensitivity_vector, rep(0, length(suivi) - length(sensitivity_vector)))
        }
      }
      
      if(length(sensitivity_vector) != length(4:ncol(dataStatus))){warning('sensitivity vector not the right length')}
      
      dataSensitivity[statusrow,4:ncol(dataStatus)]=sensitivity_vector
      
    }
  }
  return(dataSensitivity)
}


# ## TEST
# fileroot = f_fileroot("lot3", "0.20.2",1,1,1)
# dataStatus = f_read_csv(filepath_CTC, fileroot, "statusByDay")
# dataSensitivityPCR = f_dataSensitivity(dataStatus, "PCR")
# dataSensitivityRDT1 = f_dataSensitivity(dataStatus, "RDT1")
# dataSensitivityRDT2 = f_dataSensitivity(dataStatus, "RDT2")
# dataSensitivityRDT2_uniform = f_dataSensitivity(dataStatus, "RDT2", sensitivity_function = "uniform")
# dataSensitivityRDT2_perfect = f_dataSensitivity(dataStatus, "RDT2", sensitivity_function = "perfect")


################
### SYMPTOMS ###
################

### Identify the date on which individuals first present with real Covid symptoms or Covid-like symptoms (for symptomatic testing)

### NB: including Covid-like symptoms (attackrate_grippe, estimated previously using OSCOUR data, see Smith et al. BMC Medicine 2020), 
### individuals with non-symptomatic Covid infections can thus also have other Covid symptoms and be tested
### this influences testing results: same individuals potentially (but rarely) tested multiple times

f_dataSymptomIncidence = function(dataStatus, attackrate_grippe = 0.01110268, prob_severe = 0.2){

  # make empty DF to hold symptom incidence data
  dataSymptomIncidence = matrix(0, nrow = nrow(dataStatus), ncol = ncol(dataStatus))
  colnames(dataSymptomIncidence) = colnames(dataStatus)
  rownames(dataSymptomIncidence) = rownames(dataStatus)
  dataSymptomIncidence[,'id'] = as.character(dataStatus[,'id'])
  dataSymptomIncidence[,'cat'] = as.character(dataStatus[,'cat'])
  dataSymptomIncidence[,'ward'] = as.character(dataStatus[,'ward'])
  
  # loop through symptoms to add daily rate of flu symptom onset among patients
  # and to detect if first date of symptoms
  for(i in 4:(ncol(dataStatus))){
    for(j in 1:nrow(dataStatus)){
      
      # apply ILI-like symptoms
      
      loopSymptoms=0
      if(dataStatus[j,i] != '-'){
        rand = runif(1,0,1); if(rand < attackrate_grippe){
          dataSymptomIncidence[j,i] = '1'; rand2 = runif(1,0,1);
          if(rand2 < prob_severe){
            dataSymptomIncidence[j,i] = '2'
          }
        }
      }
      
      # apply actual Covid symptoms
      if(dataStatus[j,i] == '3II' & dataStatus[j,i-1] != '3II'){dataSymptomIncidence[j,i] = '1'}
      if(dataStatus[j,i] == '3III' & dataStatus[j,i-1] != '3III'){dataSymptomIncidence[j,i] = '2'}
    }
  }
  
  return(dataSymptomIncidence)
}

# ## TEST
# dataSymptoms = f_dataSymptomIncidence(dataStatus)

### quick function to calculate incidence of symptoms for any given dataSymptoms

f_symptomsIncidenceNb = function(dataSymptoms){
  nbSymptoms = length(which(dataSymptoms[,4:ncol(dataSymptoms)] != "0"))
  return(nbSymptoms)
}

# ## TEST
# dataSymptomsIncidenceNB = f_symptomsIncidenceNb(dataSymptoms)

####################
### TEST RESULTS ###
####################

### initial empty test results data.frame used throughout
df_test_results0 = data.frame('n_positive_and_tested' = 0, 'n_test' = 0, 'n_test_positive' = 0, 
                              'n_test_true_positive' = 0, 'n_test_false_positive' = 0,
                              'n_test_true_negative' = 0, 'n_test_false_negative' = 0)
if(sum(df_test_results0) != 0){warning('default df_test_results0 not properly initiated'); break()}

# Function to test selected IDs at a given time i from data_sensitivity
# and to determine whether or not the sample tested positive
# WARNING: cannot take empty vector

f_testResults = function(m_sensitivity_i_above0){
  
  if(length(m_sensitivity_i_above0) == 0){warning("trying to test empty set of IDs")}
  
  vec_rand = runif(length(m_sensitivity_i_above0)/2, 0, 1)
  
  m_test_results = data.frame(id = m_sensitivity_i_above0[,'id'], 
                              sensitivity = as.numeric(m_sensitivity_i_above0[,'sensitivity']), 
                              vec_rand = as.numeric(vec_rand))%>%
    mutate(test_result = ifelse(sensitivity > vec_rand, 1, 0))%>%
    as.matrix()
  
  return(m_test_results)
}


#############################################################################
### GENERAL TESTING FUNCTION GIVEN A SET OF IDS AND SPECIFIC COLUMN (DAY) ###
#############################################################################

f_testGivenIds = function(dataStatus, dataSensitivity, ids_tested, column){
  df_test_results = data.frame()
  # how many tested?
  n_test_i = length(ids_tested)
  
  # if no one to be tested, return all zeroes; otherwise, go through with full testing
  if(n_test_i == 0){
    
    df_test_results = data.frame('n_positive_and_tested' = 0, 'n_test' = n_test_i, 'n_test_positive' = 0, 
                                'n_test_true_positive' = 0, 'n_test_false_positive' = 0,
                                'n_test_true_negative' = 0, 'n_test_false_negative' = 0)
    ids_true_positive = c()
    list_test_results = list(df_test_results = df_test_results, ids_true_positive = ids_true_positive)
    return(list_test_results)
  }else{
    # what is status and sensitivity of those tested?
    m_status_i = dataStatus[dataStatus[,'id'] %in% ids_tested, c(1,column)]
    m_sensitivity_i = dataSensitivity[dataSensitivity[,'id'] %in% ids_tested, c(1,column)]
    
    # not recognized as matrices if only one row, so transpose character vectors
    if(class(m_status_i) == 'character'){m_status_i = t(m_status_i)}
    colnames(m_status_i) = c('id', 'status')
    if(class(m_sensitivity_i) == 'character'){m_sensitivity_i = t(m_sensitivity_i)}
    colnames(m_sensitivity_i) = c('id', 'sensitivity')
    
    ### Test!
    m_test_results_i = f_testResults(m_sensitivity_i)[,c('id', 'test_result')]
    if(class(m_test_results_i) == 'character'){m_test_results_i = t(m_test_results_i)}
    colnames(m_test_results_i) = c('id', 'test_result')
    
    ### How many positive tests?
    n_test_positive_i = sum(m_test_results_i[, 'test_result'] == "1")
    
    # calculate true positive, true negative, false positive, false negative
    df_test_results_full = merge(merge(m_status_i, m_sensitivity_i), m_test_results_i)
    df_test_results_full[,'sensitivity'] = as.numeric(as.character(df_test_results_full[,'sensitivity']))
    df_test_results_full[,'test_result'] = as.numeric(as.character(df_test_results_full[,'test_result']))
    
    df_test_results_full = df_test_results_full%>%
      mutate(true_positive = ifelse(test_result == 1 & sensitivity > 0, 1, 0))%>%
      mutate(false_positive = ifelse(test_result == 1 & status %in% c("-", "0", "4") & sensitivity == 0, 1, 0))%>%
      mutate(true_negative = ifelse(test_result == 0 & sensitivity == 0 & status %in% c("-", "0", "4"), 1, 0))%>%
      mutate(false_negative = ifelse(test_result == 0 & (sensitivity >0 | status %in% c('1', '2I', '2II', '3I', '3II', '3III')), 1, 0))
    
    # how many test positives, how many true positives 
    n_test_true_positive_i = sum(df_test_results_full$true_positive)
    n_test_false_positive_i = sum(df_test_results_full$false_positive)
    n_test_true_negative_i = sum(df_test_results_full$true_negative)
    n_test_false_negative_i = sum(df_test_results_full$false_negative)
    
    # which IDs were true positive?
    ids_true_positive = as.character(filter(df_test_results_full, true_positive == 1)$id)
    
    # how many actually positive among those tested? (back-calculated from test_true_positive + test_false_negative)
    n_tested_positive_i = n_test_true_positive_i + n_test_false_negative_i
    
    # debug testing
    if(n_test_true_positive_i + n_test_false_positive_i + n_test_true_negative_i + n_test_false_negative_i != length(ids_tested)){
      warning("Test results not equal to number of individuals tested"); break()
    }
    if(n_test_true_positive_i + n_test_false_positive_i != n_test_positive_i){
      warning("True positive + false positive don't equal total number of positive tests"); break()
    }
    if(n_test_true_positive_i > n_tested_positive_i){
      warning("More 'true positive' test results than positive samples"); break()
    }
    if(n_test_false_positive_i > 0){
      warning("Getting false positive ; I don't think I should be "); break()
    }
    
    df_test_results = data.frame('n_positive_and_tested' = n_tested_positive_i, 'n_test' = n_test_i, 'n_test_positive' = n_test_positive_i, 
                                'n_test_true_positive' = n_test_true_positive_i, 'n_test_false_positive' = n_test_false_positive_i,
                                'n_test_true_negative' = n_test_true_negative_i, 'n_test_false_negative' = n_test_false_negative_i)
    
    list_test_results = list(df_test_results = df_test_results, ids_true_positive = ids_true_positive)
    
    return(list_test_results)
  }
  warning("Broken function")
}

# ## TEST using IDs of everyone in LTCF
# f_testGivenIds(dataStatus, dataSensitivityRDT2, dataStatus[,1], 4) 

#############################################
### TEST ADMISSIONS USING COMM-ONSET DATA ###
#############################################
### identify admissions from CommCases and test them on first day
### NB: This applies to patients only (staff aren't "admitted")
### ASIDE: CommCases are all "new admissions" as reflected in TOY ADMISSION data

# First, calculate total number of admissions from admissions file to know the number of tests used for admissions
# NB: this only needs to be done once for ALL simulations (identical admissions for each, uses real admission data from Berck-sur-Mer)
admissions_all = read.csv2(paste0(filepath, "toy_admission.csv"))%>%
  mutate(Date = as.Date(firstDate))%>%
  filter(Date >= as.Date("2009-07-06") & Date <= as.Date("2009-07-20"))%>%
  filter(status == 'PA')
numAdmissionsPA = nrow(admissions_all)



f_testAdmissions = function(dataCommCases, dataStatus, dataSensitivity, testType){
  
  # already initiated df_rest_results0
  if(sum(df_test_results0) != 0){warning('default df_test_results0 not properly initiated'); break()}
  
  if(!testType %in% c('PCR', 'RDT1', 'RDT2')){warning("testType specified as neither 'PCR' nor 'RDT' so can't assess"); break()}
  
  if(testType == 'PCR'){lagTestResult = 1}
  if(testType %in% c('RDT1', 'RDT2')){lagTestResult = 0}
  
  dates_true_positive_admissions = c()
  ids_true_positive_admissions = c()
  df_test_results_admissions = df_test_results0
  
  dataCommCasesPatients = dataCommCases%>%
    mutate(cat = ifelse(str_detect(Individual, "PA-"), 'Patient', 'Staff'))%>%
    filter(cat == 'Patient')
  
  # if no community-onset cases, just return the empty list
  if(nrow(dataCommCasesPatients)==0){
    
    df_test_results_admissions$n_test = numAdmissionsPA
    df_test_results_admissions$n_test_true_negative = numAdmissionsPA
    
    list_test_results_admissions = list(df_test_results = df_test_results_admissions, 
                                   ids_true_positive = ids_true_positive_admissions,
                                   dates_true_positive = dates_true_positive_admissions)
    return(list_test_results_admissions)
    
  }else{ ### put in true negative results already for all negative admissions (PATIENTS ONLY)
    
    df_test_results_admissions$n_positive_and_tested = 0
    df_test_results_admissions$n_test = numAdmissionsPA - nrow(dataCommCasesPatients)
    df_test_results_admissions$n_test_true_negative = numAdmissionsPA - nrow(dataCommCasesPatients)
    
  }
  

  # what are the IDs and Dates of community-onset patient cases (to test individually in following loop)
  vec_IDsAdmissions = as.character(dataCommCasesPatients$Individual)
  vec_DatesAdmission = as.Date(dataCommCasesPatients$Date)
  
  
  for(index in 1:length(vec_DatesAdmission)){
    
    ids = vec_IDsAdmissions[index]
    date = vec_DatesAdmission[index]
    
    colnumber = as.numeric(4+date - as.Date("2009-07-06"))
    if(colnumber == 18){colnumber = 17}# can't add +1 to colnumber below if already the last column (18)
    
    # TEST: is dataStatus updated in time to catch admissions on same day? No! So actually use next day for status and sensitivity
    # this will return results that reflect real status on admission day
    # NB: this is only for admission testing; symptoms and daily testing are able to take appropriate same-day estimates
    df_testGivenIdsIndex = f_testGivenIds(dataStatus, dataSensitivity, ids, colnumber+1)
    df_test_results_admissions_index = df_testGivenIdsIndex$df_test_results
    ids_true_positive_admissions_index = df_testGivenIdsIndex$ids_true_positive
    if(length(ids_true_positive_admissions_index)>0){dates_true_positive_admissions_index=colnumber + lagTestResult}else{
      dates_true_positive_admissions_index=ids_true_positive_admissions_index
    }
    
    df_test_results_admissions = df_test_results_admissions+df_test_results_admissions_index
    ids_true_positive_admissions = append(ids_true_positive_admissions, ids_true_positive_admissions_index)
    dates_true_positive_admissions = append(dates_true_positive_admissions, dates_true_positive_admissions_index)
  }
  
  list_test_results_admissions = list(df_test_results = df_test_results_admissions, 
                                      ids_true_positive = ids_true_positive_admissions,
                                      dates_true_positive = dates_true_positive_admissions)
  return(list_test_results_admissions)
}


# ## TEST
# f_testAdmissions(dataCommCases, dataStatus, dataSensitivityPCR, "PCR")



###########################################################
### TEST SYMPTOMATICS BY LOOPING THROUGH DAILY SYMPTOMS ###
###########################################################

## test symptomatics; loop through daily looking for symptoms
# initialize pcr results

f_testSymptoms = function(dataStatus, dataSensitivity, dataSymptoms, testType){
  
  if(!testType %in% c('PCR', 'RDT1', 'RDT2')){warning("testType specified as neither 'PCR' nor 'RDT' so can't assess"); break()}
  
  if(testType == 'PCR'){lagTestResult = 1}
  if(testType %in% c('RDT1', 'RDT2')){lagTestResult = 0}
  
  ### INITIALIZE EMPTY RESULTS DATAFRAME ###
  # already initiated df_rest_results0
  if(sum(df_test_results0) != 0){warning('default df_test_results0 not properly initiated'); break()}
  
  dates_true_positive = c()
  ids_true_positive_symptoms = c()
  df_test_results_symptoms = df_test_results0
  
  for(simday in 4:ncol(dataSymptoms)){
    # who is tested? only symptomatics
    ids_symptoms = c(dataSymptoms[dataSymptoms[,simday] %in% c("1","2"),'id'])
    # get test results
    list_test_results_symptoms = f_testGivenIds(dataStatus, dataSensitivity, ids_symptoms, simday)
    # add these daily results to running overall results
    df_test_results_symptoms = df_test_results_symptoms + list_test_results_symptoms$df_test_results
    # select IDs of positives
    ### !! NB these include individuals repeatedly tested, use unique below to only select unique cases
    ids_true_positive_symptoms = append(ids_true_positive_symptoms, list_test_results_symptoms$ids_true_positive)
    
    # save date of positive test; date will have same index as corresponding ID saved in same list
    if(length(list_test_results_symptoms$ids_true_positive)>0){
      dates_true_positive = append(dates_true_positive, rep(simday + lagTestResult, length(list_test_results_symptoms$ids_true_positive)))}
  }
  
  
  list_test_results_final = list(df_test_results = df_test_results_symptoms, 
                                ids_true_positive = ids_true_positive_symptoms,
                                dates_true_positive = dates_true_positive)
  return(list_test_results_final)
}

# ## TEST
# f_testSymptoms(dataStatus, dataSensitivityPCR, dataSymptoms, 'PCR')

######################################
### UNIVERSAL TESTING EVERY X DAYS ###
######################################

f_testEveryoneDaily = function(dataStatus, dataSensitivity, column, testType, testTarget){
  
  if(!testType %in% c('PCR', 'RDT1', 'RDT2')){warning("testType specified as neither 'PCR' nor 'RDT' nor 'RDT2' so can't assess"); break()}
  if(!testTarget %in% c('Patients', 'Staff', 'All')){warning("testTarget not specified as patients, staff or all"); break()}
  
  if(testType == 'PCR'){lagTestResult = 1}
  if(testType %in% c('RDT1', 'RDT2')){lagTestResult = 0}
  
  # for selected day of simulation (column), select IDs to be tested and evaluate these using test

  # who is tested? everyone present, staff or patients
  if(testTarget == "Patients"){
    dataStatusPatients = dplyr::filter(dataStatus, cat == "Patient")
    ids_tested = as.character(dataStatusPatients[dataStatusPatients[,column] != c("-"),'id'])
  }
  if(testTarget == "Staff"){
    dataStatusStaff = dplyr::filter(dataStatus, cat != "Patient")
    ids_tested = as.character(dataStatusStaff[dataStatusStaff[,column] != c("-"),'id'])
  }
  if(testTarget == "All"){
    ids_tested = as.character(dataStatus[dataStatus[,column] != c("-"),'id'])
  }
  
  # get test results
  list_test_results_i = f_testGivenIds(dataStatus, dataSensitivity, ids_tested, column)

  # and format test results into standard data.frame output format
  list_test_results = list(df_test_results = list_test_results_i$df_test_results,
                           ids_true_positive = list_test_results_i$ids_true_positive,
                           dates_true_positive = rep(column+lagTestResult, length(list_test_results_i$ids_true_positive)))
  return(list_test_results)
}

# ## TEST
# f_testEveryoneDaily(dataStatus, dataSensitivityRDT2, 4, "RDT2", "Patients")

############################################################
### COMBINING RESULTS FOR MULTI-LEVEL TESTING STRATEGIES ###
############################################################

# input a list of each level of testing results, and return "combined" testing results, 
# with full numbers of tests/positive/negatives; IDs of positive cases; and "dates" (actually column numbers) of detection

f_combineTestResults = function(testResults){
  ### testResults should be a LIST of the different test results
  if(class(testResults) != 'list'){warning('Error: need a LIST of test results'); break()}
  
  # already initiated df_rest_results0
  if(sum(df_test_results0) != 0){warning('default df_test_results0 not properly initiated'); break()}
  
  # initialize test results as taking default df
  df_test_results = df_test_results0
  
  ids_true_positive = c()
  dates_true_positive = c()
  
  for(testResults_i in testResults){
    df_test_results_i = testResults_i$df_test_results
    df_test_results = df_test_results+df_test_results_i
    
    ids_true_positive_i = testResults_i$ids_true_positive
    ids_true_positive = c(ids_true_positive, ids_true_positive_i)
    
    dates_true_positive_i = testResults_i$dates_true_positive
    dates_true_positive = c(dates_true_positive, dates_true_positive_i)
  }
  
  listResults = list(df_test_results = df_test_results,
                     ids_true_positive = ids_true_positive,
                     dates_true_positive = dates_true_positive)
  
  return(listResults)
}

# ## TEST
# testResults_symptomaticTesting = f_testSymptoms(dataStatus, dataSensitivityPCR, dataSymptoms, 'PCR')
# testResults_dailyScreening = f_testEveryoneDaily(dataStatus, dataSensitivityRDT2, 10, "RDT2", "Patients")
# testResults_multiLevel = f_combineTestResults(list(testResults_symptomaticTesting, testResults_dailyScreening))

#################################################################
### IDENTIFY NOSOCOMIAL TRANSMITTERS AND DATE OF TRANSMISSION ### 
#################################################################
### NB: secondaryCases data must be read with stringsAsFactors = F

f_identifyTransmitters = function(dataSecondaryCases){
  
  if(nrow(dataSecondaryCases) == 0){
    dataTransmit = data.frame()
    return(dataTransmit)
    }
  
  dataTransmit_raw = data.frame(dataSecondaryCases)%>%
    mutate(Date = as.Date(Date))%>%
    dplyr::select(Individual, Date, Cases)
  dataTransmit_raw$Cases = as.character(dataTransmit_raw$Cases)
  
  dataTransmit_temp <- data.table(dataTransmit_raw)
  dataTransmit <- dataTransmit_temp[, list(Individual = Individual, 
                                                       Date = Date,
                                                       Cases = unlist(strsplit(Cases, "_"))), 
                                                by=1:nrow(dataTransmit_temp)]%>%
    mutate(date_transmission = as.numeric(Date - as.Date('2009-07-05')))%>%
    mutate(colnum_transmission = date_transmission+3)%>%
    dplyr::select(Individual, Cases, date_transmission, colnum_transmission)
  
  colnames(dataTransmit) <- c("id_infector", "id_infectee", "date_transmission", "colnum_transmission")
  
  return(dataTransmit)
}

# ## TEST
# dataTransmit = f_identifyTransmitters(dataSecondaryCases)

####################################################################
### REMOVE NOSOCOMIAL TRANSMITTERS AND ALL DOWNSTREAM INFECTIONS ###
####################################################################

# Take transmitor data, and the timing and success data from surveillance (i.e. which individuals identified) 
# and return all downstream infections averted by blocking nosocomial transmission from these individuals

### NB: inputted dates are in terms of COLUMN NUMBER

f_infectionsAverted = function(dataTransmit, resultsStrat){
  
  ids_detected = resultsStrat$ids_true_positive;
  colnum_surveillance = resultsStrat$dates_true_positive
  
  if(length(ids_detected) != length(colnum_surveillance)){warning("Error: number of positive test result IDs and dates don't match"); break()}
  
  ## Initialize intermediate and final outputs: IDs_averted_overall
  ids_averted_overall = c()
  ids_averted_direct = c()
  ids_averted_indirect = c()
  ids_averted_indirect2 = c()
  ids_averted_indirect3 = c()
  ids_averted_indirect4 = c()
  ids_averted_indirect5 = c()
  ids_averted_indirect6 = c()
  ids_averted_indirect7 = c()
  
  ### Loop through each ID/date combo
  for(index in 1:length(ids_detected)){
    
    ids_detected_i = ids_detected[index]
    colnum_surveillance_i = as.numeric(colnum_surveillance[index])
    
    ### NB: Assume test result is at beginning of day, so transmissions on that SAME day are avoided (hence colnum_transmission >= colnum_surveillance_i)
    inf_averted_direct = dataTransmit[which(dataTransmit$id_infector %in% ids_detected_i & dataTransmit$colnum_transmission >= colnum_surveillance_i),]
    ids_averted_direct = c(inf_averted_direct$id_infectee)
    
    if(length(ids_averted_direct) > 0){
      inf_averted_indirect = dataTransmit[which(dataTransmit$id_infector %in% ids_averted_direct & dataTransmit$colnum_transmission >= colnum_surveillance_i),]
      ids_averted_indirect = c(inf_averted_indirect$id_infectee)
      
      if(length(ids_averted_indirect) > 0){
        inf_averted_indirect2 = dataTransmit[which(dataTransmit$id_infector %in% ids_averted_indirect & dataTransmit$colnum_transmission >= colnum_surveillance_i),]
        ids_averted_indirect2 = c(inf_averted_indirect2$id_infectee)
        
        if(length(ids_averted_indirect2) > 0){
          inf_averted_indirect3 = dataTransmit[which(dataTransmit$id_infector %in% ids_averted_indirect2 & dataTransmit$colnum_transmission >= colnum_surveillance_i),]
          ids_averted_indirect3 = c(inf_averted_indirect3$id_infectee)
          
          if(length(ids_averted_indirect3) > 0){
            inf_averted_indirect4 = dataTransmit[which(dataTransmit$id_infector %in% ids_averted_indirect3 & dataTransmit$colnum_transmission >= colnum_surveillance_i),]
            ids_averted_indirect4 = c(inf_averted_indirect4$id_infectee)
            
            if(length(ids_averted_indirect4) > 0){
              inf_averted_indirect5 = dataTransmit[which(dataTransmit$id_infector %in% ids_averted_indirect4 & dataTransmit$colnum_transmission >= colnum_surveillance_i),]
              ids_averted_indirect5 = c(inf_averted_indirect5$id_infectee)
              
              if(length(ids_averted_indirect5) > 0){
                inf_averted_indirect6 = dataTransmit[which(dataTransmit$id_infector %in% ids_averted_indirect5 & dataTransmit$colnum_transmission >= colnum_surveillance_i),]
                ids_averted_indirect6 = c(inf_averted_indirect6$id_infectee)
                
                if(length(ids_averted_indirect6) > 0){
                  inf_averted_indirect7 = dataTransmit[which(dataTransmit$id_infector %in% ids_averted_indirect6 & dataTransmit$colnum_transmission >= colnum_surveillance_i),]
                  ids_averted_indirect7 = c(inf_averted_indirect7$id_infectee)
                  
                  if(length(ids_averted_indirect7) > 0){
                    warning("Re-infection in simulation, assume both prevented");
                  }
                }
              }
            }
          }
        }
      }
    }
    
    ids_averted_overall_i = c(ids_averted_direct, ids_averted_indirect, ids_averted_indirect2, ids_averted_indirect3, 
                              ids_averted_indirect4, ids_averted_indirect5, ids_averted_indirect6, ids_averted_indirect7)
    
    ids_averted_overall = unique(c(ids_averted_overall, ids_averted_overall_i))
  }
  
  return(ids_averted_overall)
}

# ## TEST
# f_infectionsAverted(dataTransmit, testResults_multiLevel)

#########################
### SURVEILLANCE LOOP ###
#########################

### Function to put together all of the above and evaluate surveillance

### For a given CTC simulation, run the full range of surveillance and do so testing_rounds times (default = 100)
f_surveillanceLoop = function(lot, immunity, IPC, SIM, Network, 
                              testing_rounds, screening_types, screening_targets, surveillance_strategies, sens_function){
  
  if(!sens_function %in% c('time-varying', 'uniform', 'perfect')){warning("sensitivity function must be 'time-carying', 'uniform' or 'perfect'"); break()}
  
  if(length(SIM)>1){warning("Putting in multiple simulations at once; not currently possible")}
  
  ### Basic print commands to orient running of function
  print(paste0(testing_rounds, " runs of surveillance per CTC simulation"))
  
  ### initialize vectors that will be looped through to produce results
  vec_testing_rounds = 1:testing_rounds
  vec_screening_types = screening_types
  vec_screening_targets = screening_targets
  vec_surveillance_strategies = surveillance_strategies
  
  ### initialize column names corresponding to surveillanceMatrix
  vec_ncols = c(
    ### which CTC simulation is this? including immunity, IPC, Network
    "immunity", "IPC", "Network", "simulationCTC", 
    ### INCIDENCE; total incidence over full simulation period; nb index cases;  nb subsequent comm-onset cases; nb within-hospital nosocomial acquisitions
    "incSarscov2", "incIndex", "incCommCases", "incNosocomial",  
    ### specifities of surveillance simulations: which surveillance simulation is this?; what is incidence of symptoms (covid + non-covid); what type of test used for daily screening?;  who is targeted by daily screening? Patients, staff or everyone?; what surveillance strategy is used?
    "testingRound", "incSymptoms", "typeDailyScreening", "targetDailyScreening", "stratSurveillance", 
    ### test results (numbers of tests, positives, negatives, etc. for PCR)
    "PCR.n_positive_and_tested", "PCR.n_test","PCR.n_test_positive","PCR.n_test_true_positive","PCR.n_test_false_positive","PCR.n_test_true_negative","PCR.n_test_false_negative",
    ### test results (numbers of tests, positives, negatives, etc. for RDT)
    "RDT.n_positive_and_tested", "RDT.n_test","RDT.n_test_positive","RDT.n_test_true_positive","RDT.n_test_false_positive","RDT.n_test_true_negative","RDT.n_test_false_negative",
    ### FINAL OUTCOME: number of infections averted through testing
    "infAverted"
  )
  
  which_cols_PCR = which(str_detect(vec_ncols, "PCR."))
  which_cols_RDT = which(str_detect(vec_ncols, "RDT."))
  
  ### initiate MATRIX for surveillance results
  m_surveillance_results = matrix(ncol = length(vec_ncols), 
                                  nrow = length(vec_testing_rounds)*length(vec_screening_types)*length(vec_screening_targets)*length(vec_surveillance_strategies))
  
  
  ### DATA
  ### Read data: prepare file root
  fileroot = f_fileroot(lot, immunity,IPC,SIM,Network)
  
  print(fileroot)
  
  # Load and process epidemiological data
  dataStatus = f_read_csv(filepath_CTC, fileroot, "statusByDay")
  dataSecondaryCases = f_read_csv(filepath_CTC, fileroot, "secondaryCases")
  dataIndex = f_read_csv(filepath_CTC, fileroot, "FirstIndex")
  dataCommCases = f_read_csv(filepath_CTC, fileroot, 'CommCases')
  dataTransmit = f_identifyTransmitters(dataSecondaryCases)
  
  n_cases_index = nrow(dataIndex)
  n_cases_commonset = nrow(dataCommCases)
  n_cases_nosocomial = nrow(dataSecondaryCases)
  n_cases_total = n_cases_nosocomial+n_cases_index+n_cases_commonset
  
  # Testing sensitivity
  dataSensitivityPCR = f_dataSensitivity(dataStatus, "PCR", sensitivity_function = sens_function)
  dataSensitivityRDT1 = f_dataSensitivity(dataStatus, "RDT1", sensitivity_function = sens_function)
  dataSensitivityRDT2 = f_dataSensitivity(dataStatus, "RDT2", sensitivity_function = sens_function)
  
  ### initiate qounter for working through rows to fill matrix
  qounter = 0
  
  ### LOOP: in this order
  # testing_rounds (normally 1:100)
  # screening_type (PCR, RDT1, RDT2)
  # screening_target (patients, staff, all)
  # surveillance_strategy (which combination of testing)
  for(round_i in vec_testing_rounds){
    
    print(paste0("On testing round ", round_i, " of ", length(vec_testing_rounds), " rounds."))
    
    # Epidemiological data: symptoms are stochastic for each testing round
    dataSymptoms = f_dataSymptomIncidence(dataStatus)
    n_symptoms = f_symptomsIncidenceNb(dataSymptoms)
    
    # Testing results: PCR
    resultsAdmissions = f_testAdmissions(dataCommCases, dataStatus, dataSensitivityPCR, 'PCR')
    resultsSymptoms = f_testSymptoms(dataStatus, dataSensitivityPCR, dataSymptoms, 'PCR')
    
    for(screening_type_i in vec_screening_types){
      
      # pick sensitivity results based on chosen testing technology
      
      if(screening_type_i == "PCR"){dataSensitivity_i = dataSensitivityPCR}
      if(screening_type_i == "RDT1"){dataSensitivity_i = dataSensitivityRDT1}
      if(screening_type_i == "RDT2"){dataSensitivity_i = dataSensitivityRDT2}
      if(!screening_type_i %in% c("PCR", "RDT1", "RDT2")){warning('Not picking any sensitivity function')}
      
      for(screening_target_i in vec_screening_targets){
        results_d01 = f_testEveryoneDaily(dataStatus, dataSensitivity_i, 1+3, screening_type_i, screening_target_i)
        results_d02 = f_testEveryoneDaily(dataStatus, dataSensitivity_i, 2+3, screening_type_i, screening_target_i)
        results_d03 = f_testEveryoneDaily(dataStatus, dataSensitivity_i, 3+3, screening_type_i, screening_target_i)
        results_d04 = f_testEveryoneDaily(dataStatus, dataSensitivity_i, 4+3, screening_type_i, screening_target_i)
        results_d05 = f_testEveryoneDaily(dataStatus, dataSensitivity_i, 5+3, screening_type_i, screening_target_i)
        results_d06 = f_testEveryoneDaily(dataStatus, dataSensitivity_i, 6+3, screening_type_i, screening_target_i)
        results_d07 = f_testEveryoneDaily(dataStatus, dataSensitivity_i, 7+3, screening_type_i, screening_target_i)
        results_d08 = f_testEveryoneDaily(dataStatus, dataSensitivity_i, 8+3, screening_type_i, screening_target_i)
        results_d09 = f_testEveryoneDaily(dataStatus, dataSensitivity_i, 9+3, screening_type_i, screening_target_i)
        
        for(strat_i in vec_surveillance_strategies){
          
          qounter = qounter + 1
          
          if(strat_i == "adm"){listResults = list(resultsAdmissions)} # test new admissions only
          if(strat_i == "sym"){listResults = list(resultsSymptoms)} # test symptomatics only
          if(strat_i == "adm_sym"){listResults = list(resultsAdmissions, resultsSymptoms)} # routine testing as defined in manuscript (admission + symptoms)
          if(strat_i == "d01"){listResults = list(results_d01)} # test everyone on day 1
          if(strat_i == "d02"){listResults = list(results_d02)} # ... day 2
          if(strat_i == "d03"){listResults = list(results_d03)} # ... day 3
          if(strat_i == "d04"){listResults = list(results_d04)} # ... day 4
          if(strat_i == "d05"){listResults = list(results_d05)} # ... day 5
          if(strat_i == "d06"){listResults = list(results_d06)} # ... day 6
          if(strat_i == "d07"){listResults = list(results_d07)} # ... day 7
          if(strat_i == "d08"){listResults = list(results_d08)} # ... day 8
          if(strat_i == "d09"){listResults = list(results_d09)} # ... day 9
          if(strat_i == "adm_sym_d01"){listResults = list(resultsAdmissions, resultsSymptoms, results_d01)} # from here down, combine above strategies as listed
          if(strat_i == "adm_sym_d02"){listResults = list(resultsAdmissions, resultsSymptoms, results_d02)}
          if(strat_i == "adm_sym_d03"){listResults = list(resultsAdmissions, resultsSymptoms, results_d03)}
          if(strat_i == "adm_sym_d04"){listResults = list(resultsAdmissions, resultsSymptoms, results_d04)}
          if(strat_i == "adm_sym_d05"){listResults = list(resultsAdmissions, resultsSymptoms, results_d05)}
          if(strat_i == "adm_sym_d06"){listResults = list(resultsAdmissions, resultsSymptoms, results_d06)}
          if(strat_i == "adm_sym_d07"){listResults = list(resultsAdmissions, resultsSymptoms, results_d07)}
          if(strat_i == "adm_sym_d08"){listResults = list(resultsAdmissions, resultsSymptoms, results_d08)}
          if(strat_i == "adm_sym_d09"){listResults = list(resultsAdmissions, resultsSymptoms, results_d09)}
          if(strat_i == "adm_sym_d01_d02"){listResults = list(resultsAdmissions, resultsSymptoms, results_d01, results_d02)}
          if(strat_i == "adm_sym_d01_d03"){listResults = list(resultsAdmissions, resultsSymptoms, results_d01, results_d03)}
          if(strat_i == "adm_sym_d01_d04"){listResults = list(resultsAdmissions, resultsSymptoms, results_d01, results_d04)}
          if(strat_i == "adm_sym_d01_d05"){listResults = list(resultsAdmissions, resultsSymptoms, results_d01, results_d05)}
          if(strat_i == "adm_sym_d01_d06"){listResults = list(resultsAdmissions, resultsSymptoms, results_d01, results_d06)}
          if(strat_i == "adm_sym_d01_d07"){listResults = list(resultsAdmissions, resultsSymptoms, results_d01, results_d07)}
          if(strat_i == "adm_sym_d01_d08"){listResults = list(resultsAdmissions, resultsSymptoms, results_d01, results_d08)}
          if(strat_i == "adm_sym_d01_d09"){listResults = list(resultsAdmissions, resultsSymptoms, results_d01, results_d09)}
          
          
          ### FINAL RESULTS: combining results using f_combineTestResults, and determine number of infections averted
          combinedResults = f_combineTestResults(listResults)
          
          infAverted = length(f_infectionsAverted(dataTransmit, combinedResults))
          
          
          ### FILLING UP MATRIX WITH RESULTS
          # if all testing is PCR, no RDT results
          if(screening_type_i == 'PCR'){
            PCR.results = combinedResults$df_test_results;
            RDT.results = df_test_results0}
          
          # if screening is RDT, then combination of PCR/RDT will depend on the strategy
          if(screening_type_i %in% c('RDT1', 'RDT2')){
            if(strat_i %in% c("adm", "sym", "adm_sym")){
              PCR.results = combinedResults$df_test_results;
              RDT.results = df_test_results0
            }
            if(strat_i %in% c(paste0("d0",1:9))){
              PCR.results = df_test_results0;
              RDT.results = combinedResults$df_test_results;
            }
            if(strat_i %in% c(paste0("adm_sym_d0", 1:9), paste0("adm_sym_d01_d0", 2:9))){
              PCR.results = f_combineTestResults(list(resultsAdmissions, resultsSymptoms))$df_test_results
              RDT.results = combinedResults$df_test_results - f_combineTestResults(list(resultsAdmissions, resultsSymptoms))$df_test_results
            }
            if(!strat_i %in% vec_surveillance_strategies){warning('strategy used not in vec_surveillance_strategies')}
          }
          
          
          ### fill up surveillance results
          m_surveillance_results[qounter, 1:c(which_cols_PCR[1]-1)] = c(
            immunity, IPC, Network, SIM, # which CTC simulation is this
            n_cases_total, # total incidence of SARS-CoV-2 up to 15 days
            n_cases_index, # total number of index cases on day=1
            n_cases_commonset, # total number of subsequent community-onset cases
            n_cases_nosocomial, # total number of within-hospital transmission
            round_i, # which round of stochastic testing is this?
            n_symptoms, # what is the incidence of symptoms (varies per stochastic testing round)
            screening_type_i, # what is the type of daily screening used?
            screening_target_i, # who is targeted by daily screening? patients, staff or all
            strat_i # what is the surveillance strategy used here?
          )
          m_surveillance_results[qounter, which_cols_PCR] = as.matrix(PCR.results)
          m_surveillance_results[qounter, which_cols_RDT] = as.matrix(RDT.results)
          
          ### and lastly, the number of infections averted
          m_surveillance_results[qounter, which(vec_ncols == "infAverted")] = infAverted
          
        }
      }
    }
  }
  # name column names and return surveillance matrix
  colnames(m_surveillance_results) = vec_ncols
  return(m_surveillance_results)
}

# ## TEST using small run, just 5 surveillance loops and 3 surveillance strategies
# f_surveillanceLoop(lot = "lot3", immunity = "0.20.2", IPC = 1, SIM = 1, Network = 1, testing_rounds = 5, screening_types = 'RDT2',
#                    screening_targets = 'Patients', surveillance_strategies = c('sym', 'd01', 'adm_sym_d01'), sens_function = 'time-varying')


####################################
### SURVEILLANCE DATA PROCESSING ###
####################################

### LOAD SURVEILLANCE OUTPUT DATA ###

# read raw surveillance data, make into df and calculate propAverted and IRR
f_read_surveillance_output = function(lot, immunity, network, ipc, sim, sensfunc){
  
  filepath_output = paste0(filepath, "Surveillance/output/Output_",lot,"/")
  
  # in lot 1, there is no sensfunc on filename
  if(lot == "lot1"){filename = paste0("m_surveillance_", immunity, "_network", network, "_IPC", ipc, "_SIM", sim,".csv")}else{
    filename = paste0("m_surveillance_", immunity, "_network", network, "_IPC", ipc, "_SIM", sim, "_SENS", sensfunc,".csv")
  }
  
  dat_surveillance = read.csv(paste0(filepath_output, filename))%>%
    data.frame()%>%
    mutate(propAverted = infAverted/incNosocomial)%>%
    mutate(incNosocomialAfterTesting = incNosocomial - infAverted)
  return(dat_surveillance)
}

# ## TEST
# f_read_surveillance_output("lot3", "0.20.2", 1, 1, 5, "time-varying")

############################################
### Function prepare surveillance output ### 
############################################

### f_false_positive, accounting for imperfect test specificity using literature estimates (see manuscript)
# test specificity put in retrospectively (not in surveillance loop)

f_false_positives_rdt = Vectorize(function(n_true_negatives, specificity = 0.997){
  n_false_positives = rbinom(1,n_true_negatives, 1-specificity)
  return(n_false_positives)
})

f_false_positives_pcr = Vectorize(function(n_true_negatives, specificity = 0.999){
  n_false_positives = rbinom(1,n_true_negatives, 1-specificity)
  return(n_false_positives)
})


### takes about a minute per scenario when including all 100 CTC simulations: putting results for every single CTC simulation into one big matrix


f_prepare_surveillance_output = function(lot, num_CTCsims, num_testingRounds,
                                         immunity, network, ipc, stratsSurveillance, typeScreening, targetScreening, sensFunc){
  
  # how many rows in matrix for each CTC simulation
  length_per_file = num_CTCsims*num_testingRounds*length(stratsSurveillance)*length(typeScreening)*length(targetScreening)
  
  m_dat = matrix(nrow = length_per_file,ncol=19)
  colnames(m_dat) = c('simulationCTC', 'incNosocomial', 'stratSurveillance', 'testingRound', 'infAverted', 'typeDailyScreening', 'targetDailyScreening',
                      'PCR.n_positive_and_tested', 'PCR.n_test', 'PCR.n_test_true_positive', 'PCR.n_test_false_positive', 'PCR.n_test_true_negative', 'PCR.n_test_false_negative',
                      'RDT.n_positive_and_tested', 'RDT.n_test', 'RDT.n_test_true_positive', 'RDT.n_test_false_positive', 'RDT.n_test_true_negative', 'RDT.n_test_false_negative')
  
  counter = 1
  
  for(n_sims in 1:num_CTCsims){
    print(paste0("CTC simulation ", n_sims))
    
    dat = f_read_surveillance_output(lot, immunity, network, ipc, n_sims, sensFunc)%>%
      filter(typeDailyScreening %in% typeScreening, targetDailyScreening %in% targetScreening, stratSurveillance %in% stratsSurveillance)%>%
      # add false positives to RDT, and adjust true negatives and overall # of positive vs. negative tests
      mutate(RDT.n_test_false_positive = f_false_positives_rdt(RDT.n_test_true_negative),
             RDT.n_test_true_negative = RDT.n_test_true_negative - RDT.n_test_false_positive,
             RDT.n_test_positive = RDT.n_test_positive + RDT.n_test_false_positive)%>%
      # and again, add false positives to PCR
      mutate(PCR.n_test_false_positive = f_false_positives_pcr(PCR.n_test_true_negative),
             PCR.n_test_true_negative = PCR.n_test_true_negative - PCR.n_test_false_positive,
             PCR.n_test_positive = PCR.n_test_positive + PCR.n_test_false_positive)%>%
      dplyr::select(simulationCTC, incNosocomial, stratSurveillance, testingRound, infAverted, typeDailyScreening, targetDailyScreening,
                    PCR.n_positive_and_tested, PCR.n_test, PCR.n_test_true_positive, PCR.n_test_false_positive, PCR.n_test_true_negative, PCR.n_test_false_negative,
                    RDT.n_positive_and_tested, RDT.n_test, RDT.n_test_true_positive, RDT.n_test_false_positive, RDT.n_test_true_negative, RDT.n_test_false_negative)
    
    m_dat[counter:(counter+(length_per_file/num_CTCsims)-1),1:ncol(dat)] = as.matrix(dat)
    
    counter = counter + nrow(dat)
  }
  
  df_dat = data.frame(m_dat, stringsAsFactors = F)
  df_dat[,c(1,2,4,5,8:ncol(df_dat))] <- sapply(df_dat[,c(1,2,4,5,8:ncol(df_dat))], as.numeric)
  
  df_dat_updated = df_dat%>%  
    mutate(infAverted = as.numeric(as.character(infAverted)))%>%
    mutate(incNosocomial = as.numeric(as.character(incNosocomial)))%>%
    mutate(propAverted = ifelse(incNosocomial == 0, 0, infAverted/incNosocomial))%>%
    mutate(PCR.n_negative_and_tested = PCR.n_test - PCR.n_positive_and_tested,
           RDT.n_negative_and_tested = RDT.n_test - RDT.n_positive_and_tested)%>%
    mutate(PCR.sensitivity_realized = PCR.n_test_true_positive/PCR.n_positive_and_tested,
           PCR.specificity_realized = PCR.n_test_true_negative/PCR.n_negative_and_tested,
           RDT.sensitivity_realized = RDT.n_test_true_positive/RDT.n_positive_and_tested,
           RDT.specificity_realized = RDT.n_test_true_negative/RDT.n_negative_and_tested,
           PCR.PPV = PCR.n_test_true_positive/(PCR.n_test_true_positive + PCR.n_test_false_positive),
           PCR.NPV = PCR.n_test_true_negative/(PCR.n_test_true_negative + PCR.n_test_false_negative),
           RDT.PPV = RDT.n_test_true_positive/(RDT.n_test_true_positive + RDT.n_test_false_positive),
           RDT.NPV = RDT.n_test_true_negative/(RDT.n_test_true_negative + RDT.n_test_false_negative))%>%
    mutate(sensitivity_function = sensFunc)
  

  return(df_dat_updated)
  
}

### TESTS
### In this case these have to really reflect what went into these 'lots' of simulations

# test_strats = c('adm_sym', 'd01', 'adm_sym_d01_d02')
# test_strats2 = c('adm_sym', 'adm_sym_d01', 'adm_sym_d01_d02', 'adm_sym_d01_d03', 'adm_sym_d01_d04', 'adm_sym_d01_d05', 'adm_sym_d01_d06', 'adm_sym_d01_d07', 'adm_sym_d01_d08', 'adm_sym_d01_d09')

# test1 = f_prepare_surveillance_output("lot3", 5, 100, "0.20.2", 1, 1, test_strats, c('PCR','RDT2'), c('All','Patients'), "time-varying")
# test2 = f_prepare_surveillance_output("lot3", 5, 100, "0.20.2", 1, 1, test_strats2, 'RDT2', 'All', "uniform")

#######################################################################
### Functions to calculate outcomes from prepared surveillance data ###
#######################################################################

### Nb: a range of uncertainty measures were considered for each outcome:
# quantiles, standard error using classical formulae for different outcomes, standard deviation of all 10,000 simulations, standard deviation of means for each CTC simulation
# many of these problematic due to bimodal nature of infection prevention data
# ultimately bootstrap confidence intervals were determined to be most robust and consistent measure of uncertainty across outcomes

###########################
### BOOTSTRAP FUNCTIONS ###
###########################

# functions to implement bootstrap function to calculate mean, lower and upper CIs assuming normal distribution of means

f_boot = function(vec){
  if(unique(vec)%>%length() == 1){if(is.na(vec%>%unique())){return(NA)}}
  if(anyNA(vec)){vec = c(na.omit(vec))}
  if(unique(vec)%>%length() == 1){return(unique(vec))}
  boot_out = boot::boot.ci(boot.out = one.boot(vec, mean, R=10^2), type = "norm")$t0
  return(boot_out)
}

f_boot_lower = function(vec){
  if(unique(vec)%>%length() == 1){if(is.na(vec%>%unique())){return(NA)}}
  if(anyNA(vec)){vec = c(na.omit(vec))}
  if(unique(vec)%>%length() == 1){return(unique(vec))}
  boot_out_lower = boot::boot.ci(boot.out = one.boot(vec, mean, R=10^2), type = "norm")$normal[2]
  return(boot_out_lower)
}

f_boot_upper = function(vec){
  if(unique(vec)%>%length() == 1){if(is.na(vec%>%unique())){return(NA)}}
  if(anyNA(vec)){vec = c(na.omit(vec))}
  if(unique(vec)%>%length() == 1){return(unique(vec))}
  boot_out_upper = boot::boot.ci(boot.out = one.boot(vec, mean, R=10^2), type = "norm")$normal[3]
  return(boot_out_upper)
}

####################################
### NUMBER OF INFECTIONS AVERTED ###
####################################


# NB: confidence intervals for IRRs "total" calculated following http://www.pmean.com/07/CountData.html
# CI for count data: SE = square root of the number of counts, so simply # +/- 1.96 * SE

f_infAverted = function(dataSurveillanceOutput){
  
  # quantiles, total (assuming count data from a Poisson distribution) and single mean
  dataSurveillanceOutputInfAverted = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(infAverted_025 = quantile(infAverted, 0.025),
                     infAverted_250 = quantile(infAverted, 0.25),
                     infAverted_500 = quantile(infAverted, 0.5),
                     infAverted_750 = quantile(infAverted, 0.75),
                     infAverted_975 = quantile(infAverted, 0.975),
                     infAverted_total = sum(infAverted)/n(),
                     infAverted_total_SD = sqrt(abs(mean(infAverted))),
                     infAverted_singlemean = mean(infAverted),
                     infAverted_singlemean_SD = sd(infAverted))%>%
    mutate(infAverted_total_upper = infAverted_total+1.96*infAverted_total_SD,
           infAverted_total_lower = infAverted_total-1.96*infAverted_total_SD,
           infAverted_singlemean_upper = infAverted_singlemean+1.96*infAverted_singlemean_SD,
           infAverted_singlemean_lower = infAverted_singlemean-1.96*infAverted_singlemean_SD)
  
  # distribution of means
  dataSurveillanceOutputInfAverted_Mean = dataSurveillanceOutput%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(infAvertedMean = mean(infAverted))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(infAverted_distrmean = mean(infAvertedMean),
           infAverted_distrmean_SEM = sd(infAvertedMean))%>%
    mutate(infAverted_distrmean_upper = infAverted_distrmean+1.96*infAverted_distrmean_SEM,
           infAverted_distrmean_lower = infAverted_distrmean-1.96*infAverted_distrmean_SEM)
  
  # bootstrap
  dataSurveillanceOutputInfAverted_bootstrap = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(infAverted_bootstrap = f_boot(infAverted),
              infAverted_bootstrap_lower = f_boot_lower(infAverted),
              infAverted_bootstrap_upper = f_boot_upper(infAverted))
  
  dataSurveillanceOutputInfAverted_Final = left_join(dataSurveillanceOutputInfAverted, 
                                                     left_join(dataSurveillanceOutputInfAverted_Mean, dataSurveillanceOutputInfAverted_bootstrap,
                                                     by = c("stratSurveillance", "typeDailyScreening", "targetDailyScreening", "sensitivity_function")))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputInfAverted_Final)
}

# # TEST
# f_infAverted(test1)
# f_infAverted(test2)

#############################################################
### NUMBER OF INFECTIONS AVERTED RELATIVE TO BASELINE PCR ###
#############################################################

### BASELINE = results from PCR testing
### RELATIVE = results from screening relative to baseline

# NB: as above, confidence intervals for IRR totals calculated following http://www.pmean.com/07/CountData.html
# CI for count data: SE = square root of the number of counts, so simply # +/- 1.96 * SE

f_infAverted_relative = function(dataSurveillanceOutput){
  
  defaultTargetDailyScreening = first(dataSurveillanceOutput$targetDailyScreening)
  defaultTypeDailyScreening = first(dataSurveillanceOutput$typeDailyScreening)
  
  dataSurveillanceOutput_AdmSym = dplyr::filter(dataSurveillanceOutput, stratSurveillance == 'adm_sym',
                                                targetDailyScreening == defaultTargetDailyScreening,
                                                typeDailyScreening == defaultTypeDailyScreening)%>%
    mutate(infAvertedBaseline = infAverted)%>%
    dplyr::select(simulationCTC, testingRound, sensitivity_function, infAvertedBaseline)
  
  dataSurveillanceOutputRelative = left_join(filter(dataSurveillanceOutput, stratSurveillance != 'adm_sym'), dataSurveillanceOutput_AdmSym,
                                             by = c("simulationCTC", "testingRound", "sensitivity_function"))
  
  # quantiles, total (assuming count data from a Poisson distribution) and single mean
  dataSurveillanceOutputInfAvertedRelative = dataSurveillanceOutputRelative%>%
    mutate(infAvertedRelative = infAverted - infAvertedBaseline)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(infAvertedRelative_025 = quantile(infAvertedRelative, 0.025),
                     infAvertedRelative_250 = quantile(infAvertedRelative, 0.25),
                     infAvertedRelative_500 = quantile(infAvertedRelative, 0.5),
                     infAvertedRelative_750 = quantile(infAvertedRelative, 0.75),
                     infAvertedRelative_975 = quantile(infAvertedRelative, 0.975),
                     infAvertedRelative_total = sum(infAvertedRelative)/n(),
                     infAvertedRelative_total_SD = sqrt(abs(mean(infAvertedRelative))),
                     infAvertedRelative_singlemean = mean(infAvertedRelative),
                     infAvertedRelative_singlemean_SD = sd(infAvertedRelative))%>%
    mutate(infAvertedRelative_total_upper = infAvertedRelative_total+1.96*infAvertedRelative_total_SD,
           infAvertedRelative_total_lower = infAvertedRelative_total-1.96*infAvertedRelative_total_SD,
           infAvertedRelative_singlemean_upper = infAvertedRelative_singlemean+1.96*infAvertedRelative_singlemean_SD,
           infAvertedRelative_singlemean_lower = infAvertedRelative_singlemean-1.96*infAvertedRelative_singlemean_SD)
  
  # distribution of means
  dataSurveillanceOutputInfAvertedRelative_Mean = dataSurveillanceOutputRelative%>%
    mutate(infAvertedRelative = infAverted - infAvertedBaseline)%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(infAvertedRelativeMean = mean(infAvertedRelative))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(infAvertedRelative_distrmean = mean(infAvertedRelativeMean),
                     infAvertedRelative_distrmean_SEM = sd(infAvertedRelativeMean))%>%
    mutate(infAvertedRelative_distrmean_upper = infAvertedRelative_distrmean+1.96*infAvertedRelative_distrmean_SEM,
           infAvertedRelative_distrmean_lower = infAvertedRelative_distrmean-1.96*infAvertedRelative_distrmean_SEM)
  
  # bootstrap
  dataSurveillanceOutputInfAvertedRelative_bootstrap = dataSurveillanceOutputRelative%>%
    mutate(infAvertedRelative = infAverted - infAvertedBaseline)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(infAvertedRelative_bootstrap = f_boot(infAvertedRelative),
              infAvertedRelative_bootstrap_lower = f_boot_lower(infAvertedRelative),
              infAvertedRelative_bootstrap_upper = f_boot_upper(infAvertedRelative))
  
  dataSurveillanceOutputInfAvertedRelative_Final = left_join(dataSurveillanceOutputInfAvertedRelative, 
                                                             left_join(dataSurveillanceOutputInfAvertedRelative_Mean, dataSurveillanceOutputInfAvertedRelative_bootstrap,
                                                             by = c("stratSurveillance", "typeDailyScreening", "targetDailyScreening", "sensitivity_function")))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputInfAvertedRelative_Final)
}

# # TEST
# f_infAverted_relative(test1)
# f_infAverted_relative(test2)

###########
### IRR ###
###########

# Incidence rates for included testing strategies (relative to no testing)
# NB: for 'total' this takes sum of incidence across ALL simulations as one outcome,
# and sum of incidence across all simulations with intervention as alternative 'treatment' outcome
# for this outcome, confidence intervals for IRRs calculated following Rothman KJ, Greenland S. Modern epidemiology. 3rd ed. Philadelphia: Lippincott Williams & Wilkins, 2008.

# IRR = (A1/T1)/(A0/T0)
# SD(ln(IR)) = (1/A1 + 1/A0)^0.5
# 95%CI = exp(ln(IR) +- 1.96*SD)

f_irr = function(dataSurveillanceOutput){
  
  # # calculated separately (not actually needed -- say personDays for simulations with and without intervention)
  # personDays = 5916
  
  # Quantiles across all simulations simultaneously
  dataSurveillanceOutputIRRoverall = dataSurveillanceOutput%>%
    mutate(incAfterTesting = incNosocomial - infAverted)%>%
    mutate(incidenceRate = ifelse(incNosocomial == 0, 0, incAfterTesting / incNosocomial))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRateOverall_025 = quantile(incidenceRate, 0.025),
                     incidenceRateOverall_250 = quantile(incidenceRate, 0.25),
                     incidenceRateOverall_500 = quantile(incidenceRate, 0.5),
                     incidenceRateOverall_750 = quantile(incidenceRate, 0.75),
                     incidenceRateOverall_975 = quantile(incidenceRate, 0.975))
  
  # Quantiles across all simulations simultaneously (outbreaks)
  dataSurveillanceOutputIRRoverallOutbreaks = dataSurveillanceOutput%>%
    filter(incNosocomial > 0)%>%
    mutate(incAfterTesting = incNosocomial - infAverted)%>%
    mutate(incidenceRate = incAfterTesting / incNosocomial)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRateOverallOutbreaks_025 = quantile(incidenceRate, 0.025),
                     incidenceRateOverallOutbreaks_250 = quantile(incidenceRate, 0.25),
                     incidenceRateOverallOutbreaks_500 = quantile(incidenceRate, 0.5),
                     incidenceRateOverallOutbreaks_750 = quantile(incidenceRate, 0.75),
                     incidenceRateOverallOutbreaks_975 = quantile(incidenceRate, 0.975))
  
  # overall total as one simulation
  dataSurveillanceOutputIRR_total = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incNosocomial_sum = sum(incNosocomial),
                     infAverted_sum = sum(infAverted))%>%
    mutate(incAfterTesting_sum = incNosocomial_sum - infAverted_sum)%>%
    mutate(incidenceRate_total = incAfterTesting_sum / incNosocomial_sum)%>%
    mutate(incidenceRate_total_SD = (1/incAfterTesting_sum + 1/incNosocomial_sum)^0.5)%>%
    mutate(incidenceRate_total_upper = exp(log(incidenceRate_total) + 1.96*incidenceRate_total_SD))%>%
    mutate(incidenceRate_total_lower = exp(log(incidenceRate_total) - 1.96*incidenceRate_total_SD))
  
  # single means (ONLY INCLUDE OUTBREAKS)
  dataSurveillanceOutputIRR_singlemean = dataSurveillanceOutput%>%
    filter(incNosocomial>0)%>%
    mutate(incidenceRate = 1- (infAverted/incNosocomial))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRate_singlemean = mean(incidenceRate),
                     incidenceRate_singlemean_SD = sd(incidenceRate))%>%
    mutate(incidenceRate_singlemean_upper = incidenceRate_singlemean + 1.96*incidenceRate_singlemean_SD)%>%
    mutate(incidenceRate_singlemean_lower = incidenceRate_singlemean - 1.96*incidenceRate_singlemean_SD)
  
  # distribution of means (ONLY INCLUDE OUTBREAKS)
  dataSurveillanceOutputIRR_distrmean = dataSurveillanceOutput%>%
    filter(incNosocomial>0)%>%
    mutate(incidenceRate = 1- (infAverted/incNosocomial))%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRateMean = mean(incidenceRate))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRate_distrmean = mean(incidenceRateMean),
                     incidenceRate_distrmean_SEM = sd(incidenceRateMean))%>%
    mutate(incidenceRate_distrmean_upper = incidenceRate_distrmean + 1.96*incidenceRate_distrmean_SEM)%>%
    mutate(incidenceRate_distrmean_lower = incidenceRate_distrmean - 1.96*incidenceRate_distrmean_SEM)
  
  # bootstrap (ONLY INCLUDE OUTBREAKS)
  dataSurveillanceOutputIRR_bootstrap = dataSurveillanceOutput%>%
    filter(incNosocomial>0)%>%
    mutate(incidenceRate = 1- (infAverted/incNosocomial))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(incidenceRate_bootstrap = f_boot(incidenceRate),
              incidenceRate_bootstrap_lower = f_boot_lower(incidenceRate),
              incidenceRate_bootstrap_upper = f_boot_upper(incidenceRate))
  
  
  
  dataSurveillanceOutputIRR = left_join(dataSurveillanceOutputIRRoverall,
                                        left_join(dataSurveillanceOutputIRRoverallOutbreaks,
                                                  left_join(dataSurveillanceOutputIRR_total, 
                                                            left_join(dataSurveillanceOutputIRR_singlemean, 
                                                                      left_join(dataSurveillanceOutputIRR_distrmean, dataSurveillanceOutputIRR_bootstrap)))))%>%
    as.data.frame()
  
  
  return(dataSurveillanceOutputIRR)
  
}

# # TEST
# f_irr(test1)
# f_irr(test2)



###############################
### IRR RELATIVE TO ADM/SYM ###
###############################

### BASELINE = results from PCR testing
### RELATIVE = results from screening relative to baseline (only include simulations with )

# NB: if PCR detects all cases, assume IRR for screening = 1 (no additional benefit to screening)

f_irr_relative = function(dataSurveillanceOutput){
  
  defaultTargetDailyScreening = first(dataSurveillanceOutput$targetDailyScreening)
  defaultTypeDailyScreening = first(dataSurveillanceOutput$typeDailyScreening)
  
  dataSurveillanceOutput_AdmSym = dplyr::filter(dataSurveillanceOutput, stratSurveillance == 'adm_sym',
                                                targetDailyScreening == defaultTargetDailyScreening,
                                                typeDailyScreening == defaultTypeDailyScreening)%>%
    mutate(incNosocomialBaseline = incNosocomial - infAverted)%>%
    dplyr::select(simulationCTC, testingRound, sensitivity_function, incNosocomialBaseline)
  
  dataSurveillanceOutputRelative = left_join(filter(dataSurveillanceOutput, stratSurveillance != 'adm_sym'), dataSurveillanceOutput_AdmSym, 
                                             by = c("simulationCTC", "testingRound", "sensitivity_function"))
  
  
  # Quantiles across all simulations simultaneously
  dataSurveillanceOutputIRRoverall = dataSurveillanceOutputRelative%>%
    mutate(incAfterScreening = incNosocomial - infAverted)%>%
    mutate(incidenceRateRelative = ifelse(incNosocomialBaseline == 0, 1, incAfterScreening / incNosocomialBaseline))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRateOverall_025 = quantile(incidenceRateRelative, 0.025),
                     incidenceRateOverall_250 = quantile(incidenceRateRelative, 0.25),
                     incidenceRateOverall_500 = quantile(incidenceRateRelative, 0.5),
                     incidenceRateOverall_750 = quantile(incidenceRateRelative, 0.75),
                     incidenceRateOverall_975 = quantile(incidenceRateRelative, 0.975))
  
  # Quantiles across all simulations simultaneously (outbreaks [after PCR!])
  dataSurveillanceOutputIRRoverallOutbreaks = dataSurveillanceOutputRelative%>%
    filter(incNosocomialBaseline > 0)%>%
    mutate(incAfterScreening = incNosocomial - infAverted)%>%
    mutate(incidenceRateRelative = ifelse(incNosocomialBaseline == 0, 1, incAfterScreening / incNosocomialBaseline))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRateOverallOutbreaks_025 = quantile(incidenceRateRelative, 0.025),
                     incidenceRateOverallOutbreaks_250 = quantile(incidenceRateRelative, 0.25),
                     incidenceRateOverallOutbreaks_500 = quantile(incidenceRateRelative, 0.5),
                     incidenceRateOverallOutbreaks_750 = quantile(incidenceRateRelative, 0.75),
                     incidenceRateOverallOutbreaks_975 = quantile(incidenceRateRelative, 0.975))
  
  
  
  # overall total as one simulation
  dataSurveillanceOutputIRR_total = dataSurveillanceOutputRelative%>%
    mutate(incAfterScreening = incNosocomial - infAverted)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incNosocomialBaseline_sum = sum(incNosocomialBaseline),
                     incAfterScreening_sum = sum(incAfterScreening))%>%
    mutate(incidenceRate_total = incAfterScreening_sum / incNosocomialBaseline_sum)%>%
    mutate(incidenceRate_total_SD = (1/incAfterScreening_sum + 1/incNosocomialBaseline_sum)^0.5)%>%
    mutate(incidenceRate_total_upper = exp(log(incidenceRate_total) + 1.96*incidenceRate_total_SD))%>%
    mutate(incidenceRate_total_lower = exp(log(incidenceRate_total) - 1.96*incidenceRate_total_SD))
  
  # single means (ONLY INCLUDE OUTBREAKS, BUT COUNT IRR = 1 WHEN PCR ALREADY ELIMINATES ALL CASES)
  dataSurveillanceOutputIRR_singlemean = dataSurveillanceOutputRelative%>%
    filter(incNosocomialBaseline > 0)%>%
    mutate(incAfterScreening = incNosocomial - infAverted,
           incidenceRate = ifelse(incNosocomialBaseline == 0, 1, incAfterScreening/incNosocomialBaseline))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRate_singlemean = mean(incidenceRate),
                     incidenceRate_singlemean_SD = sd(incidenceRate))%>%
    mutate(incidenceRate_singlemean_upper = incidenceRate_singlemean + 1.96*incidenceRate_singlemean_SD,
           incidenceRate_singlemean_lower = incidenceRate_singlemean - 1.96*incidenceRate_singlemean_SD)
  
  # distribution of means (ONLY INCLUDE OUTBREAKS, BUT COUNT IRR = 1 WHEN PCR ALREADY ELIMINATES ALL CASES)
  dataSurveillanceOutputIRR_distrmean = dataSurveillanceOutputRelative%>%
    filter(incNosocomialBaseline > 0)%>%
    mutate(incAfterScreening = incNosocomial - infAverted,
           incidenceRate = ifelse(incNosocomialBaseline == 0, 1, incAfterScreening/incNosocomialBaseline))%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRateMean = mean(incidenceRate))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incidenceRate_distrmean = mean(incidenceRateMean),
                     incidenceRate_distrmean_SEM = sd(incidenceRateMean))%>%
    mutate(incidenceRate_distrmean_upper = incidenceRate_distrmean + 1.96*incidenceRate_distrmean_SEM,
           incidenceRate_distrmean_lower = incidenceRate_distrmean - 1.96*incidenceRate_distrmean_SEM)
  
  # bootstrap (ONLY INCLUDE OUTBREAKS)
  dataSurveillanceOutputIRR_bootstrap = dataSurveillanceOutputRelative%>%
    filter(incNosocomialBaseline>0)%>%
    mutate(incAfterScreening = incNosocomial - infAverted,
           incidenceRate = ifelse(incNosocomialBaseline == 0, 1, incAfterScreening/incNosocomialBaseline))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(incidenceRate_bootstrap = f_boot(incidenceRate),
              incidenceRate_bootstrap_lower = f_boot_lower(incidenceRate),
              incidenceRate_bootstrap_upper = f_boot_upper(incidenceRate))
  
  
  dataSurveillanceOutputIRR = left_join(dataSurveillanceOutputIRRoverall, 
                                        left_join(dataSurveillanceOutputIRRoverallOutbreaks,
                                                  left_join(dataSurveillanceOutputIRR_total,
                                                            left_join(dataSurveillanceOutputIRR_singlemean, 
                                                                      left_join(dataSurveillanceOutputIRR_distrmean, dataSurveillanceOutputIRR_bootstrap)))))%>%
    as.data.frame()
  
  
  return(dataSurveillanceOutputIRR)
  
}

# # TEST
# f_irr_relative(test1)
# f_irr_relative(test2)


########################################
### Proportion of infections averted ###
########################################

f_propAverted = function(dataSurveillanceOutput){
  
  # Quantiles across all simulations simultaneously
  dataSurveillanceOutputPropAvertedOverall = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAvertedOverall_025 = quantile(propAverted, 0.025),
                     propAvertedOverall_250 = quantile(propAverted, 0.25),
                     propAvertedOverall_500 = quantile(propAverted, 0.5),
                     propAvertedOverall_750 = quantile(propAverted, 0.75),
                     propAvertedOverall_975 = quantile(propAverted, 0.975))
  
  # Quantiles across all simulations simultaneously (outbreaks only)
  dataSurveillanceOutputPropAvertedOverallOutbreaks = dataSurveillanceOutput%>%
    dplyr::filter(incNosocomial > 0)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAvertedOverallOutbreaks_025 = quantile(propAverted, 0.025),
                     propAvertedOverallOutbreaks_250 = quantile(propAverted, 0.25),
                     propAvertedOverallOutbreaks_500 = quantile(propAverted, 0.5),
                     propAvertedOverallOutbreaks_750 = quantile(propAverted, 0.75),
                     propAvertedOverallOutbreaks_975 = quantile(propAverted, 0.975),
                     numOutbreaks = length(unique(simulationCTC)))
  
  # overall total as one simulation
  dataSurveillanceOutputPropAverted_total = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incNosocomial_sum = sum(incNosocomial),
                     infAverted_sum = sum(infAverted))%>%
    mutate(propAverted_total = infAverted_sum / incNosocomial_sum)%>%
    # calculate 95% CIs of proportion using formula p = p +- 1.96 * sqrt((p*(1-p))/n), where n = number of cases across 100 CTC simulations (NOT # of simulations as assumed previously)
    mutate(propAverted_total_SD = sqrt((propAverted_total*(1-propAverted_total))/(incNosocomial_sum)))%>%
    mutate(propAverted_total_lower = propAverted_total - 1.96*propAverted_total_SD)%>%
    mutate(propAverted_total_upper = propAverted_total + 1.96*propAverted_total_SD)
  
  # single means (ONLY INCLUDE OUTBREAKS)
  dataSurveillanceOutputPropAverted_singlemean = dataSurveillanceOutput%>%
    filter(incNosocomial > 0)%>%
    mutate(propAverted = infAverted/incNosocomial)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAverted_singlemean = mean(propAverted),
                     propAverted_singlemean_SD = sd(propAverted))%>%
    mutate(propAverted_singlemean_upper = propAverted_singlemean + 1.96 * propAverted_singlemean_SD,
           propAverted_singlemean_lower = propAverted_singlemean - 1.96 * propAverted_singlemean_SD)

  # distribution of means (ONLY INCLUDE OUTBREAKS)
  dataSurveillanceOutputPropAverted_distrmean = dataSurveillanceOutput%>%
    filter(incNosocomial > 0)%>%
    mutate(propAverted = infAverted/incNosocomial)%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAvertedMean = mean(propAverted))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAverted_distrmean = mean(propAvertedMean),
                     propAverted_distrmean_SEM = sd(propAvertedMean))%>%
    mutate(propAverted_distrmean_upper = propAverted_distrmean + 1.96 * propAverted_distrmean_SEM,
           propAverted_distrmean_lower = propAverted_distrmean - 1.96 * propAverted_distrmean_SEM)
  
  # bootstrap (ONLY INCLUDE OUTBREAKS)
  dataSurveillanceOutputPropAverted_bootstrap = dataSurveillanceOutput%>%
    filter(incNosocomial>0)%>%
    mutate(propAverted = infAverted/incNosocomial)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(propAverted_bootstrap = f_boot(propAverted),
              propAverted_bootstrap_lower = f_boot_lower(propAverted),
              propAverted_bootstrap_upper = f_boot_upper(propAverted))


  dataSurveillanceOutputPropAverted = left_join(dataSurveillanceOutputPropAvertedOverall,
                                                left_join(dataSurveillanceOutputPropAvertedOverallOutbreaks,
                                                          left_join(dataSurveillanceOutputPropAverted_total, 
                                                                    left_join(dataSurveillanceOutputPropAverted_singlemean, 
                                                                              left_join(dataSurveillanceOutputPropAverted_distrmean, dataSurveillanceOutputPropAverted_bootstrap)))))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputPropAverted)
}

# # TEST
# f_propAverted(test1)
# f_propAverted(test2)


############################################################
### Proportion of infections averted RELATIVE to adm/sym ###
############################################################

f_propAverted_relative = function(dataSurveillanceOutput){
  
  defaultTargetDailyScreening = first(dataSurveillanceOutput$targetDailyScreening)
  defaultTypeDailyScreening = first(dataSurveillanceOutput$typeDailyScreening)
  
  dataSurveillanceOutput_AdmSym = dplyr::filter(dataSurveillanceOutput, stratSurveillance == 'adm_sym',
                                                targetDailyScreening == defaultTargetDailyScreening,
                                                typeDailyScreening == defaultTypeDailyScreening)%>%
    mutate(infAvertedBaseline = infAverted,
           incNosocomialAfterBaseline = incNosocomial - infAverted)%>%
    dplyr::select(simulationCTC, testingRound, incNosocomialAfterBaseline, infAvertedBaseline, sensitivity_function)
  
  dataSurveillanceOutputRelative = left_join(filter(dataSurveillanceOutput, stratSurveillance != 'adm_sym'), dataSurveillanceOutput_AdmSym)%>%
    mutate(infAvertedAfterBaseline = infAverted - infAvertedBaseline)
  
  # Quantiles across all simulations simultaneously
  dataSurveillanceOutputPropAvertedOverall = dataSurveillanceOutputRelative%>%
    mutate(infAvertedScreening = infAverted - infAvertedBaseline,
           incNosocomialAfterBaseline = incNosocomial - infAvertedBaseline,
           propAvertedOverall = infAvertedScreening/incNosocomialAfterBaseline)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAvertedOverall_025 = quantile(propAvertedOverall, 0.025, na.rm = T),
                     propAvertedOverall_250 = quantile(propAvertedOverall, 0.25, na.rm = T),
                     propAvertedOverall_500 = quantile(propAvertedOverall, 0.5, na.rm = T),
                     propAvertedOverall_750 = quantile(propAvertedOverall, 0.75, na.rm = T),
                     propAvertedOverall_975 = quantile(propAvertedOverall, 0.975, na.rm = T))
  
  # Quantiles across all simulations simultaneously (outbreaks only)
  dataSurveillanceOutputPropAvertedOutbreaks = dataSurveillanceOutputRelative%>%
    mutate(infAvertedScreening = infAverted - infAvertedBaseline,
           incNosocomialAfterBaseline = incNosocomial - infAvertedBaseline)%>%
    filter(incNosocomialAfterBaseline>0)%>%
    mutate(propAvertedOutbreaks = infAvertedScreening/incNosocomialAfterBaseline)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAvertedOverallOutbreaks_025 = quantile(propAvertedOutbreaks, 0.025),
                     propAvertedOverallOutbreaks_250 = quantile(propAvertedOutbreaks, 0.25),
                     propAvertedOverallOutbreaks_500 = quantile(propAvertedOutbreaks, 0.5),
                     propAvertedOverallOutbreaks_750 = quantile(propAvertedOutbreaks, 0.75),
                     propAvertedOverallOutbreaks_975 = quantile(propAvertedOutbreaks, 0.975))
  
  
  # overall total as one simulation
  dataSurveillanceOutputPropAvertedTotal = dataSurveillanceOutputRelative%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(incNosocomialAfterBaseline_sum = sum(incNosocomialAfterBaseline),
                     infAvertedAfterBaseline_sum = sum(infAvertedAfterBaseline))%>%
    mutate(propAverted_total = infAvertedAfterBaseline_sum / incNosocomialAfterBaseline_sum)%>%
    # calculate 95% CIs of proportion using formula p = p +- 1.96 * sqrt((p*(1-p))/n), where n = number of cases across 100 CTC simulations (NOT # of simulations as assumed previously)
    mutate(propAverted_total_SD = sqrt((abs(propAverted_total)*(1-abs(propAverted_total)))/(incNosocomialAfterBaseline_sum)))%>%
    mutate(propAverted_total_lower = propAverted_total - 1.96*propAverted_total_SD)%>%
    mutate(propAverted_total_upper = propAverted_total + 1.96*propAverted_total_SD)
  
  
  
  # single means (ONLY INCLUDE OUTBREAKS, EVEN AFTER BASELINE TESTING)
  dataSurveillanceOutputPropAverted_singlemean = dataSurveillanceOutputRelative%>%
    filter(incNosocomialAfterBaseline > 0)%>%
    mutate(propAverted = infAvertedAfterBaseline/incNosocomialAfterBaseline)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAverted_singlemean = mean(propAverted),
                     propAverted_singlemean_SD = sd(propAverted))%>%
    mutate(propAverted_singlemean_upper = propAverted_singlemean + 1.96 * propAverted_singlemean_SD,
           propAverted_singlemean_lower = propAverted_singlemean - 1.96 * propAverted_singlemean_SD)
  
  # distribution of means (ONLY INCLUDE OUTBREAKS, EVEN AFTER BASELINE TESTING)
  dataSurveillanceOutputPropAverted_distrmean = dataSurveillanceOutputRelative%>%
    filter(incNosocomialAfterBaseline > 0)%>%
    mutate(propAverted = infAvertedAfterBaseline/incNosocomialAfterBaseline)%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAvertedMean = mean(propAverted))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(propAverted_distrmean = mean(propAvertedMean),
                     propAverted_distrmean_SEM = sd(propAvertedMean))%>%
    mutate(propAverted_distrmean_upper = propAverted_distrmean + 1.96 * propAverted_distrmean_SEM,
           propAverted_distrmean_lower = propAverted_distrmean - 1.96 * propAverted_distrmean_SEM)
  
  # bootstrap (ONLY INCLUDE OUTBREAKS)
  dataSurveillanceOutputPropAverted_bootstrap = dataSurveillanceOutputRelative%>%
    filter(incNosocomialAfterBaseline>0)%>%
    mutate(propAverted = infAvertedAfterBaseline/incNosocomialAfterBaseline)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(propAverted_bootstrap = f_boot(propAverted),
              propAverted_bootstrap_lower = f_boot_lower(propAverted),
              propAverted_bootstrap_upper = f_boot_upper(propAverted))
  
  
  dataSurveillanceOutputPropAverted = left_join(dataSurveillanceOutputPropAvertedOverall,
                                                left_join(dataSurveillanceOutputPropAvertedOutbreaks, 
                                                          left_join(dataSurveillanceOutputPropAvertedTotal,
                                                                    left_join(dataSurveillanceOutputPropAverted_singlemean, 
                                                                              left_join(dataSurveillanceOutputPropAverted_distrmean, dataSurveillanceOutputPropAverted_bootstrap)))))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputPropAverted)
}

# # TEST
# f_propAverted_relative(test1)
# f_propAverted_relative(test2)


################################################
### Proportion of simulations with outbreaks ###
################################################

f_probOutbreaks = function(dataSurveillanceOutput){
  
  dataSurveillanceOutputProbOutbreaksPrep = dataSurveillanceOutput%>%
    mutate(incAfterTesting = incNosocomial - infAverted)%>%
    mutate(outbreakAfterTesting_1 = ifelse(incAfterTesting >= 1, 1, 0))%>%
    mutate(outbreakAfterTesting_2 = ifelse(incAfterTesting >= 2, 1, 0))%>%
    mutate(outbreakAfterTesting_3 = ifelse(incAfterTesting >= 3, 1, 0))%>%
    mutate(outbreakAfterTesting_4 = ifelse(incAfterTesting >= 4, 1, 0))%>%
    mutate(outbreakAfterTesting_5 = ifelse(incAfterTesting >= 5, 1, 0))%>%
    mutate(outbreakAfterTesting_10 = ifelse(incAfterTesting >= 10, 1, 0))%>%
    mutate(outbreakAfterTesting_20 = ifelse(incAfterTesting >= 20, 1, 0))
  
  dataSurveillanceOutputProbOutbreaksIndividual = dataSurveillanceOutputProbOutbreaksPrep%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening)%>%
    dplyr::summarise(probOutbreakAfterTesting_1 = sum(outbreakAfterTesting_1)/length(outbreakAfterTesting_1),
                     probOutbreakAfterTesting_2 = sum(outbreakAfterTesting_2)/length(outbreakAfterTesting_2),
                     probOutbreakAfterTesting_3 = sum(outbreakAfterTesting_3)/length(outbreakAfterTesting_3),
                     probOutbreakAfterTesting_4 = sum(outbreakAfterTesting_4)/length(outbreakAfterTesting_4),
                     probOutbreakAfterTesting_5 = sum(outbreakAfterTesting_5)/length(outbreakAfterTesting_5),
                     probOutbreakAfterTesting_10 = sum(outbreakAfterTesting_10)/length(outbreakAfterTesting_10),
                     probOutbreakAfterTesting_20 = sum(outbreakAfterTesting_20)/length(outbreakAfterTesting_20))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening)%>%
    dplyr::summarise(probOutbreakAfterTesting_1_025 = quantile(probOutbreakAfterTesting_1, 0.025),
                     probOutbreakAfterTesting_1_250 = quantile(probOutbreakAfterTesting_1, 0.25),
                     probOutbreakAfterTesting_1_500 = quantile(probOutbreakAfterTesting_1, 0.5),
                     probOutbreakAfterTesting_1_750 = quantile(probOutbreakAfterTesting_1, 0.75),
                     probOutbreakAfterTesting_1_975 = quantile(probOutbreakAfterTesting_1, 0.975),
                     probOutbreakAfterTesting_2_025 = quantile(probOutbreakAfterTesting_2, 0.025),
                     probOutbreakAfterTesting_2_250 = quantile(probOutbreakAfterTesting_2, 0.25),
                     probOutbreakAfterTesting_2_500 = quantile(probOutbreakAfterTesting_2, 0.5),
                     probOutbreakAfterTesting_2_750 = quantile(probOutbreakAfterTesting_2, 0.75),
                     probOutbreakAfterTesting_2_975 = quantile(probOutbreakAfterTesting_2, 0.975),
                     probOutbreakAfterTesting_3_025 = quantile(probOutbreakAfterTesting_3, 0.025),
                     probOutbreakAfterTesting_3_250 = quantile(probOutbreakAfterTesting_3, 0.25),
                     probOutbreakAfterTesting_3_500 = quantile(probOutbreakAfterTesting_3, 0.5),
                     probOutbreakAfterTesting_3_750 = quantile(probOutbreakAfterTesting_3, 0.75),
                     probOutbreakAfterTesting_3_975 = quantile(probOutbreakAfterTesting_3, 0.975),
                     probOutbreakAfterTesting_4_025 = quantile(probOutbreakAfterTesting_4, 0.025),
                     probOutbreakAfterTesting_4_250 = quantile(probOutbreakAfterTesting_4, 0.25),
                     probOutbreakAfterTesting_4_500 = quantile(probOutbreakAfterTesting_4, 0.5),
                     probOutbreakAfterTesting_4_750 = quantile(probOutbreakAfterTesting_4, 0.75),
                     probOutbreakAfterTesting_4_975 = quantile(probOutbreakAfterTesting_4, 0.975),
                     probOutbreakAfterTesting_5_025 = quantile(probOutbreakAfterTesting_5, 0.025),
                     probOutbreakAfterTesting_5_250 = quantile(probOutbreakAfterTesting_5, 0.25),
                     probOutbreakAfterTesting_5_500 = quantile(probOutbreakAfterTesting_5, 0.5),
                     probOutbreakAfterTesting_5_750 = quantile(probOutbreakAfterTesting_5, 0.75),
                     probOutbreakAfterTesting_5_975 = quantile(probOutbreakAfterTesting_5, 0.975),
                     probOutbreakAfterTesting_10_025 = quantile(probOutbreakAfterTesting_10, 0.025),
                     probOutbreakAfterTesting_10_250 = quantile(probOutbreakAfterTesting_10, 0.25),
                     probOutbreakAfterTesting_10_500 = quantile(probOutbreakAfterTesting_10, 0.5),
                     probOutbreakAfterTesting_10_750 = quantile(probOutbreakAfterTesting_10, 0.75),
                     probOutbreakAfterTesting_10_975 = quantile(probOutbreakAfterTesting_10, 0.975),
                     probOutbreakAfterTesting_20_025 = quantile(probOutbreakAfterTesting_20, 0.025),
                     probOutbreakAfterTesting_20_250 = quantile(probOutbreakAfterTesting_20, 0.25),
                     probOutbreakAfterTesting_20_500 = quantile(probOutbreakAfterTesting_20, 0.5),
                     probOutbreakAfterTesting_20_750 = quantile(probOutbreakAfterTesting_20, 0.75),
                     probOutbreakAfterTesting_20_975 = quantile(probOutbreakAfterTesting_20, 0.975))
  
  dataSurveillanceOutputProbOutbreaksOverall = dataSurveillanceOutputProbOutbreaksPrep%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening)%>%
    dplyr::summarise(probOutbreakAfterTesting_1 = sum(outbreakAfterTesting_1)/length(outbreakAfterTesting_1),
                     probOutbreakAfterTesting_2 = sum(outbreakAfterTesting_2)/length(outbreakAfterTesting_2),
                     probOutbreakAfterTesting_3 = sum(outbreakAfterTesting_3)/length(outbreakAfterTesting_3),
                     probOutbreakAfterTesting_4 = sum(outbreakAfterTesting_4)/length(outbreakAfterTesting_4),
                     probOutbreakAfterTesting_5 = sum(outbreakAfterTesting_5)/length(outbreakAfterTesting_5),
                     probOutbreakAfterTesting_10 = sum(outbreakAfterTesting_10)/length(outbreakAfterTesting_10),
                     probOutbreakAfterTesting_20 = sum(outbreakAfterTesting_20)/length(outbreakAfterTesting_20))%>%
    # calculate 95% CIs of proportion using formula p = p +- 1.96 * sqrt((p*(1-p))/n), where n = number of simulations
    mutate(probOutbreakAfterTesting_1_lower = probOutbreakAfterTesting_1 - 1.96*sqrt((probOutbreakAfterTesting_1*(1-probOutbreakAfterTesting_1))/100))%>%
    mutate(probOutbreakAfterTesting_1_upper = probOutbreakAfterTesting_1 + 1.96*sqrt((probOutbreakAfterTesting_1*(1-probOutbreakAfterTesting_1))/100))%>%
    mutate(probOutbreakAfterTesting_2_lower = probOutbreakAfterTesting_2 - 1.96*sqrt((probOutbreakAfterTesting_2*(1-probOutbreakAfterTesting_2))/100))%>%
    mutate(probOutbreakAfterTesting_2_upper = probOutbreakAfterTesting_2 + 1.96*sqrt((probOutbreakAfterTesting_2*(1-probOutbreakAfterTesting_2))/100))%>%
    mutate(probOutbreakAfterTesting_3_lower = probOutbreakAfterTesting_3 - 1.96*sqrt((probOutbreakAfterTesting_3*(1-probOutbreakAfterTesting_3))/100))%>%
    mutate(probOutbreakAfterTesting_3_upper = probOutbreakAfterTesting_3 + 1.96*sqrt((probOutbreakAfterTesting_3*(1-probOutbreakAfterTesting_3))/100))%>%
    mutate(probOutbreakAfterTesting_4_lower = probOutbreakAfterTesting_4 - 1.96*sqrt((probOutbreakAfterTesting_4*(1-probOutbreakAfterTesting_4))/100))%>%
    mutate(probOutbreakAfterTesting_4_upper = probOutbreakAfterTesting_4 + 1.96*sqrt((probOutbreakAfterTesting_4*(1-probOutbreakAfterTesting_4))/100))%>%
    mutate(probOutbreakAfterTesting_5_lower = probOutbreakAfterTesting_5 - 1.96*sqrt((probOutbreakAfterTesting_5*(1-probOutbreakAfterTesting_5))/100))%>%
    mutate(probOutbreakAfterTesting_5_upper = probOutbreakAfterTesting_5 + 1.96*sqrt((probOutbreakAfterTesting_5*(1-probOutbreakAfterTesting_5))/100))%>%
    mutate(probOutbreakAfterTesting_10_lower = probOutbreakAfterTesting_10 - 1.96*sqrt((probOutbreakAfterTesting_10*(1-probOutbreakAfterTesting_10))/100))%>%
    mutate(probOutbreakAfterTesting_10_upper = probOutbreakAfterTesting_10 + 1.96*sqrt((probOutbreakAfterTesting_10*(1-probOutbreakAfterTesting_10))/100))%>%
    mutate(probOutbreakAfterTesting_20_lower = probOutbreakAfterTesting_20 - 1.96*sqrt((probOutbreakAfterTesting_20*(1-probOutbreakAfterTesting_20))/100))%>%
    mutate(probOutbreakAfterTesting_20_upper = probOutbreakAfterTesting_20 + 1.96*sqrt((probOutbreakAfterTesting_20*(1-probOutbreakAfterTesting_20))/100))
  
  dataSurveillanceOutputProbOutbreaks = merge(dataSurveillanceOutputProbOutbreaksOverall, dataSurveillanceOutputProbOutbreaksIndividual)
  
  return(dataSurveillanceOutputProbOutbreaks)
  
}

# # TEST
# f_probOutbreaks(test1)
# f_probOutbreaks(test2)



###########################################################
### PLOT PREP: Proportion of simulations with outbreaks ###
####################################@######################

# Function to summarize data, characterize outbreaks as:
# (i) any_trans [any transmission], (ii) any_outbreak [any outbreak, 5+ cases], (iii) any_outbreak_large [large outbreak, 20+ cases]

f_probOutbreaks_plotReady = function(dataProbOutbreaks){
  dataProbOutbreaks_plotReady = dplyr::select(dataProbOutbreaks, stratSurveillance, typeDailyScreening, targetDailyScreening,
                                                  probOutbreakAfterTesting_1, probOutbreakAfterTesting_1_lower, probOutbreakAfterTesting_1_upper,
                                                  probOutbreakAfterTesting_1_500, probOutbreakAfterTesting_1_250, probOutbreakAfterTesting_1_750,
                                                  probOutbreakAfterTesting_5, probOutbreakAfterTesting_5_lower, probOutbreakAfterTesting_5_upper,
                                                  probOutbreakAfterTesting_5_500, probOutbreakAfterTesting_5_250, probOutbreakAfterTesting_5_750,
                                                  probOutbreakAfterTesting_20, probOutbreakAfterTesting_20_lower, probOutbreakAfterTesting_20_upper,
                                                  probOutbreakAfterTesting_20_500, probOutbreakAfterTesting_20_250, probOutbreakAfterTesting_20_750)%>%
    pivot_longer(-c(stratSurveillance, typeDailyScreening, targetDailyScreening), names_to = "measure", values_to = "value")%>%
    mutate(method = ifelse(grepl('_500', measure) | grepl('_250', measure) | grepl('_750', measure), 'quantiles', 'mean'))%>%
    mutate(quantile = ifelse(grepl('_lower' , measure) | grepl('_250' , measure), 'lower', ifelse(grepl('_upper', measure) | grepl('_750', measure), 'upper', 'median')))%>%
    mutate(threshold = ifelse(grepl('Testing_20', measure), 'any_outbreak_large', ifelse(grepl('Testing_5', measure), 'any_outbreak', 'any_trans')))%>%
    dplyr::select(-measure)%>%
    pivot_wider(names_from = 'quantile', values_from = 'value')
  
  return(dataProbOutbreaks_plotReady)
}




############################
### TEST POSITIVITY RATE ###
############################


f_testPositivity = function(dataSurveillanceOutput){
  
  dataSurveillanceOutputTestPositivityOverall = dataSurveillanceOutput%>%
    mutate(PCR.positivity_rate = (PCR.n_test_true_positive + PCR.n_test_false_positive)/PCR.n_test,
           RDT.positivity_rate = (RDT.n_test_true_positive + RDT.n_test_false_positive)/RDT.n_test)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.positivity_rate_025 = quantile(PCR.positivity_rate, 0.025, na.rm = T),
                     PCR.positivity_rate_250 = quantile(PCR.positivity_rate, 0.25, na.rm = T),
                     PCR.positivity_rate_500 = quantile(PCR.positivity_rate, 0.5, na.rm = T),
                     PCR.positivity_rate_750 = quantile(PCR.positivity_rate, 0.75, na.rm = T),
                     PCR.positivity_rate_975 = quantile(PCR.positivity_rate, 0.975, na.rm = T),
                     RDT.positivity_rate_025 = quantile(RDT.positivity_rate, 0.025, na.rm = T),
                     RDT.positivity_rate_250 = quantile(RDT.positivity_rate, 0.25, na.rm = T),
                     RDT.positivity_rate_500 = quantile(RDT.positivity_rate, 0.5, na.rm = T),
                     RDT.positivity_rate_750 = quantile(RDT.positivity_rate, 0.75, na.rm = T),
                     RDT.positivity_rate_975 = quantile(RDT.positivity_rate, 0.975, na.rm = T))
  
  
  # overall total as one simulation
  dataSurveillanceOutputTestPositivityTotal = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(PCR.n_test_sum = sum(PCR.n_test),
              PCR.n_positive_sum = sum(PCR.n_test_true_positive)+sum(PCR.n_test_false_positive),
              RDT.n_test_sum = sum(RDT.n_test),
              RDT.n_positive_sum = sum(RDT.n_test_true_positive)+sum(RDT.n_test_false_positive))%>%
    mutate(PCR.positivity_rate_total = PCR.n_positive_sum/PCR.n_test_sum,
           PCR.positivity_rate_total_SD = sqrt((PCR.positivity_rate_total*(1-PCR.positivity_rate_total))/PCR.n_test_sum),
           PCR.positivity_rate_total_lower = PCR.positivity_rate_total - 1.96 * PCR.positivity_rate_total_SD,
           PCR.positivity_rate_total_upper = PCR.positivity_rate_total + 1.96 * PCR.positivity_rate_total_SD,
           RDT.positivity_rate_total = RDT.n_positive_sum/RDT.n_test_sum,
           RDT.positivity_rate_total_SD = sqrt((RDT.positivity_rate_total*(1-RDT.positivity_rate_total))/RDT.n_test_sum),
           RDT.positivity_rate_total_lower = RDT.positivity_rate_total - 1.96 * RDT.positivity_rate_total_SD,
           RDT.positivity_rate_total_upper = RDT.positivity_rate_total + 1.96 * RDT.positivity_rate_total_SD)
  
  # single means
  dataSurveillanceOutputTestPositivity_singlemean = dataSurveillanceOutput%>%
    mutate(PCR.positivity_rate = (PCR.n_test_true_positive + PCR.n_test_false_positive)/PCR.n_test,
           RDT.positivity_rate = (RDT.n_test_true_positive + RDT.n_test_false_positive)/RDT.n_test)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.positivity_rate_singlemean = mean(PCR.positivity_rate),
                     PCR.positivity_rate_singlemean_SD = sd(PCR.positivity_rate),
                     RDT.positivity_rate_singlemean = mean(RDT.positivity_rate),
                     RDT.positivity_rate_singlemean_SD = sd(RDT.positivity_rate))%>%
    mutate(PCR.positivity_rate_singlemean_upper = PCR.positivity_rate_singlemean + 1.96 * PCR.positivity_rate_singlemean_SD,
           PCR.positivity_rate_singlemean_lower = PCR.positivity_rate_singlemean - 1.96 * PCR.positivity_rate_singlemean_SD,
           RDT.positivity_rate_singlemean_upper = RDT.positivity_rate_singlemean + 1.96 * RDT.positivity_rate_singlemean_SD,
           RDT.positivity_rate_singlemean_lower = RDT.positivity_rate_singlemean - 1.96 * RDT.positivity_rate_singlemean_SD)
  
  # distribution of means
  dataSurveillanceOutputTestPositivity_distrmean = dataSurveillanceOutput%>%
    mutate(PCR.positivity_rate = (PCR.n_test_true_positive + PCR.n_test_false_positive)/PCR.n_test,
           RDT.positivity_rate = (RDT.n_test_true_positive + RDT.n_test_false_positive)/RDT.n_test)%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.positivity_rateMean = mean(PCR.positivity_rate),
                     RDT.positivity_rateMean = mean(RDT.positivity_rate),)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.positivity_rate_distrmean = mean(PCR.positivity_rateMean),
                     PCR.positivity_rate_distrmean_SEM = sd(PCR.positivity_rateMean),
                     RDT.positivity_rate_distrmean = mean(RDT.positivity_rateMean),
                     RDT.positivity_rate_distrmean_SEM = sd(RDT.positivity_rateMean))%>%
    mutate(PCR.positivity_rate_distrmean_upper = PCR.positivity_rate_distrmean + 1.96 * PCR.positivity_rate_distrmean_SEM,
           PCR.positivity_rate_distrmean_lower = PCR.positivity_rate_distrmean - 1.96 * PCR.positivity_rate_distrmean_SEM,
           RDT.positivity_rate_distrmean_upper = RDT.positivity_rate_distrmean + 1.96 * RDT.positivity_rate_distrmean_SEM,
           RDT.positivity_rate_distrmean_lower = RDT.positivity_rate_distrmean - 1.96 * RDT.positivity_rate_distrmean_SEM)
  
  # bootstrap
  dataSurveillanceOutputTestPositivity_bootstrap = dataSurveillanceOutput%>%
    mutate(PCR.positivity_rate = (PCR.n_test_true_positive + PCR.n_test_false_positive)/PCR.n_test,
           RDT.positivity_rate = (RDT.n_test_true_positive + RDT.n_test_false_positive)/RDT.n_test)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(PCR.positivity_rate_bootstrap = f_boot(PCR.positivity_rate),
              PCR.positivity_rate_bootstrap_lower = f_boot_lower(PCR.positivity_rate),
              PCR.positivity_rate_bootstrap_upper = f_boot_upper(PCR.positivity_rate),
              RDT.positivity_rate_bootstrap = f_boot(RDT.positivity_rate),
              RDT.positivity_rate_bootstrap_lower = f_boot_lower(RDT.positivity_rate),
              RDT.positivity_rate_bootstrap_upper = f_boot_upper(RDT.positivity_rate))
  
  dataSurveillanceOutputTestPositivity = left_join(dataSurveillanceOutputTestPositivityOverall, 
                                                   left_join(dataSurveillanceOutputTestPositivityTotal,
                                                             left_join(dataSurveillanceOutputTestPositivity_singlemean, 
                                                                       left_join(dataSurveillanceOutputTestPositivity_distrmean, dataSurveillanceOutputTestPositivity_bootstrap))))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputTestPositivity)
}

# # TEST
# f_testPositivity(test1)
# f_testPositivity(test2)



###################
### SENSITIVITY ###
###################

f_sensitivityRealized = function(dataSurveillanceOutput){
  
  # Quantiles
  dataSurveillanceOutputSensitivityOverall = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.sensitivity_realized_025 = quantile(PCR.sensitivity_realized, 0.025, na.rm = T),
                     PCR.sensitivity_realized_250 = quantile(PCR.sensitivity_realized, 0.25, na.rm = T),
                     PCR.sensitivity_realized_500 = quantile(PCR.sensitivity_realized, 0.5, na.rm = T),
                     PCR.sensitivity_realized_750 = quantile(PCR.sensitivity_realized, 0.75, na.rm = T),
                     PCR.sensitivity_realized_975 = quantile(PCR.sensitivity_realized, 0.975, na.rm = T),
                     RDT.sensitivity_realized_025 = quantile(RDT.sensitivity_realized, 0.025, na.rm = T),
                     RDT.sensitivity_realized_250 = quantile(RDT.sensitivity_realized, 0.25, na.rm = T),
                     RDT.sensitivity_realized_500 = quantile(RDT.sensitivity_realized, 0.5, na.rm = T),
                     RDT.sensitivity_realized_750 = quantile(RDT.sensitivity_realized, 0.75, na.rm = T),
                     RDT.sensitivity_realized_975 = quantile(RDT.sensitivity_realized, 0.975, na.rm = T))
  
  
  # overall total as one simulation
  # single means (ONLY INCLUDE OUTBREAKS)
  # distribution of means (ONLY INCLUDE OUTBREAKS)
  
  # Totals
  dataSurveillanceOutputSensitivityTotal = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(PCR.n_test_true_positive_sum = sum(PCR.n_test_true_positive),
              PCR.n_positive_and_tested_sum = sum(PCR.n_positive_and_tested),
              RDT.n_test_true_positive_sum = sum(RDT.n_test_true_positive),
              RDT.n_positive_and_tested_sum = sum(RDT.n_positive_and_tested))%>%
    mutate(PCR.sensitivity_total = PCR.n_test_true_positive_sum/PCR.n_positive_and_tested_sum,
           PCR.sensitivity_total_SD = sqrt((PCR.sensitivity_total*(1-PCR.sensitivity_total))/PCR.n_positive_and_tested_sum),
           PCR.sensitivity_total_lower = PCR.sensitivity_total - 1.96 * PCR.sensitivity_total_SD,
           PCR.sensitivity_total_upper = PCR.sensitivity_total + 1.96 * PCR.sensitivity_total_SD,
           RDT.sensitivity_total = RDT.n_test_true_positive_sum/RDT.n_positive_and_tested_sum,
           RDT.sensitivity_total_SD = sqrt((RDT.sensitivity_total*(1-RDT.sensitivity_total))/RDT.n_positive_and_tested_sum),
           RDT.sensitivity_total_lower = RDT.sensitivity_total - 1.96 * RDT.sensitivity_total_SD,
           RDT.sensitivity_total_upper = RDT.sensitivity_total + 1.96 * RDT.sensitivity_total_SD)
  
  
  # single means
  dataSurveillanceOutputSensitivity_singlemean = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.sensitivity_singlemean = mean(PCR.sensitivity_realized, na.rm=T),
                     PCR.sensitivity_singlemean_SD = sd(PCR.sensitivity_realized, na.rm=T),
                     RDT.sensitivity_singlemean = mean(RDT.sensitivity_realized, na.rm=T),
                     RDT.sensitivity_singlemean_SD = sd(RDT.sensitivity_realized, na.rm=T))%>%
    mutate(PCR.sensitivity_singlemean_upper = PCR.sensitivity_singlemean + 1.96 * PCR.sensitivity_singlemean_SD,
           PCR.sensitivity_singlemean_lower = PCR.sensitivity_singlemean - 1.96 * PCR.sensitivity_singlemean_SD,
           RDT.sensitivity_singlemean_upper = RDT.sensitivity_singlemean + 1.96 * RDT.sensitivity_singlemean_SD,
           RDT.sensitivity_singlemean_lower = RDT.sensitivity_singlemean - 1.96 * RDT.sensitivity_singlemean_SD)
  
  # distribution of means
  dataSurveillanceOutputSensitivity_distrmean = dataSurveillanceOutput%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.sensitivityMean = mean(PCR.sensitivity_realized, na.rm=T),
                     RDT.sensitivityMean = mean(RDT.sensitivity_realized, na.rm=T))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.sensitivity_distrmean = mean(PCR.sensitivityMean, na.rm=T),
                     PCR.sensitivity_distrmean_SEM = sd(PCR.sensitivityMean, na.rm=T),
                     RDT.sensitivity_distrmean = mean(RDT.sensitivityMean, na.rm=T),
                     RDT.sensitivity_distrmean_SEM = sd(RDT.sensitivityMean, na.rm=T))%>%
    mutate(PCR.sensitivity_distrmean_upper = PCR.sensitivity_distrmean + 1.96 * PCR.sensitivity_distrmean_SEM,
           PCR.sensitivity_distrmean_lower = PCR.sensitivity_distrmean - 1.96 * PCR.sensitivity_distrmean_SEM,
           RDT.sensitivity_distrmean_upper = RDT.sensitivity_distrmean + 1.96 * RDT.sensitivity_distrmean_SEM,
           RDT.sensitivity_distrmean_lower = RDT.sensitivity_distrmean - 1.96 * RDT.sensitivity_distrmean_SEM)
  
  # bootstrap
  dataSurveillanceOutputSensitivity_bootstrap = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.sensitivity_bootstrap = f_boot(PCR.sensitivity_realized),
                     PCR.sensitivity_bootstrap_lower = f_boot_lower(PCR.sensitivity_realized),
                     PCR.sensitivity_bootstrap_upper = f_boot_upper(PCR.sensitivity_realized),
                     RDT.sensitivity_bootstrap = f_boot(RDT.sensitivity_realized),
                     RDT.sensitivity_bootstrap_lower = f_boot_lower(RDT.sensitivity_realized),
                     RDT.sensitivity_bootstrap_upper = f_boot_upper(RDT.sensitivity_realized))
  
  # Merge and return
  
  dataSurveillanceOutputSensitivity = left_join(dataSurveillanceOutputSensitivityOverall, 
                                                left_join(dataSurveillanceOutputSensitivityTotal,
                                                          left_join(dataSurveillanceOutputSensitivity_singlemean, 
                                                                    left_join(dataSurveillanceOutputSensitivity_distrmean, dataSurveillanceOutputSensitivity_bootstrap))))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputSensitivity)
}


# # Test
# f_sensitivityRealized(test1)
# f_sensitivityRealized(test2)


###################
### SPECIFICITY ###
###################

f_specificityRealized = function(dataSurveillanceOutput){
  
  # Quantiles
  dataSurveillanceOutputSpecificityOverall = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.specificity_realized_025 = quantile(PCR.specificity_realized, 0.025, na.rm = T),
                     PCR.specificity_realized_250 = quantile(PCR.specificity_realized, 0.25, na.rm = T),
                     PCR.specificity_realized_500 = quantile(PCR.specificity_realized, 0.5, na.rm = T),
                     PCR.specificity_realized_750 = quantile(PCR.specificity_realized, 0.75, na.rm = T),
                     PCR.specificity_realized_975 = quantile(PCR.specificity_realized, 0.975, na.rm = T),
                     RDT.specificity_realized_025 = quantile(RDT.specificity_realized, 0.025, na.rm = T),
                     RDT.specificity_realized_250 = quantile(RDT.specificity_realized, 0.25, na.rm = T),
                     RDT.specificity_realized_500 = quantile(RDT.specificity_realized, 0.5, na.rm = T),
                     RDT.specificity_realized_750 = quantile(RDT.specificity_realized, 0.75, na.rm = T),
                     RDT.specificity_realized_975 = quantile(RDT.specificity_realized, 0.975, na.rm = T))
  
  
  # overall total as one simulation
  # single means (ONLY INCLUDE OUTBREAKS)
  # distribution of means (ONLY INCLUDE OUTBREAKS)
  
  # Totals
  dataSurveillanceOutputSpecificityTotal = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(PCR.n_test_true_negative_sum = sum(PCR.n_test_true_negative),
              PCR.n_negative_and_tested_sum = sum(PCR.n_negative_and_tested),
              RDT.n_test_true_negative_sum = sum(RDT.n_test_true_negative),
              RDT.n_negative_and_tested_sum = sum(RDT.n_negative_and_tested))%>%
    mutate(PCR.specificity_total = PCR.n_test_true_negative_sum/PCR.n_negative_and_tested_sum,
           PCR.specificity_total_SD = sqrt((PCR.specificity_total*(1-PCR.specificity_total))/PCR.n_negative_and_tested_sum),
           PCR.specificity_total_lower = PCR.specificity_total - 1.96 * PCR.specificity_total_SD,
           PCR.specificity_total_upper = PCR.specificity_total + 1.96 * PCR.specificity_total_SD,
           RDT.specificity_total = RDT.n_test_true_negative_sum/RDT.n_negative_and_tested_sum,
           RDT.specificity_total_SD = sqrt((RDT.specificity_total*(1-RDT.specificity_total))/RDT.n_negative_and_tested_sum),
           RDT.specificity_total_lower = RDT.specificity_total - 1.96 * RDT.specificity_total_SD,
           RDT.specificity_total_upper = RDT.specificity_total + 1.96 * RDT.specificity_total_SD)
  
  # single means
  dataSurveillanceOutputSpecificity_singlemean = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.specificity_singlemean = mean(PCR.specificity_realized, na.rm=T),
                     PCR.specificity_singlemean_SD = sd(PCR.specificity_realized, na.rm=T),
                     RDT.specificity_singlemean = mean(RDT.specificity_realized, na.rm=T),
                     RDT.specificity_singlemean_SD = sd(RDT.specificity_realized, na.rm=T))%>%
    mutate(PCR.specificity_singlemean_upper = PCR.specificity_singlemean + 1.96 * PCR.specificity_singlemean_SD,
           PCR.specificity_singlemean_lower = PCR.specificity_singlemean - 1.96 * PCR.specificity_singlemean_SD,
           RDT.specificity_singlemean_upper = RDT.specificity_singlemean + 1.96 * RDT.specificity_singlemean_SD,
           RDT.specificity_singlemean_lower = RDT.specificity_singlemean - 1.96 * RDT.specificity_singlemean_SD)
  
  # distribution of means
  dataSurveillanceOutputSpecificity_distrmean = dataSurveillanceOutput%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.specificityMean = mean(PCR.specificity_realized, na.rm=T),
                     RDT.specificityMean = mean(RDT.specificity_realized, na.rm=T))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.specificity_distrmean = mean(PCR.specificityMean, na.rm=T),
                     PCR.specificity_distrmean_SEM = sd(PCR.specificityMean, na.rm=T),
                     RDT.specificity_distrmean = mean(RDT.specificityMean, na.rm=T),
                     RDT.specificity_distrmean_SEM = sd(RDT.specificityMean, na.rm=T))%>%
    mutate(PCR.specificity_distrmean_upper = PCR.specificity_distrmean + 1.96 * PCR.specificity_distrmean_SEM,
           PCR.specificity_distrmean_lower = PCR.specificity_distrmean - 1.96 * PCR.specificity_distrmean_SEM,
           RDT.specificity_distrmean_upper = RDT.specificity_distrmean + 1.96 * RDT.specificity_distrmean_SEM,
           RDT.specificity_distrmean_lower = RDT.specificity_distrmean - 1.96 * RDT.specificity_distrmean_SEM)
  
  # bootstrap
  dataSurveillanceOutputSpecificity_bootstrap = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.specificity_bootstrap = f_boot(PCR.specificity_realized),
                     PCR.specificity_bootstrap_lower = f_boot_lower(PCR.specificity_realized),
                     PCR.specificity_bootstrap_upper = f_boot_upper(PCR.specificity_realized),
                     RDT.specificity_bootstrap = f_boot(RDT.specificity_realized),
                     RDT.specificity_bootstrap_lower = f_boot_lower(RDT.specificity_realized),
                     RDT.specificity_bootstrap_upper = f_boot_upper(RDT.specificity_realized))
  
  # Merge and return
  
  dataSurveillanceOutputSpecificity = left_join(dataSurveillanceOutputSpecificityOverall,
                                                left_join(dataSurveillanceOutputSpecificityTotal, 
                                                          left_join(dataSurveillanceOutputSpecificity_singlemean, 
                                                                    left_join(dataSurveillanceOutputSpecificity_distrmean, dataSurveillanceOutputSpecificity_bootstrap))))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputSpecificity)
}


# # Test
# f_specificityRealized(test1)
# f_specificityRealized(test2)


###########
### PPV ### Positive predictive value
###########

f_PPV = function(dataSurveillanceOutput){
  
  # Quantiles
  dataSurveillanceOutputPPVOverall = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.PPV_025 = quantile(PCR.PPV, 0.025, na.rm = T),
                     PCR.PPV_250 = quantile(PCR.PPV, 0.25, na.rm = T),
                     PCR.PPV_500 = quantile(PCR.PPV, 0.5, na.rm = T),
                     PCR.PPV_750 = quantile(PCR.PPV, 0.75, na.rm = T),
                     PCR.PPV_975 = quantile(PCR.PPV, 0.975, na.rm = T),
                     RDT.PPV_025 = quantile(RDT.PPV, 0.025, na.rm = T),
                     RDT.PPV_250 = quantile(RDT.PPV, 0.25, na.rm = T),
                     RDT.PPV_500 = quantile(RDT.PPV, 0.5, na.rm = T),
                     RDT.PPV_750 = quantile(RDT.PPV, 0.75, na.rm = T),
                     RDT.PPV_975 = quantile(RDT.PPV, 0.975, na.rm = T))
  
  
  # overall total as one simulation
  # single means (ONLY INCLUDE OUTBREAKS)
  # distribution of means (ONLY INCLUDE OUTBREAKS)
  
  # Totals
  dataSurveillanceOutputPPVTotal = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(PCR.n_test_true_positive_sum = sum(PCR.n_test_true_positive),
              PCR.n_test_positive_sum = sum(PCR.n_test_true_positive) + sum(PCR.n_test_false_positive),
              RDT.n_test_true_positive_sum = sum(RDT.n_test_true_positive),
              RDT.n_test_positive_sum = sum(RDT.n_test_true_positive) + sum(RDT.n_test_false_positive))%>%
    mutate(PCR.PPV_total = PCR.n_test_true_positive_sum/PCR.n_test_positive_sum,
           PCR.PPV_total_SD = sqrt((PCR.PPV_total*(1-PCR.PPV_total))/PCR.n_test_positive_sum),
           PCR.PPV_total_lower = PCR.PPV_total - 1.96 * PCR.PPV_total_SD,
           PCR.PPV_total_upper = PCR.PPV_total + 1.96 * PCR.PPV_total_SD,
           RDT.PPV_total = RDT.n_test_true_positive_sum/RDT.n_test_positive_sum,
           RDT.PPV_total_SD = sqrt((RDT.PPV_total*(1-RDT.PPV_total))/RDT.n_test_positive_sum),
           RDT.PPV_total_lower = RDT.PPV_total - 1.96 * RDT.PPV_total_SD,
           RDT.PPV_total_upper = RDT.PPV_total + 1.96 * RDT.PPV_total_SD)
  
  # single means
  dataSurveillanceOutputPPV_singlemean = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.PPV_singlemean = mean(PCR.PPV, na.rm=T),
                     PCR.PPV_singlemean_SD = sd(PCR.PPV, na.rm=T),
                     RDT.PPV_singlemean = mean(RDT.PPV, na.rm=T),
                     RDT.PPV_singlemean_SD = sd(RDT.PPV, na.rm=T))%>%
    mutate(PCR.PPV_singlemean_upper = PCR.PPV_singlemean + 1.96 * PCR.PPV_singlemean_SD,
           PCR.PPV_singlemean_lower = PCR.PPV_singlemean - 1.96 * PCR.PPV_singlemean_SD,
           RDT.PPV_singlemean_upper = RDT.PPV_singlemean + 1.96 * RDT.PPV_singlemean_SD,
           RDT.PPV_singlemean_lower = RDT.PPV_singlemean - 1.96 * RDT.PPV_singlemean_SD)
  
  # distribution of means
  dataSurveillanceOutputPPV_distrmean = dataSurveillanceOutput%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.PPVMean = mean(PCR.PPV, na.rm=T),
                     RDT.PPVMean = mean(RDT.PPV, na.rm=T))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.PPV_distrmean = mean(PCR.PPVMean, na.rm=T),
                     PCR.PPV_distrmean_SEM = sd(PCR.PPVMean, na.rm=T),
                     RDT.PPV_distrmean = mean(RDT.PPVMean, na.rm=T),
                     RDT.PPV_distrmean_SEM = sd(RDT.PPVMean, na.rm=T))%>%
    mutate(PCR.PPV_distrmean_upper = PCR.PPV_distrmean + 1.96 * PCR.PPV_distrmean_SEM,
           PCR.PPV_distrmean_lower = PCR.PPV_distrmean - 1.96 * PCR.PPV_distrmean_SEM,
           RDT.PPV_distrmean_upper = RDT.PPV_distrmean + 1.96 * RDT.PPV_distrmean_SEM,
           RDT.PPV_distrmean_lower = RDT.PPV_distrmean - 1.96 * RDT.PPV_distrmean_SEM)
  
  # bootstrap (need to filter out 0 denominators (simulations with no positive tests))
  dataSurveillanceOutputPPV_bootstrap.PCR = dataSurveillanceOutput%>%
    filter(PCR.n_test_true_positive + PCR.n_test_false_positive > 0)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.PPV_bootstrap = f_boot(PCR.PPV),
                     PCR.PPV_bootstrap_lower = f_boot_lower(PCR.PPV),
                     PCR.PPV_bootstrap_upper = f_boot_upper(PCR.PPV))
  
  dataSurveillanceOutputPPV_bootstrap.RDT = dataSurveillanceOutput%>%
    filter(RDT.n_test_true_positive + RDT.n_test_false_positive > 0)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(RDT.PPV_bootstrap = f_boot(RDT.PPV),
                     RDT.PPV_bootstrap_lower = f_boot_lower(RDT.PPV),
                     RDT.PPV_bootstrap_upper = f_boot_upper(RDT.PPV))
  
  dataSurveillanceOutputPPV_bootstrap = left_join(dataSurveillanceOutputPPV_bootstrap.PCR, dataSurveillanceOutputPPV_bootstrap.RDT)
  
  
  # Merge and return
  
  dataSurveillanceOutputPPV = left_join(dataSurveillanceOutputPPVTotal, 
                                        left_join(dataSurveillanceOutputPPVOverall,
                                                  left_join(dataSurveillanceOutputPPV_singlemean, 
                                                            left_join(dataSurveillanceOutputPPV_distrmean, dataSurveillanceOutputPPV_bootstrap))))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputPPV)
}

# # Test
# f_PPV(test1)
# f_PPV(test2)

###########
### NPV ### Negative predictive value
###########

f_NPV = function(dataSurveillanceOutput){
  
  # Quantiles
  dataSurveillanceOutputNPVOverall = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.NPV_025 = quantile(PCR.NPV, 0.025, na.rm = T),
                     PCR.NPV_250 = quantile(PCR.NPV, 0.25, na.rm = T),
                     PCR.NPV_500 = quantile(PCR.NPV, 0.5, na.rm = T),
                     PCR.NPV_750 = quantile(PCR.NPV, 0.75, na.rm = T),
                     PCR.NPV_975 = quantile(PCR.NPV, 0.975, na.rm = T),
                     RDT.NPV_025 = quantile(RDT.NPV, 0.025, na.rm = T),
                     RDT.NPV_250 = quantile(RDT.NPV, 0.25, na.rm = T),
                     RDT.NPV_500 = quantile(RDT.NPV, 0.5, na.rm = T),
                     RDT.NPV_750 = quantile(RDT.NPV, 0.75, na.rm = T),
                     RDT.NPV_975 = quantile(RDT.NPV, 0.975, na.rm = T))
  
  
  # overall total as one simulation
  # single means (ONLY INCLUDE OUTBREAKS)
  # distribution of means (ONLY INCLUDE OUTBREAKS)
  
  # Totals
  dataSurveillanceOutputNPVTotal = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(PCR.n_test_true_negative_sum = sum(PCR.n_test_true_negative),
              PCR.n_test_negative_sum = sum(PCR.n_test_true_negative) + sum(PCR.n_test_false_negative),
              RDT.n_test_true_negative_sum = sum(RDT.n_test_true_negative),
              RDT.n_test_negative_sum = sum(RDT.n_test_true_negative) + sum(RDT.n_test_false_negative))%>%
    mutate(PCR.NPV_total = PCR.n_test_true_negative_sum/PCR.n_test_negative_sum,
           PCR.NPV_total_SD = sqrt((PCR.NPV_total*(1-PCR.NPV_total))/PCR.n_test_negative_sum),
           PCR.NPV_total_lower = PCR.NPV_total - 1.96 * PCR.NPV_total_SD,
           PCR.NPV_total_upper = PCR.NPV_total + 1.96 * PCR.NPV_total_SD,
           RDT.NPV_total = RDT.n_test_true_negative_sum/RDT.n_test_negative_sum,
           RDT.NPV_total_SD = sqrt((RDT.NPV_total*(1-RDT.NPV_total))/RDT.n_test_negative_sum),
           RDT.NPV_total_lower = RDT.NPV_total - 1.96 * RDT.NPV_total_SD,
           RDT.NPV_total_upper = RDT.NPV_total + 1.96 * RDT.NPV_total_SD)
  
  # single means
  dataSurveillanceOutputNPV_singlemean = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.NPV_singlemean = mean(PCR.NPV, na.rm=T),
                     PCR.NPV_singlemean_SD = sd(PCR.NPV, na.rm=T),
                     RDT.NPV_singlemean = mean(RDT.NPV, na.rm=T),
                     RDT.NPV_singlemean_SD = sd(RDT.NPV, na.rm=T))%>%
    mutate(PCR.NPV_singlemean_upper = PCR.NPV_singlemean + 1.96 * PCR.NPV_singlemean_SD,
           PCR.NPV_singlemean_lower = PCR.NPV_singlemean - 1.96 * PCR.NPV_singlemean_SD,
           RDT.NPV_singlemean_upper = RDT.NPV_singlemean + 1.96 * RDT.NPV_singlemean_SD,
           RDT.NPV_singlemean_lower = RDT.NPV_singlemean - 1.96 * RDT.NPV_singlemean_SD)
  
  # distribution of means
  dataSurveillanceOutputNPV_distrmean = dataSurveillanceOutput%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.NPVMean = mean(PCR.NPV, na.rm=T),
                     RDT.NPVMean = mean(RDT.NPV, na.rm=T))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.NPV_distrmean = mean(PCR.NPVMean, na.rm=T),
                     PCR.NPV_distrmean_SEM = sd(PCR.NPVMean, na.rm=T),
                     RDT.NPV_distrmean = mean(RDT.NPVMean, na.rm=T),
                     RDT.NPV_distrmean_SEM = sd(RDT.NPVMean, na.rm=T))%>%
    mutate(PCR.NPV_distrmean_upper = PCR.NPV_distrmean + 1.96 * PCR.NPV_distrmean_SEM,
           PCR.NPV_distrmean_lower = PCR.NPV_distrmean - 1.96 * PCR.NPV_distrmean_SEM,
           RDT.NPV_distrmean_upper = RDT.NPV_distrmean + 1.96 * RDT.NPV_distrmean_SEM,
           RDT.NPV_distrmean_lower = RDT.NPV_distrmean - 1.96 * RDT.NPV_distrmean_SEM)
  
  # bootstrap
  dataSurveillanceOutputNPV_bootstrap.PCR = dataSurveillanceOutput%>%
    filter(PCR.n_test_true_negative + PCR.n_test_false_negative > 0)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.NPV_bootstrap = f_boot(PCR.NPV),
                     PCR.NPV_bootstrap_lower = f_boot_lower(PCR.NPV),
                     PCR.NPV_bootstrap_upper = f_boot_upper(PCR.NPV))
  
  dataSurveillanceOutputNPV_bootstrap.RDT = dataSurveillanceOutput%>%
    filter(RDT.n_test_true_negative + RDT.n_test_false_negative > 0)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(RDT.NPV_bootstrap = f_boot(RDT.NPV),
                     RDT.NPV_bootstrap_lower = f_boot_lower(RDT.NPV),
                     RDT.NPV_bootstrap_upper = f_boot_upper(RDT.NPV))
  
  dataSurveillanceOutputNPV_bootstrap = left_join(dataSurveillanceOutputNPV_bootstrap.PCR, dataSurveillanceOutputNPV_bootstrap.RDT)
  
  
  # Merge and return
  dataSurveillanceOutputNPV = left_join(dataSurveillanceOutputNPVTotal, 
                                        left_join(dataSurveillanceOutputNPVOverall,
                                                  left_join(dataSurveillanceOutputNPV_singlemean, 
                                                            left_join(dataSurveillanceOutputNPV_distrmean, dataSurveillanceOutputNPV_bootstrap))))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputNPV)
}


# # Test
# f_NPV(test1)
# f_NPV(test2)



###############################
### CASES DETECTED PER TEST ###
###############################

f_casesDetectedPerTest = function(dataSurveillanceOutput){
  
  dataSurveillanceOutputCasesDetectedPerTestOverall = dataSurveillanceOutput%>%
    mutate(PCR.detections_per_test = PCR.n_test_true_positive/PCR.n_test,
           RDT.detections_per_test = RDT.n_test_true_positive/RDT.n_test)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.detections_per_test_025 = quantile(PCR.detections_per_test, 0.025, na.rm = T),
                     PCR.detections_per_test_250 = quantile(PCR.detections_per_test, 0.25, na.rm = T),
                     PCR.detections_per_test_500 = quantile(PCR.detections_per_test, 0.5, na.rm = T),
                     PCR.detections_per_test_750 = quantile(PCR.detections_per_test, 0.75, na.rm = T),
                     PCR.detections_per_test_975 = quantile(PCR.detections_per_test, 0.975, na.rm = T),
                     RDT.detections_per_test_025 = quantile(RDT.detections_per_test, 0.025, na.rm = T),
                     RDT.detections_per_test_250 = quantile(RDT.detections_per_test, 0.25, na.rm = T),
                     RDT.detections_per_test_500 = quantile(RDT.detections_per_test, 0.5, na.rm = T),
                     RDT.detections_per_test_750 = quantile(RDT.detections_per_test, 0.75, na.rm = T),
                     RDT.detections_per_test_975 = quantile(RDT.detections_per_test, 0.975, na.rm = T))
  
  
  # overall total as one simulation
  # single means (ONLY INCLUDE OUTBREAKS)
  # distribution of means (ONLY INCLUDE OUTBREAKS)
  
  dataSurveillanceOutputCasesDetectedPerTestTotal = dataSurveillanceOutput%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(PCR.n_test_sum = sum(PCR.n_test),
              PCR.n_test_true_positive_sum = sum(PCR.n_test_true_positive),
              RDT.n_test_sum = sum(RDT.n_test),
              RDT.n_test_true_positive_sum = sum(RDT.n_test_true_positive))%>%
    mutate(PCR.detections_per_test_total = PCR.n_test_true_positive_sum/PCR.n_test_sum,
           PCR.detections_per_test_total_SD = sqrt(PCR.detections_per_test_total*(1-PCR.detections_per_test_total)/PCR.n_test_sum),
           PCR.detections_per_test_total_lower = PCR.detections_per_test_total - 1.96 * PCR.detections_per_test_total_SD,
           PCR.detections_per_test_total_upper = PCR.detections_per_test_total + 1.96 * PCR.detections_per_test_total_SD,
           RDT.detections_per_test_total = RDT.n_test_true_positive_sum/RDT.n_test_sum,
           RDT.detections_per_test_total_SD = sqrt(RDT.detections_per_test_total*(1-RDT.detections_per_test_total)/RDT.n_test_sum),
           RDT.detections_per_test_total_lower = RDT.detections_per_test_total - 1.96 * RDT.detections_per_test_total_SD,
           RDT.detections_per_test_total_upper = RDT.detections_per_test_total + 1.96 * RDT.detections_per_test_total_SD)
  
  # single means
  dataSurveillanceOutputCasesDetectedPerTest_singlemean = dataSurveillanceOutput%>%
    mutate(PCR.detections_per_test = PCR.n_test_true_positive/PCR.n_test,
           RDT.detections_per_test = RDT.n_test_true_positive/RDT.n_test)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.detections_per_test_singlemean = mean(PCR.detections_per_test, na.rm=T),
                     PCR.detections_per_test_singlemean_SD = sd(PCR.detections_per_test, na.rm=T),
                     RDT.detections_per_test_singlemean = mean(RDT.detections_per_test, na.rm=T),
                     RDT.detections_per_test_singlemean_SD = sd(RDT.detections_per_test, na.rm=T))%>%
    mutate(PCR.detections_per_test_singlemean_upper = PCR.detections_per_test_singlemean + 1.96 * PCR.detections_per_test_singlemean_SD,
           PCR.detections_per_test_singlemean_lower = PCR.detections_per_test_singlemean - 1.96 * PCR.detections_per_test_singlemean_SD,
           RDT.detections_per_test_singlemean_upper = RDT.detections_per_test_singlemean + 1.96 * RDT.detections_per_test_singlemean_SD,
           RDT.detections_per_test_singlemean_lower = RDT.detections_per_test_singlemean - 1.96 * RDT.detections_per_test_singlemean_SD)
  
  # distribution of means
  dataSurveillanceOutputCasesDetectedPerTest_distrmean = dataSurveillanceOutput%>%
    mutate(PCR.detections_per_test = PCR.n_test_true_positive/PCR.n_test,
           RDT.detections_per_test = RDT.n_test_true_positive/RDT.n_test)%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.detections_per_testMean = mean(PCR.detections_per_test, na.rm=T),
                     RDT.detections_per_testMean = mean(RDT.detections_per_test, na.rm=T))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.detections_per_test_distrmean = mean(PCR.detections_per_testMean, na.rm=T),
                     PCR.detections_per_test_distrmean_SEM = sd(PCR.detections_per_testMean, na.rm=T),
                     RDT.detections_per_test_distrmean = mean(RDT.detections_per_testMean, na.rm=T),
                     RDT.detections_per_test_distrmean_SEM = sd(RDT.detections_per_testMean, na.rm=T))%>%
    mutate(PCR.detections_per_test_distrmean_upper = PCR.detections_per_test_distrmean + 1.96 * PCR.detections_per_test_distrmean_SEM,
           PCR.detections_per_test_distrmean_lower = PCR.detections_per_test_distrmean - 1.96 * PCR.detections_per_test_distrmean_SEM,
           RDT.detections_per_test_distrmean_upper = RDT.detections_per_test_distrmean + 1.96 * RDT.detections_per_test_distrmean_SEM,
           RDT.detections_per_test_distrmean_lower = RDT.detections_per_test_distrmean - 1.96 * RDT.detections_per_test_distrmean_SEM)
  
  
  # bootstrap
  dataSurveillanceOutputCasesDetectedPerTest_bootstrap = dataSurveillanceOutput%>%
    mutate(PCR.detections_per_test = PCR.n_test_true_positive/PCR.n_test,
           RDT.detections_per_test = RDT.n_test_true_positive/RDT.n_test)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.detections_per_test_bootstrap = f_boot(PCR.detections_per_test),
                     PCR.detections_per_test_bootstrap_lower = f_boot_lower(PCR.detections_per_test),
                     PCR.detections_per_test_bootstrap_upper = f_boot_upper(PCR.detections_per_test),
                     RDT.detections_per_test_bootstrap = f_boot(RDT.detections_per_test),
                     RDT.detections_per_test_bootstrap_lower = f_boot_lower(RDT.detections_per_test),
                     RDT.detections_per_test_bootstrap_upper = f_boot_upper(RDT.detections_per_test))
  
  # Merge and return
  dataSurveillanceOutputCasesDetectedPerTest = left_join(dataSurveillanceOutputCasesDetectedPerTestOverall, 
                                                         left_join(dataSurveillanceOutputCasesDetectedPerTestTotal,
                                                                   left_join(dataSurveillanceOutputCasesDetectedPerTest_singlemean,
                                                                             left_join(dataSurveillanceOutputCasesDetectedPerTest_distrmean,
                                                                                       dataSurveillanceOutputCasesDetectedPerTest_bootstrap))))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputCasesDetectedPerTest)
  
}

# # TEST
# f_casesDetectedPerTest(test1)
# f_casesDetectedPerTest(test2)



##############################
### CASES AVERTED PER TEST ###
##############################

f_casesAvertedPerTest = function(dataSurveillanceOutput){
  
  ### only PCR: all {screeningType == 'PCR'} OR {stratSurveillance %in% c('adm', 'sym', 'adm_sym')}
  dataSurveillanceOutput_PCRonly = dplyr::filter(dataSurveillanceOutput, typeDailyScreening == 'PCR' | stratSurveillance %in% c('adm', 'sym', 'adm_sym'))%>%
    mutate(PCR.infAverted = infAverted)%>%
    mutate(RDT.infAverted = 0)
  
  
  ### only RDT: all {screeningType %in% c('RDT1', 'RDT2') AND stratSurveillance %in% c(d01, d02, d03, d04, d05, d06, d07)} 
  dataSurveillanceOutput_RDTonly = dplyr::filter(dataSurveillanceOutput, typeDailyScreening %in% c('RDT1', 'RDT2') & stratSurveillance %in% paste0('d0',1:9))%>%
    mutate(PCR.infAverted = 0)%>%
    mutate(RDT.infAverted = infAverted)
  
  ### mixed PCR/RDT: all {screeningType %in% c('RDT1', 'RDT2') AND stratSurveillance %in% c(adm_sym_dX) for all X}
  
  defaultTargetDailyScreening = first(dataSurveillanceOutput$targetDailyScreening)
  defaultTypeDailyScreening = first(dataSurveillanceOutput$typeDailyScreening)
  
  # Isolate PCR to calculate infections averted from PCR testing
  dataSurveillanceOutput_mixedPCR.RDT_justPCR = dplyr::filter(dataSurveillanceOutput,
                                                              targetDailyScreening == defaultTargetDailyScreening,
                                                              typeDailyScreening == defaultTypeDailyScreening,
                                                              stratSurveillance == 'adm_sym')%>%
    mutate(PCR.infAverted = infAverted)%>%
    dplyr::select(simulationCTC, testingRound, PCR.n_test, PCR.infAverted, sensitivity_function)
  
  # Remove PCR-detected infections to isolate infections averted from RDT testing
  dataSurveillanceOutput_mixedPCR.RDT_noPCR = dplyr::filter(dataSurveillanceOutput, 
                                                            typeDailyScreening %in% c('RDT1', 'RDT2') & !stratSurveillance %in% c('adm', 'sym', 'adm_sym', paste0('d0',1:9)))%>%
    select(-c(PCR.n_test))
  
  # Calcualte infections averted from RDT testing
  dataSurveillanceOutput_mixedPCR.RDT = left_join(dataSurveillanceOutput_mixedPCR.RDT_justPCR, 
                                                  dataSurveillanceOutput_mixedPCR.RDT_noPCR, 
                                                  by = c('simulationCTC', 'testingRound', 'sensitivity_function'))%>%
    mutate(RDT.infAverted = infAverted - PCR.infAverted)
  
  
  ### PUT ALL TOGETHER
  dataSurveillanceOutputTestsPerInfAverted = rbind(dataSurveillanceOutput_PCRonly, dataSurveillanceOutput_RDTonly, dataSurveillanceOutput_mixedPCR.RDT)%>%
    mutate(PCR.averted_per_test = PCR.infAverted/PCR.n_test,
           RDT.averted_per_test = RDT.infAverted/RDT.n_test)
  
  dataSurveillanceOutputCasesAvertedPerTestOverall = dataSurveillanceOutputTestsPerInfAverted%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.averted_per_test_025 = quantile(PCR.averted_per_test, 0.025, na.rm = T),
                     PCR.averted_per_test_250 = quantile(PCR.averted_per_test, 0.25, na.rm = T),
                     PCR.averted_per_test_500 = quantile(PCR.averted_per_test, 0.5, na.rm = T),
                     PCR.averted_per_test_750 = quantile(PCR.averted_per_test, 0.75, na.rm = T),
                     PCR.averted_per_test_975 = quantile(PCR.averted_per_test, 0.975, na.rm = T),
                     RDT.averted_per_test_025 = quantile(RDT.averted_per_test, 0.025, na.rm = T),
                     RDT.averted_per_test_250 = quantile(RDT.averted_per_test, 0.25, na.rm = T),
                     RDT.averted_per_test_500 = quantile(RDT.averted_per_test, 0.5, na.rm = T),
                     RDT.averted_per_test_750 = quantile(RDT.averted_per_test, 0.75, na.rm = T),
                     RDT.averted_per_test_975 = quantile(RDT.averted_per_test, 0.975, na.rm = T))
  
  
  # overall total as one simulation
  # single means (ONLY INCLUDE OUTBREAKS)
  # distribution of means (ONLY INCLUDE OUTBREAKS)
  
  dataSurveillanceOutputCasesAvertedPerTestTotal = dataSurveillanceOutputTestsPerInfAverted%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    summarise(PCR.n_test_sum = sum(PCR.n_test),
              PCR.n_infAverted_sum = sum(PCR.infAverted),
              RDT.n_test_sum = sum(RDT.n_test),
              RDT.n_infAverted_sum = sum(RDT.infAverted))%>%
    mutate(PCR.averted_per_test_total = PCR.n_infAverted_sum/PCR.n_test_sum,
           PCR.averted_per_test_total_SD = sqrt(PCR.averted_per_test_total*(1-PCR.averted_per_test_total)/PCR.n_test_sum),
           PCR.averted_per_test_total_lower = PCR.averted_per_test_total - 1.96 * PCR.averted_per_test_total_SD,
           PCR.averted_per_test_total_upper = PCR.averted_per_test_total + 1.96 * PCR.averted_per_test_total_SD,
           RDT.averted_per_test_total = RDT.n_infAverted_sum/RDT.n_test_sum,
           RDT.averted_per_test_total_SD = sqrt(RDT.averted_per_test_total*(1-RDT.averted_per_test_total)/RDT.n_test_sum),
           RDT.averted_per_test_total_lower = RDT.averted_per_test_total - 1.96 * RDT.averted_per_test_total_SD,
           RDT.averted_per_test_total_upper = RDT.averted_per_test_total + 1.96 * RDT.averted_per_test_total_SD)
  
  # single means
  dataSurveillanceOutputCasesAvertedPerTest_singlemean = dataSurveillanceOutputTestsPerInfAverted%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.averted_per_test_singlemean = mean(PCR.averted_per_test, na.rm=T),
                     PCR.averted_per_test_singlemean_SD = sd(PCR.averted_per_test, na.rm=T),
                     RDT.averted_per_test_singlemean = mean(RDT.averted_per_test, na.rm=T),
                     RDT.averted_per_test_singlemean_SD = sd(RDT.averted_per_test, na.rm=T))%>%
    mutate(PCR.averted_per_test_singlemean_upper = PCR.averted_per_test_singlemean + 1.96 * PCR.averted_per_test_singlemean_SD,
           PCR.averted_per_test_singlemean_lower = PCR.averted_per_test_singlemean - 1.96 * PCR.averted_per_test_singlemean_SD,
           RDT.averted_per_test_singlemean_upper = RDT.averted_per_test_singlemean + 1.96 * RDT.averted_per_test_singlemean_SD,
           RDT.averted_per_test_singlemean_lower = RDT.averted_per_test_singlemean - 1.96 * RDT.averted_per_test_singlemean_SD)
  
  # distribution of means
  dataSurveillanceOutputCasesAvertedPerTest_distrmean = dataSurveillanceOutputTestsPerInfAverted%>%
    group_by(simulationCTC, stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.averted_per_testMean = mean(PCR.averted_per_test, na.rm=T),
                     RDT.averted_per_testMean = mean(RDT.averted_per_test, na.rm=T))%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.averted_per_test_distrmean = mean(PCR.averted_per_testMean, na.rm=T),
                     PCR.averted_per_test_distrmean_SEM = sd(PCR.averted_per_testMean, na.rm=T),
                     RDT.averted_per_test_distrmean = mean(RDT.averted_per_testMean, na.rm=T),
                     RDT.averted_per_test_distrmean_SEM = sd(RDT.averted_per_testMean, na.rm=T))%>%
    mutate(PCR.averted_per_test_distrmean_upper = PCR.averted_per_test_distrmean + 1.96 * PCR.averted_per_test_distrmean_SEM,
           PCR.averted_per_test_distrmean_lower = PCR.averted_per_test_distrmean - 1.96 * PCR.averted_per_test_distrmean_SEM,
           RDT.averted_per_test_distrmean_upper = RDT.averted_per_test_distrmean + 1.96 * RDT.averted_per_test_distrmean_SEM,
           RDT.averted_per_test_distrmean_lower = RDT.averted_per_test_distrmean - 1.96 * RDT.averted_per_test_distrmean_SEM)
  
  # bootstrap
  dataSurveillanceOutputCasesAvertedPerTest_bootstrap = dataSurveillanceOutputTestsPerInfAverted%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function)%>%
    dplyr::summarise(PCR.averted_per_test_bootstrap = f_boot(PCR.averted_per_test),
                     PCR.averted_per_test_bootstrap_lower = f_boot_lower(PCR.averted_per_test),
                     PCR.averted_per_test_bootstrap_upper = f_boot_upper(PCR.averted_per_test),
                     RDT.averted_per_test_bootstrap = f_boot(RDT.averted_per_test),
                     RDT.averted_per_test_bootstrap_lower = f_boot_lower(RDT.averted_per_test),
                     RDT.averted_per_test_bootstrap_upper = f_boot_upper(RDT.averted_per_test))
  
  # Merge and return
  dataSurveillanceOutputCasesAvertedPerTest = left_join(dataSurveillanceOutputCasesAvertedPerTestOverall, 
                                                        left_join(dataSurveillanceOutputCasesAvertedPerTestTotal,
                                                                  left_join(dataSurveillanceOutputCasesAvertedPerTest_singlemean, 
                                                                            left_join(dataSurveillanceOutputCasesAvertedPerTest_distrmean,
                                                                                      dataSurveillanceOutputCasesAvertedPerTest_bootstrap))))%>%
    as.data.frame()
  
  return(dataSurveillanceOutputCasesAvertedPerTest)
  
}

# # TEST
# f_casesAvertedPerTest(test1)
# f_casesAvertedPerTest(test2)


#############################
### COST PER CASE AVERTED ###
#############################

f_costPerCaseAverted = function(dataSurveillanceOutput, cost_rdt, cost_pcr){
  
  dataCost = dataSurveillanceOutput%>%
    dplyr::select(simulationCTC, testingRound, sensitivity_function, stratSurveillance, 
                  typeDailyScreening, targetDailyScreening, infAverted, PCR.n_test, RDT.n_test)%>%
    mutate(cost_pcr = cost_pcr, cost_rdt = cost_rdt, 
           cost = PCR.n_test * cost_pcr + RDT.n_test * cost_rdt,
           caseAverted_per_cost = infAverted/cost)%>%
    group_by(stratSurveillance, typeDailyScreening, targetDailyScreening, sensitivity_function, cost_pcr, cost_rdt)%>%
    dplyr::summarise(caseAverted_per_cost_bootstrap = f_boot(caseAverted_per_cost),
                     caseAverted_per_cost_bootstrap_lower = f_boot_lower(caseAverted_per_cost),
                     caseAverted_per_cost_bootstrap_upper = f_boot_upper(caseAverted_per_cost))%>%
    mutate(cost_per_caseAverted = 1/caseAverted_per_cost_bootstrap,
           cost_per_caseAverted_lower = 1/caseAverted_per_cost_bootstrap_upper,
           cost_per_caseAverted_upper = 1/caseAverted_per_cost_bootstrap_lower)
  
  return(dataCost)
}

# # TEST
# f_costPerCaseAverted(test1, cost_rdt = 5, cost_pcr = 50)
# f_costPerCaseAverted(test2, cost_rdt = 5, cost_pcr = 50)