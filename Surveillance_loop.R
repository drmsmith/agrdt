source("functions.R")

###################################################################################
### LOOPING THROUGH SIMULATIONS AND TESTING STRATEGIES TO SIMULATE SURVEILLANCE ###
###################################################################################

### Output Lot 3: baseline "low incidence" scenario
lot = "lot3"
lot_x = "lot3"

v_LTCF_x = 1#:3 # 3 baseline LTCFs considered
v_sens_function_x = c('time-varying','uniform', 'perfect') # 3 sensitivity functions considered
files_CTC = 1#:100 # range of CTC outbreaks to evaluate (from min 1 to max 100)
testing_rounds_x = 100 # number of surveillance rounds per surveillance strategy per outbreak (100 in manuscript)


vec_screening_types = c("PCR", "RDT1", "RDT2") # types of test considered for screening
vec_screening_targets = c("Patients", "Staff", "All") # types of individuals considered for screening
vec_surveillance_strategies = c( 
  "adm","sym","adm_sym",
  "d01","d02","d03","d04","d05","d06","d07","d08", "d09",
  "adm_sym_d01","adm_sym_d02","adm_sym_d03","adm_sym_d04","adm_sym_d05","adm_sym_d06","adm_sym_d07","adm_sym_d08","adm_sym_d09",
  "adm_sym_d01_d02","adm_sym_d01_d03","adm_sym_d01_d04","adm_sym_d01_d05","adm_sym_d01_d06","adm_sym_d01_d07","adm_sym_d01_d08","adm_sym_d01_d09"
) # surveillance strategies considered

##################
### LAUNCH !!! ###
##################

### put together plot names with associated lots
filepath_lot = paste0("Output_",lot, "/")
filepath_output = paste0(filepath_noel, "Surveillance/output/", filepath_lot)

# Set seed so stochastic results are reproducible
set.seed(1991)

# Loop through all factors
for(sens_function_x in v_sens_function_x){
  for(LTCF_x in v_LTCF_x){
    for(SIM_CTC in files_CTC){
      
      if(LTCF_x == 1){immunity_x = "0.20.2"; IPC_x = 1; Network_x = 1} 
      if(LTCF_x == 2){immunity_x = "0.20.2"; IPC_x = 1; Network_x = 2}
      if(LTCF_x == 3){immunity_x = "0.50.5"; IPC_x = 2; Network_x = 2}
      
      print(paste0('Executing surveillance loop \n Sensitivity ',sens_function_x,'; immunity ', immunity_x, '; IPC ', IPC_x, '; network ', Network_x, '; CTC simulation ', SIM_CTC))
      print(paste0("On CTC simulation ", SIM_CTC, " of ", last(files_CTC)))
      
      m_surveillance_results = f_surveillanceLoop(lot = lot_x, immunity = immunity_x, IPC = IPC_x, SIM = SIM_CTC, Network = Network_x,
                                                  testing_rounds = testing_rounds_x, screening_types = vec_screening_types, 
                                                  screening_targets = vec_screening_targets, surveillance_strategies = vec_surveillance_strategies,
                                                  sens_function = sens_function_x)
      
      write.csv(m_surveillance_results, file = paste0(filepath_output, "m_surveillance_",
                                                      immunity_x,
                                                      "_network", Network_x,
                                                      "_IPC", IPC_x,
                                                      "_SIM", SIM_CTC, "_SENS", sens_function_x,".csv"))
      # alternative for compressed file
      #write.csv(mtcars, file=gzfile("mtcars.csv.gz"))
    }
  }
}





