### all the necessities for plotting aesthetics and formatting

#cols_classes = c('black','#E31B23', '#fdae61', '#005CAB')
cols_classes = c('black', '#005CAB', '#E31B23', '#fdae61')
cols_types = c('black', '#008837', '#7b3294')
cols_onset = c("#00BFC4", "#F8766D")
cols_infection_status = c('#d7191c', '#fdae61',"#FFFFBF", '#abdda4', '#2b83ba')
cols_infection_status_blue = c('#bdc9e1','#74a9cf','#2b8cbe','#045a8d')
cols_basePCR = c("#B16243FF", "#5CC8D7FF")
cols_bootstrap = rev(c("#542788", "#e08214", "black"))
cols_typescreening = c('#c51b7d', '#4d9221', 'black')
cols_sensfuncs = c("#a6611a",'black', "#018571")
cols_uncertainty = c("#e08214", "#8073ac", "black")
cols_implementation = c('#d95f02', '#984ea3', 'black')
cols_rdt_pcr = c('red','black')
cols_ltcf = c('#7fcdbb','#1d91c0', '#253494')
cols_strats3 = c("#ffa600", "#ff6361", "darkgrey", "#58508d")

shapes_cat = c(17,16,8)

vec_3ltcfs_labels = c('LTCF 1 (high risk)', 'LTCF 2 (moderate risk)', 'LTCF 3 (low risk)')
vec_3ltcfs_labels_control = c('LTCF 1 (low control)', 'LTCF 2 (moderate control)', 'LTCF 3 (high control)')

all_strats_3ltcfs = c("adm", "sym", "adm_sym", 
                      "d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09",
                      "adm_sym_d01", "adm_sym_d02", "adm_sym_d03", "adm_sym_d04", "adm_sym_d05",
                      "adm_sym_d06", "adm_sym_d07", "adm_sym_d08", "adm_sym_d09",
                      "adm_sym_d01_d02", "adm_sym_d01_d03", "adm_sym_d01_d04", "adm_sym_d01_d05", 
                      "adm_sym_d01_d06", "adm_sym_d01_d07", "adm_sym_d01_d08", "adm_sym_d01_d09")
all_strats_3ltcfs_labels = c("Admissions (routine PCR)", "Symptoms (routine PCR)",'Baseline (routine PCR)', 
                             paste0('Only Ag-RDT (day ',1:9, ')'),
                             paste0('Ag-RDT (day ', 1:9, ')'),
                             paste0('Ag-RDT (days 1 & ',2:9,')'))
all_strats_3ltcfs_labels_UPDATED = c("alternative RT-PCR 2 (only admissions)", "alternative RT-PCR 1 (only symptoms)",'routine RT-PCR (admissions + symptoms)', 
                                     paste0('Ag-RDT screening (day ',1:9, ')'),
                                     paste0('routine RT-PCR + Ag-RDT screening (day ', 1:9, ')'),
                                     paste0('routine RT-PCR + Ag-RDT screening (days 1 & ',2:9,')'))


### labels for super-spreading categories
vec_superspreading_labels = c('low spreader', 'super spreader', 'non spreader')
vec_superspreading_labels_long = c('low spreader (1-2 individuals)', 'super spreader (≥3 individuals)', 'non spreader (0 individuals)')
