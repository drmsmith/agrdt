# agrdt
An evaluation of rapid diagnostic testing (Ag-RDT) for SARS-CoV-2 transmission prevention in long-term care.

This R code supports Smith DRM et al., (2021). Rapid antigen testing as a reactive response to surges in nosocomial SARS-CoV-2 outbreak risk.

Available at: https://www.nature.com/articles/s41467-021-27845-w

All code was developed, tested and run using Rv3.6.0

# about
Ouptut data from stochastic outbreak simulations using CTCmodeler are provided. Epidemiological outcomes from these data are calculated. Surveillance interventions are simulated, applied retrospectively to the outbreak data. Summary outcomes from simulated surveillance are calculated, and plots are rendered. Correspondings files and folders are described below. We recommend working through them in order to understand the data and successive steps of operations underlying simulations and outcomes.

# Rproject
* All files are associated with an R project (for stable working directory)
  * agrdt.Rproj

# CTCmodeler outbreak simulations
* (1) Model_herdI_incidenceLow.zip:
  *  Folder containing output data from outbreak simulations. Sub-folders contain corresponding types of output data, for all outbreak simulations across 3 included LTCFs. For simplicity, only data for the baseline "low community incidence" scenario are provided; throughout, this scenario and its corresponding lot of simulations are referred to as "lot3"
  *  Guide to interpret filenames, using e.g. CommCases0.20.2Christmasnetwork1IPC1low_SIM_1_Bacteria_Coronaviridae_Alphaletovirus_Covid19___)
		* CommCases: type of output data saved herein (others are FirstIndex, ManuportageTransmission, PEreplacement, secondaryCases, statusByDay)
		* 0.20.2: baseline 20% immunization rate in patients, staff (alternative: 0.50.5)
		* network1: baseline contact network (alternative: network2 is social distancing intervention)
		* IPC1: baseline IPC rate, 0.36 (alternative: IPC2 is with face masks, 0.8) 
		* low: baseline low incidence scenario
		* SIM_1: stochastic outbreak ID (range 1:100)
	* (1.1) CommCases/: community-onset "introductions"
	* (1.2) FirstIndex/: community-onset "index cases" 
	* (1.3) ManuportageTransmission/: cases of vectors acquiring transient SARS-CoV-2 carriage
	* (1.4) PEreplacement/: log of staff sick leave
	* (1.5) secondaryCases/: transmission chains
	* (1.6) statusByDay/: daily SARS-CoV-2 infection status for all individuals

# functions
* (2) functions.R
  * contains all functions for analysis and surveillance
 * (2.1) toyadmission.csv: a file describing all patient admissions to the LTCF, needed for functions that simulate testing of newly admitted individuals  

# diagnostic sensitivity curves
* (3.1) kucirka_raw.csv: downloaded from https://www.acpjournals.org/doi/suppl/10.7326/M20-1495
* (3.2) sensitivity_curves.R: using kucirka_raw.csv, generates the 3 alternative sensitivity curves listed below, as described in article
* (3.3) kucirka_adjusted.csv: sensitivity curves for main analysis 
* (3.4) kucirka_perfect.csv: sensitivity curves for sensitivity analysis (perfect sensitivity)
* (3.5) kucirka_uniform.csv: sensitivity curves for sensitivity analysis (uniform sensitivity)

# analysis_CTC
* (4) analysis_CTC/: a folder containing analysis of epidemiological data from CTC simulations (prior to surveillance simulation)
  * (4.1) incidence.R: calculates incidence, including contribution of super-spreading to transmission and acquisition
  * (4.2) prevalence.R: calculates prevalence
  * (4.3) sensitivity_dynamics.R: evaluate dynamics of true-positive rate, combining infection data and sensitivity curves
  * (4.4) data_analysis_CTC/lot3/: corresponding CTC analysis outcomes are saved here 

# surveillance loop
* (5) surveillance_loop.R: R file that executes functions from functions.R upon files from Model_herdI_incidenceLow to simulate surveillance
  * results are saved to Surveillance/output/Output_lot3/
  * follows naming convention for CTC files: files for LTCF1 end in 0202_1_1, LTCF2 in 0202_2_1, and LTCF3 in 0505_2_2
  * NB: this program was originally designed for parallel dispatch on a computing cluster; to run all simulations considered and ultimately used to produce results (below), estimated run-time on a normal desktop computer is roughly 100 hours. See file to modify initialization conditions to run small demo runs as desired.

# raw surveillance output
* (6) Surveillance/output/Output_lot3/: folder containing CSV output files from surveillance loop
  * each corresponds to a different CTC outbreak and sensitivity assumption (time-varying, uniform, or perfect)
  * due to space limitations, simulated surveillance data provided here are limited to the first 5 CTC simulations for each LTCF
  * however, if desired, surveillance for all outbreaks can be reproduced locally using surveillance loop as described above

# processing surveillance output
* (7) process_surveillance_data_lot3.R: a file with two key functions
  * (1) first, it reads and consolidates raw surveillance output from Surveillance/output/Output_lot3/: 
  * for each LTCF, all surveillance results are combined into one dataframe; these can be saved into a folder called Surveillance/output/Output_prepared_lot3/ (not possible here; files are roughly 50-75MB, and exceed present space limit)
  * (2) second, it calculates surveillance outcomes: 
  * by executing corresponding functions from functions.R upon the "prepared" data from Surveillance/output/Output_prepared_lot3/ and saving results into the folder Surveillance/output/Outcomes_lot3/
  * NB: although final "prepared data" could not fit here, all final outcomes files are provided, corresponding to the 3 LTCFs, and hence all results and figures for the baseline low incidence scenario presented in the manuscript can be reproduced

# figures
* (8) plots/: a folder containing two R files to render plots
  * (8.1) paper_figure_aesthetics.R: a file with labels and aesthetic arguments for plots
  * (8.2) paper_figures.R: a file to render figures, using data from data_analysis_CTC/lot3/ and Surveillance/output/Outcomes_lot3/

# contact
David Smith \
Institut Pasteur / Inserm / UVSQ \
david.smith@pasteur.fr \
davidrobertmundysmith@gmail.com
