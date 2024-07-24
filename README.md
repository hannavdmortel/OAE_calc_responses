# OAE_calc_responses
'Unifying framework for assessing sensitivity of marine calcifiers to ocean alkalinity enhancement categorizes responses and identifies biological thresholds - importance of precautionary principle'
Authors: Nina Bednaršek, Hanna van de Mortel, Greg Pelletier, Marisol García-Reyes, Richard A. Feely, Andrew G. Dickson

Correspondence regarding code: vandemortel.hanna@gmail.com

OAE_calc_responses includes the code publicly available from this paper, separated into two folders: Contour_plots and TA_addition_and_regressions.

Contour_plots allows you to recreate figure 1 of the paper as well as the supplementary figure 1. It is also possible to add alkalinity sources other than NaOH and Na2CO3 if the effect on TA and DIC is known. 

The folder TA_addition_and_regressions contains two scripts. 
They both import a file (df1) containing at least a column with species names, calcification rate and two carbonate system parameters. 
Naming of columns and variables used in increased_TA_carb_system and plots_per_species can be changed accordingly. 
increased_TA_carb_system should be run before plots_per_species. The former computes the current condition baselines, conceptually adds alkalinity in the form of NaOH and Na2CO3 from this baseline and computes at what alkalinity addition a pH of 9 is exceeded.
plots_per_species calculating the calcification rate vs TA:DIC regressions and visualizes them, as well as calculating the species-specific NaOH and Na2CO3 thresholds & inflection points. 
