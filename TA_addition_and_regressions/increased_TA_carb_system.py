##########################################################
## Unifying framework for assessing sensitivity of marine
## calcifiers to ocean alkalinity enhancement categorizes
## responses and identifies biological thresholds -
## importance of precautionary principle

## Compute current condition baselines, conceptually add
## alkalinity (in the form of NaOH and Na2CO3) and compute
## at what alkalinity addition a pH of 9 is exceeded
##########################################################
## Author: H. van de Mortel, M. Garcia-Reyes
## Version: 1.0
## Maintainer: H. van de Mortel
## Email: vandemortel.hanna@gmail.com
##########################################################

import pandas as pd
import numpy as np
import PyCO2SYS as pyco2
import re
from scipy.stats import zscore

# Can consult PyCO2SYS manual:
# https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/


filepath = "C:/path/to/your/folder/"

# Import data:
# Should include (at minimum) two columns with carb parameters to solve 
# the carbonate system (in this case TA and DIC), a column with species names 
# and a column with calcification rate data. If species and/or rate units
# are all uniform, the for-loop can be simplified.
df1 = pd.read_excel(filepath + "your_data_file.xlsx", 
                    sheet_name='Sheet1')

#%% Variables and grouping

#Group df by species
grouped = df1.groupby('Species (compiled)')

# Define different rate units, so they can be handled separately in for-loop
# (These are all separate column headers in df1; change accordingly)
rate1 = 'Calc rate [mmol/g/hr] (compiled)'
rate2 = 'Calc rate [mmol/m**2/hr] (compiled)'
rate3 = 'Calc rate [mmol/m**3/hr] (compiled)'
rate4 = 'Calc rate [mmol/#/hr] (compiled)'
rate5 = 'Calc rate [mmol/hr] (compiled)'
rate6 = 'Calc rate [mmol/#] (compiled)'
rate7 = 'Calc rate [mmol/cm**2] (compiled)'
rate8 = 'Calc rate [%] (compiled)'
rate9 = 'Calc rate [mg] (compiled)'
rate10 = 'Calc rate [%/hr] (compiled)'
rate11 = 'Calc rate [#/hr] (compiled)'
rate12 = 'Calc rate [1/hr] (compiled)'
rate13 = 'Calc rate (shell length) [mm] (compiled)'
rate14 = 'Calc rate (shell thickness) [µm] (compiled)'

# Create list for for-loop
rates = [rate1, rate2, rate3, rate4, rate5, rate6, rate7, rate8, 
         rate9, rate10, rate11, rate12, rate13, rate14]

#%% For-loop
# Generates baselines for each species, adds alkalinity from the baseline
# and computes at which alkalinity addition a pH of 9 is exceeded

# Create an empty list to store DataFrames
results_list = []

#Estimation silicate and phosphate concentrations (minor influence)
silicate = 50 #umol/kg (= 30 uM)
phosphate = 0.5 #umol/kg (= 0.05 mg/L)

# Variable to track the increment at which pH exceeds 9
increment_exceeds_pH9_1 = None
increment_exceeds_pH9_2 = None
found_exceed_pH9_1 = False
found_exceed_pH9_2 = False

# Loop through species
for species, group in grouped:
    # Loop through different rate units
    for rate in rates:          
        if rate in group.columns and not group[rate].dropna().empty:
            # Remove nan values and make sure data is all numeric
            rate_group = group[group[rate].notna()]
            rate_group.loc[:, rate] = pd.to_numeric(rate_group[rate], errors='coerce')

            # Define parameters for PyCO2SYS
            # Can insert fixed values for Sal and Temp, or use specific 
            # Sal and Temp data
            parms = dict(salinity=rate_group['Sal'].mean(), temperature = 20,
                         total_silicate=silicate, total_phosphate=phosphate,
                         opt_pH_scale=1, #total scale
                         opt_k_carbonic=16 #Sulpis, 2020
                         #Sulfuric acid dissociation = Dickson 1990
                         #Boron dissociation = Uppstrom 1974
                         )
            
            #Compute current baseline values with avg sal for this species 
            # and T = 20, using a pH of 8.1 and pCO2 of 425 ppm
            baseline = pyco2.sys(**parms, 
                                par1 = 8.1, par1_type = 3,
                                par2 = 425, par2_type = 4)
            
            # Store baseline values as variables
            TA_baseline = baseline['alkalinity']
            DIC_baseline = baseline['dic']
            pH_baseline = baseline['pH']
            
            # Save variables in results_dict
            results_dict = {
            'Group': rate_group['Group'].iloc[0],
            'Species': species,
            'Rate': rate,
            'Baseline': TA_baseline/DIC_baseline,
            'TA (check)': TA_baseline,
            'DIC (check)': DIC_baseline,
            'pH (check)': pH_baseline}
        
            # Add TA in steps (up to 10,000 umol/kg) as NaOH and Na2CO3
            # (Can also add in the form of another alkalinity source, if
            # the effect on TA:DIC is known)
            for i in [10] + list(range(50, 10001, 50)):
                #########NaOH#########
                # Compute carbonate system for each addition of NaOH
                results_NaOH = pyco2.sys(**parms, 
                                    par1 = TA_baseline+i, par1_type = 1,
                                    par2 = DIC_baseline, par2_type = 2)

                # Populate results_dict with carb system variables of interest
                # Define carb system variables of interest and corresponding 
                # result fields
                fields1 = [
                    ('TA:DIC', 'alkalinity', lambda r: r['alkalinity'] / r['dic']),
                    ('pH', 'pH', lambda r: r['pH']),
                    ('Ω', 'saturation_aragonite', lambda r: r['saturation_aragonite']),
                    ('H+ (free)', 'hydrogen_free', lambda r: r['hydrogen_free']),
                    ('HCO3', 'bicarbonate', lambda r: r['bicarbonate']),
                    ('CO3', 'carbonate', lambda r: r['carbonate']),
                    ('CO2', 'aqueous_CO2', lambda r: r['aqueous_CO2']),
                    ('pCO2', 'pCO2', lambda r: r['pCO2'])
                ]

                # Populate results_dict using a loop
                for field_suffix, result_key, result_extractor in fields1:
                    column_name = f'NaOH {field_suffix} +{i} umol/kg'
                    results_dict[column_name] = result_extractor(results_NaOH)

                #########Na2CO3#########
                results_Na2CO3 = pyco2.sys(**parms, 
                                par1 = TA_baseline+i, par1_type = 1,
                                par2 = DIC_baseline+(0.5*i), par2_type = 2)
                # Populate results_dict with carb system variables of interest
                # Define carb system variables of interest and corresponding 
                # result fields
                fields2 = [
                    ('TA:DIC', 'alkalinity', lambda r: r['alkalinity'] / r['dic']),
                    ('pH', 'pH', lambda r: r['pH']),
                    ('Ω', 'saturation_aragonite', lambda r: r['saturation_aragonite']),
                    ('H+ (free)', 'hydrogen_free', lambda r: r['hydrogen_free']),
                    ('HCO3', 'bicarbonate', lambda r: r['bicarbonate']),
                    ('CO3', 'carbonate', lambda r: r['carbonate']),
                    ('CO2', 'aqueous_CO2', lambda r: r['aqueous_CO2']),
                    ('pCO2', 'pCO2', lambda r: r['pCO2'])
                ]

                # Populate results_dict using a loop
                for field_suffix, result_key, result_extractor in fields2:
                    column_name = f'Na2CO3 {field_suffix} +{i} umol/kg'
                    results_dict[column_name] = result_extractor(results_Na2CO3)
                    
                # Compute at which TA additions pH exceeds 9
                if results_NaOH['pH'] > 9:
                    if found_exceed_pH9_1 == False:
                        increment_exceeds_pH9_1 = i
                        found_exceed_pH9_1 = True
                        break
                
                if results_Na2CO3['pH'] > 9:
                    if found_exceed_pH9_2 == False:
                        increment_exceeds_pH9_2 = i
                        found_exceed_pH9_2 = True
                        break
            
            # Add column for TA increment pH 9 exceedance
            results_dict['+NaOH [umol/kg] > pH 9'] = increment_exceeds_pH9_1
            results_dict['+Na2CO3 [umol/kg] > pH 9'] = increment_exceeds_pH9_2

            # Append results into results_dict
            results_list.append(results_dict)   

# Concatenate all DataFrames in the list into a single DataFrame
final_results = pd.DataFrame(results_list)
final_results = final_results.sort_values(by=['Group', 'Species'], ascending=[True, True])

# Save df to Excel file
output_file_path = filepath+"name_file.xlsx"
final_results.to_excel(output_file_path, index=False)
