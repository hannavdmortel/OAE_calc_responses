##########################################################
## Unifying framework for assessing sensitivity of marine
## calcifiers to ocean alkalinity enhancement categorizes
## responses and identifies biological thresholds -
## importance of precautionary principle

## Contour plots for Fig. 1 (experimental data) and
## GLODAP data (Supplementary Fig. 1)
##########################################################
## Author: H. van de Mortel, G. Pelletier
## Version: 1.0
## Maintainer: H. van de Mortel
## Email: vandemortel.hanna@gmail.com
##########################################################

import PyCO2SYS as pyco2
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import PowerNorm

# Can consult PyCO2SYS manual:
# https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/

# Define filepath
filepath = "C:/path/to/your/folder/"

# Choose variables for CO2SYS calculations
sal = 34.68
temp = 16
silicate = 50  # umol/kg (= 30 uM)
phosphate = 0.5  # umol/kg (= 0.05 mg/L)

# Input parameters and constants for CO2SYS
parms = dict(
    salinity=sal,
    temperature=temp,
    total_silicate=silicate,
    total_phosphate=phosphate,
    opt_pH_scale=1, #total scale
    opt_k_carbonic=16 #Sulpis, 2020
    #Sulfuric acid dissociation = Dickson 1990
    #Boron dissociation = Uppstrom 1974
    )

# Functions to calculate pH, omega_ar & pCO2 using TA and 
# DIC and parameters defined above
# (Can change par1 and par2 types to calculate carb system 
# based on other parameters)
def calculate_pH(TA, DIC):
    results = pyco2.sys(
        **parms,
        par1=TA, par1_type=1,
        par2=DIC, par2_type=2
    )
    return results['pH_total']

def calculate_omega_ar(TA, DIC):
    results = pyco2.sys(
        **parms,
        par1=TA, par1_type=1,
        par2=DIC, par2_type=2
    )
    return results['saturation_aragonite']

def calculate_pCO2(TA, DIC):
    results = pyco2.sys(
        **parms,
        par1=TA, par1_type=1,
        par2=DIC, par2_type=2
    )
    return results['pCO2']

#%% Run for experimental data

# Import carbonate system data
# Requires TA and DIC 
df1 = pd.read_excel(filepath + "your_data_file.xlsx", 
                    sheet_name='Sheet1')

# Since we only want to use data for which we also have calc rate,
# salinity and temperature data, we filter out where this is missing
# (Can remove calc_rate_mask & temp_sal_mask if not relevant for you)

# Identify 'Calc rate' columns and create a mask for rows where
# any 'Calc rate' column has non-NaN values
calc_rate_cols = [col for col in df1.columns if 'Calc rate' in col]
calc_rate_mask = df1[calc_rate_cols].notna().any(axis=1)

# Create a mask for rows where Temp and Sal are not NaN
temp_sal_mask = df1[['Temp [°C]', 'Sal']].notna().all(axis=1)

# Create a mask for rows where 'TA' and 'DIC' are not NaN
ta_dic_mask = df1[['AT [µmol/kg] (compiled)', 
                   'DIC [µmol/kg] (compiled)']].notna().all(axis=1)

# Combine the masks and apply to filter df
final_mask = calc_rate_mask & ta_dic_mask & temp_sal_mask
filtered_df1 = df1[final_mask]

# Make sure values are numeric
filtered_df1['DIC'] = pd.to_numeric(filtered_df1['DIC [µmol/kg] (compiled)'], 
                                    errors='coerce')
filtered_df1['TA'] = pd.to_numeric(filtered_df1['AT [µmol/kg] (compiled)'], 
                                   errors='coerce')

# Save experimental DIC and TA data as x and y
x = filtered_df1['DIC']
y = filtered_df1['TA']

#%% Run for GLODAP data instead of experimental

# Import GLODAP data (0-50 m)
# (uncomment to replace definitions of df1, x and y)
# df1 = pd.read_excel(filepath+'GLODAP_DATA', sheet_name='0-50m')

# x = df1['DIC']
# y = df1['TA']

#%%
# Define limits of x and y-axes
xmin=1000
xmax=3200
ymin=1000
ymax=3200

# Create a grid over these specified limits
DIC_values = np.linspace(xmin, xmax, 100)
TA_values = np.linspace(ymin, ymax, 100)
DIC, TA = np.meshgrid(DIC_values, TA_values)

# Calculate variables for each grid point
pH = calculate_pH(TA, DIC)
omega_ar = calculate_omega_ar(TA, DIC)
pCO2 = calculate_pCO2(TA, DIC)

# Use PowerNorm to emphasize lower omega_ar values
gamma_omega_ar = 0.5  # Adjust gamma to emphasize lower values more
norm_omega_ar = PowerNorm(gamma=gamma_omega_ar, vmin=np.min(omega_ar), 
                          vmax=np.max(omega_ar))

# Use PowerNorm to emphasize lower pCO2 values
gamma_pCO2 = 0.5  # Adjust gamma to emphasize lower values more
norm_pCO2 = PowerNorm(gamma=gamma_pCO2, vmin=np.min(pCO2), vmax=np.max(pCO2))

#%% Create figure
fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=300)

# Scatter experimental data
for ax in axes:
    ax.scatter(x, y, c='xkcd:pink', alpha=0.5, zorder=3, edgecolor='none')
    # Set axes limits
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # Set x-axis label
    ax.set_xlabel('DIC [µmol/kg]', fontsize=15)

# Subplot 1 - pH
contour1 = axes[0].contourf(
    DIC, TA, pH, levels=np.arange(pH.min(), pH.max(), 0.01), cmap='viridis')
contour1.set_clim(6, 9);
contour1_lines = axes[0].contour(
    DIC, TA, pH, levels=np.array([7, 7.5, 8, 8.4, 8.6, 8.8, 9]), 
    colors='black', linewidths=2)
axes[0].clabel(contour1_lines, inline=True, fontsize=14, fmt='%1.1f')
axes[0].set_ylabel('TA [µmol/kg]', fontsize=15)
axes[0].set_title('$pH_{tot}$', fontsize=16)

# Subplot 2 - omega_ar
contour2 = axes[1].contourf(
    DIC, TA, omega_ar, levels=np.arange(omega_ar.min(), omega_ar.max(), 0.01), 
    cmap='viridis', norm=norm_omega_ar)
contour2.set_clim(0, 12);
contour2_lines = axes[1].contour(
    DIC, TA, omega_ar, levels=np.array([1,3,5,7,9,11]), 
    colors='black', linewidths=2)
axes[1].clabel(contour2_lines, inline=True, fontsize=14, fmt='%1.0f')
axes[1].set_title('$Ω_{ar}$', fontsize=16)

# Subplot 3 - pCO2
contour3 = axes[2].contourf(
    DIC, TA, pCO2, levels=np.arange(pCO2.min(), pCO2.max(), 2), 
    cmap='viridis_r', norm=norm_pCO2)
contour3.set_clim(0, 2000);
contour3_lines = axes[2].contour(
    DIC, TA, pCO2, levels=np.array([50, 100, 200, 500, 2000]), 
    colors='black', linewidths=2)
axes[2].clabel(contour3_lines, inline=True, fontsize=14, fmt='%1.0f')
axes[2].set_title('$pCO_{2}$', fontsize=16)

# Add NaOH and Na2CO3 lines
start_point1 = 2034  # DIC start point
start_point2 = 2303  # TA start point
end_point1 = 2034    # DIC end point
end_point2 = 3200    # TA end point

for ax in axes:
    # NaOH vertical line (DIC=1950)
    ax.plot([start_point1, end_point1], [start_point2, end_point2], 
            color='white', linewidth=4, zorder=10)
    ax.annotate('NaOH', xy=(start_point1-50, (start_point2 + end_point2) / 2), 
                color='black', fontsize=13, zorder=15,
                va='bottom', ha='right', rotation=90,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))

    # Add and annotate Na₂CO₃ line (starts at DIC=1950, TA=2200 and ends at TA=4000)
    end_point_na2co3 = 0.5 * (end_point2 - start_point2) + start_point1
    ax.plot([start_point1, end_point_na2co3], [start_point2, end_point2], color='white', linewidth=4, linestyle='--', zorder=10)
    ax.annotate('Na₂CO₃', xy=(((start_point1 + end_point_na2co3) / 2) + 100, 
                              (start_point2 + end_point2) / 2), 
                color='black', fontsize=13, zorder=15, va='bottom', ha='right', 
                rotation=55, #adjust rotation as needed
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))


# Annotate plots
axes[0].annotate('a.', xy=(3050, 1150), color='white', fontsize=21)
axes[0].annotate('pH<7', xy=(2500, 1600), color='white', fontsize=21)
axes[1].annotate('b.', xy=(3050, 1150), color='white', fontsize=21)
axes[1].annotate('$Ω_{ar}$<1', xy=(2500, 1600), color='white', fontsize=21)
axes[2].annotate('c.', xy=(3050, 1150), color='white', fontsize=21)
axes[2].annotate('$pCO_{2}$>2000', xy=(2350, 1600), color='white', fontsize=21)

# Save figure
plt.tight_layout()
plt.savefig('C:/insert/figure/location/here.png', bbox_inches='tight')
plt.show()
