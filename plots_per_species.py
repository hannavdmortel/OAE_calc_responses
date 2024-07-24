##########################################################
## Unifying framework for assessing sensitivity of marine
## calcifiers to ocean alkalinity enhancement categorizes
## responses and identifies biological thresholds -
## importance of precautionary principle

## Calculating calc rate vs TA:DIC regressions and
## visualizing them, as well as calculating the species-
## specific NaOH and Na2CO3 thresholds & inflection points
##########################################################
## Author: H. van de Mortel, G. Pelletier, M. Garcia-Reyes
## (Functions by G. Pelletier)
## Version: 1.0
## Maintainer: H. van de Mortel
## Email: vandemortel.hanna@gmail.com
##########################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures
from scipy.stats import zscore
from numpy import exp, linspace
import lmfit as fit
from lmfit.models import ExpressionModel
from scipy import stats
from statsmodels.tools.eval_measures import rmse
import re

#%% Import data
filepath = "C:/path/to/your/folder/"

# Import data:
# Should include (at minimum) two columns with carb parameters to solve 
# the carbonate system (in this case TA and DIC), a column with species names 
# and a column with calcification rate data. If species and/or rate units
# are all uniform, the for-loop can be simplified.
df1 = pd.read_excel(filepath + "your_data_file.xlsx", 
                    sheet_name='Sheet1')

#Import added alkalinity data
df2 = pd.read_excel(filepath + "output_from_increased_TA_carb_system.xlsx", 
                    sheet_name='Sheet1')

#Choose until what value you want to run plots and compute thresholds
max_thresh = 501 #10001

#Define variables for simplicity
pH = 'pH_tot (compiled)'
TA = 'AT [µmol/kg] (compiled)'
DIC = 'DIC [µmol/kg] (compiled)'
species2 = 'Species (compiled)'

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

# Create list of rates
rates = [rate1, rate2, rate3, rate4, rate5, rate6, rate7, rate8, 
         rate9, rate10, rate11, rate12, rate13, rate14]

# Make sure all rate columns are numeric
for rate in rates:
    df1[rate] = pd.to_numeric(df1[rate], errors='coerce')

#Convert numeric columns to numeric type and handle non-numeric values
df1[pH] = pd.to_numeric(df1[pH], errors='coerce')
df1[DIC] = pd.to_numeric(df1[DIC], errors='coerce')
df1[TA] = pd.to_numeric(df1[TA], errors='coerce')

#Calculate TA:DIC column
nan_mask = np.logical_or(np.isnan(df1[DIC]), np.isnan(df1[TA]))
df1['TA:DIC'] = np.where(nan_mask, np.nan, df1[TA] / df1[DIC])

#Group by species for plotting
grouped = df1.groupby(species2)

#Define a colormap and the number of unique colors for plotting from different studies
num_unique_colors = 10
cmap = plt.get_cmap('Accent')

#Choose colors for +TA data
cmap_NaOH = ['#810f7c', '#8856a7', '#8c96c6', '#9ebcda', '#bfd3e6', '#edf8fb','#810f7c', '#8856a7', '#8c96c6', '#9ebcda', '#bfd3e6', '#edf8fb']
cmap_Na2CO3 = ['#253494', '#2c7fb8', '#41b6c4', '#7fcdbb', '#c7e9b4', '#ffffcc', '#253494', '#2c7fb8', '#41b6c4', '#7fcdbb', '#c7e9b4', '#ffffcc']

#Create lists for variables to be filled later
s_species, s_rate, pqs, pls, s_iyval, pes = list(), list(), list(), list(), list(), list()
perc = [10,25,50,75,90]
#%% Make lists for best fits dataframe
best_fits_studies = [] #studies
best_fits_n = [] #number of datapoints
best_fits_groups = [] #functional group
best_fits_species = [] #species
best_fits_rate = [] #rate unit
best_fits = [] #response (linear, parabolic or exponential)
best_fits_p = [] #p-value
best_fits_r2 = [] #R^2
best_fits_rmse = [] #RMSE
best_fits_pre_industrial= [] #pre-industrial calcification rate
best_fits_current = [] #current calcification rate (baseline)
best_fits_pH_9 = [] #

#Make lists for inflection points dataframe (for parabolic responders)
inflection_species = []
inflection_rate = []
inflection_point = []

#Make lists for thresholds
NaOH_threshold_species = []
NaOH_threshold_groups = []
NaOH_threshold_rate = []
NaOH_threshold_amount = []
NaOH_thresholds = []

Na2CO3_threshold_species = []
Na2CO3_threshold_groups = []
Na2CO3_threshold_rate = []
Na2CO3_threshold_amount = []
Na2CO3_thresholds = []

#%% Functions

# Function to use lmfit to find the best-fit parameters of the 3-parameter 
# exponential regression for threshold+ or threshold- response
def fit_exp3(x,y,  nper):
    # - - -
    # function to calculate the best-fit nonlinear regression parameters 
    # of the 3-parameter exponential regression for threshold+ or threshold- response
    # using lmfit (https://lmfit.github.io//lmfit-py/index.html)
    # Greg Pelletier (gjpelletier@gmail.com)
    # - - -

    # INPUT
    # x = vector of observed x
    # y = vector ov observed y
    #
    # OUTPUT
    # param = vector of best-fit parameters of the nonlinear regression function

    # find y stats to use for intial parameter estimates
    ymin = np.nanmin(y)
    yn = np.percentile(y,nper)
    ymed = np.median(y)
    ymax = np.nanmax(y,0)

    mod = ExpressionModel('b1 + b2 * exp(b3 * x)')
    pars = fit.Parameters()

    # - - -
    # This is where you set the initial values of the regression parameters b1, 
    # b2, b3 in pars.
    # The initial guess of parameters set to ypred = b1 = ymed (np.median(y)) 
    # usually is the best guess for intial b1
    # Initial values of b2=0 and b3=0 usually works best
    # In other words, the following initial values of b1, b2, b3 as follows 
    # usually works:
    # b1 = ymed
    # b2 = 0
    # b3 = 0
    # Sometimes this does not work. 
    # You can tell that it does not work when the best fit is a linear straight 
    # line through the observations instead of a curving exponential line.
    # In that case, try using either b1=y25th or b1=y75th instead of b1=med
    # Usually one of those for initial b1 will find a better fit with a curving 
    # exponential thrshold line when the b1=ymed does not work
    # In this example, b1=ymed does not work, but b1=y25th does work
    # Uncomment one of the following three pars['b1'] lines. Try ymed first. 
    # If that does not work then try y25th and y75th, and pick whichever works best
    # - - -
    #if yb=='ymed':
    #    pars['b1'] = fit.Parameter(name='b1', value=ymed, min=-np.inf, max=np.inf)
    #elif yb=='y25th':
    #    pars['b1'] = fit.Parameter(name='b1', value=y25th, min=-np.inf, max=np.inf)
    #else:
    #    pars['b1'] = fit.Parameter(name='b1', value=y75th, min=-np.inf, max=np.inf)
        
    pars['b1'] = fit.Parameter(name='b1', value=yn, min=-np.inf, max=np.inf)
    pars['b2'] = fit.Parameter(name='b2', value=0, min=-np.inf, max=np.inf)
    pars['b3'] = fit.Parameter(name='b3', value=0, min=-np.inf, max=np.inf)

    # extract parameter values from the best fit output and make user-defined function for the regression
    out = mod.fit(y, pars, x=x)
    d = out.params
    b1 = d['b1'].value
    b2 = d['b2'].value 
    b3 = d['b3'].value 
    param = [b1, b2, b3]
    
    return out, param


# function to calculate the confidence interval or prediction interval
def prediction_interval(f,param,xpred,x,y,alpha,prediction):
    # - - -
    # function to calculate the confidence interval or prediction interval
    # for any user-defined regression function.
    # Greg Pelletier (gjpelletier@gmail.com)
    # - - -

    # INPUT
    # f = user-defined regression @ function to predict y given inputs if param 
    # and x values (xpred and x)
    # 	For example, if using the 3-parameter nonlinear regression exponential 
    #   threshold function, then f = lambda param,xval : param[0] + param[1] * 
    #   exp(param[2] * xval)
    # param = vector of best-fit parameters of the regression function
    # xpred = vector of x values to evaluate predicted y values 
    # (e.g. xpred=linspace(min(x),max(x))
    # x = vector of observed x
    # y = vector ov observed y
    # alpha = probability value for the prediction interval (e.g. alpha=0.05 is the 95# prediction interval)
    # prediction = True or False, where True= prediction interval of additional measurements, and False= confidence iterval of original observations
    #
    # OUTPUT
    # ypred = predicted y values at xpred
    # ciLo = lower prediction interval or confidence interval for each value in xpred
    # ciUp = upper prediction interval or confidence interval for each value in xpred
    # p-value of the regression F-test comparing MSregression/MSresidual

    # use 2-tailed t statistic
    pLo = alpha/2
    pUp = 1-alpha/2
    nobs = np.size(x)
    nu = nobs-np.size(param)    # degrees of freedom = number of samples - number of parameters

    # critical t values for the upper and lower CI
    tcrit = stats.t.ppf([pLo, pUp], nu)
    
    # predicted ypred at xpred and yhat at x
    ypred = f(param,xpred)
    yhat = f(param,x)
    
    # residual sum of squares
    SSresidual = 0
    for i in range(nobs):
        SSresidual = SSresidual + (y[i] - yhat[i]) ** 2
    # residual mean square
    MSresidual = SSresidual / (np.size(x)-np.size(param))
    
    # syx = standard error of the estimate 
    syx = np.sqrt(MSresidual)

    # mean value of observed x
    xmean = np.mean(x)

    # sum of squares for x
    SSx = np.sum(x **2) - np.sum(x) **2 / nobs

    # calculate F-statistic and p-value of the regression
    SStotal = np.sum(y **2) - np.sum(y) **2 / nobs
    SSregression = SStotal - SSresidual
    # MSregression = SSregression
    MSregression = SSregression / (np.size(param)-1)
    Fstat = MSregression / MSresidual
    dfn = np.size(param) - 1                # df numerator = degrees of freedom for model = number of model parameters - 1
    dfd = np.size(x) - np.size(param)       # df denomenator = degrees of freedom of the residual = nobs - nparam
    p = 1-stats.f.cdf(Fstat, dfn, dfd)      # p-value of F test statistic 

    # calculate rsquared (ordinary and adjusted)
    rsquared = SSregression / SStotal                                           # ordinary rsquared
    adj_rsquared = 1-(1-rsquared)*(np.size(x)-1)/(np.size(x)-np.size(param)-1)  # adjusted rsquared

    # calculate sqrterm
    npred = np.size(xpred)
    sqrterm = np.empty(npred)
    for i in range(npred):
        if prediction:
            # prediction interval for additional observations
            sqrterm[i] = np.sqrt(1 + 1/nobs + (xpred[i] - xmean)**2 / SSx)
        else:
            # confidence interval of original observations
            sqrterm[i] = np.sqrt(1/nobs + (xpred[i] - xmean)**2 / SSx)

    ciLo = ypred + tcrit[0] * syx * sqrterm
    ciUp = ypred + tcrit[1] * syx * sqrterm

    return ypred, ciLo, ciUp, p, rsquared

# - - -
# function to calculate the adjusted r-squared
def adjusted_rsquared(rsquared,n,p):
    # calculate the adjusted rsquared from input of the following:
    # rsquared
    # n = number of samples
    # p = number of parameters
    adj_rsquared = 1-(1-rsquared)*(n-1)/(n-p-1)
    return adj_rsquared

#%% For-loop
# Loop through species
for species, group in grouped:
    # Loop through different rate units
    for rate in rates:                    
        rate_group = group[group[rate].notna()]
        rate_group[rate] = pd.to_numeric(rate_group[rate], errors='coerce')

        # Determine the number of unique 'File Name' values in the current species group
        unique_files = rate_group['Study'].nunique()
        unique_file_names = ''  # Initialize an empty string to hold combined study names

        # Generate a list of colors based on the number of unique 'File Name' values
        colors = [cmap(i) for i in np.linspace(0, 1, num_unique_colors)[:unique_files]]
        
        #Create mask to identify NaN
        x_ = rate_group['TA:DIC']
        y_ = rate_group[rate]
        nan_mask = np.logical_or(pd.isnull(x_), pd.isnull(y_))
        x = x_[~nan_mask]
        y = y_[~nan_mask]
            
        # Check if the data contains valid values to perform regressions
        if len(y) > 3:
            #Remove outliers (up to 6 std deviations)
            z_scores = zscore(y)
            z_threshold = 6
            outlier_mask = np.abs(z_scores) < z_threshold
            x = x[outlier_mask]
            y = y[outlier_mask]

            #Filter df2 based on species
            df2_ = df2[(df2['Species'] == species) & (df2['Rate'] == rate)]

            ##### Import added TA values #####
            # Define the column prefixes and suffixes
            prefixes = ['NaOH TA:DIC +', 'Na2CO3 TA:DIC +']
            suffix = ' umol/kg'
            
            # Create the array for NaOH and Na2CO3 TA:DIC up until 10,000 umol/kg
            # (This is for calculating thresholds)
            naoh_columns = [df2_[f'{prefixes[0]}{i}{suffix}'] for i in [10] + list(range(50, max_thresh, 50))]
            naco3_columns = [df2_[f'{prefixes[1]}{i}{suffix}'] for i in [10] + list(range(50, max_thresh, 50))]
      
            # Concatenate both arrays into x1
            x1 = np.array(naoh_columns + naco3_columns)
            
            # Create labels for TA added
            labels =[10] + list(range(50, max_thresh, 50))
            
            # Append study names
            for i, study in enumerate(rate_group['Study'][~nan_mask].unique()):
                if i > 0:
                     unique_file_names += ', '  
                unique_file_names += study
             
            # Convert x to array and create linspace for regressions
            x_array = x.to_numpy()
            xpred_l = linspace(x.min(), x.max(), 100)
            
            # Create figure
            fig = plt.figure(figsize=(6,5),dpi=200)

            ##### REGRESSIONS #####
            # Define models and stats
            
            ### Linear ###
            polynomial_features = PolynomialFeatures(degree=1)
            xp = polynomial_features.fit_transform(x_array.reshape(-1, 1))
            model_lin = sm.OLS(y, xp).fit()
            lp=model_lin.f_pvalue
            lrsqu=model_lin.rsquared
            lparams = model_lin.params
            fl = lambda lparams, xval : lparams[0] + lparams[1]*xval 
            ypred_lin = fl(lparams, x)
            lrmse = rmse(y, ypred_lin)
            
            ### Quadratic ###
            polynomial_features = PolynomialFeatures(degree=2)
            xp = polynomial_features.fit_transform(x_array.reshape(-1, 1))
            model_quad = sm.OLS(y, xp).fit()
            qp=model_quad.f_pvalue
            qrsqu=model_quad.rsquared
            qparams = model_quad.params
            fq = lambda qparams, xval : qparams[0] + qparams[1]*xval + qparams[2]*xval**2
            ypred_quad = fq(qparams, x)
            qrmse = rmse(y, ypred_quad)
            
            ### Exponential ###
            # Uses function due to statsmodels exponential model not 
            # functioning correctly for our data
            #s_func.append()
            ps, rsqs = list(), list()
            for j in perc:
                out, param = fit_exp3(x,y, j)
                #print(out.fit_report())

                f = lambda param,xval : param[0] + param[1] * exp(param[2] * xval)

                # use function to find ypred, ciLo, ciUp, and p-value
                # on the range of the data
                #xpred_l = linspace(x.min(), x.max(), 100)
                ypred_l, ciLo, ciUp, p, rsq = prediction_interval(f,param,xpred_l,x.values,y.values,.1,True)
                #print('p-value= ',p)
                ps.append(p)
                rsqs.append(rsq)
            xp = np.min(ps)
            mi = ps.index(xp)
            xqsqu = rsqs[mi]
            
            # recalculate exponential
            out, xparam = fit_exp3(x,y, perc[mi])
            #print(out.fit_report())
            fe = lambda xparam,xval : xparam[0] + xparam[1] * exp(xparam[2] * xval)
            # use function to find ypred, ciLo, ciUp, and p-value
            # on the range of the data
            ypred_exp = fe(xparam, x)
            ypred_exp_l, ciLo, ciUp, xp, rsq = prediction_interval(fe,xparam,xpred_l,x.values,y.values,.1,True)
            xrmse = rmse(y, ypred_exp)
                
            # Define p-value threshold
            p_thr = 0.05
            
            # If any of the model p-values are below the threshold:
            if (lp<=p_thr) | (qp<=p_thr) | (xp<=p_thr):                
                studies_handles = []
                #Plot the data points with different colors for each 'File Name' value
                for i, study in enumerate(rate_group['Study'][~nan_mask].unique()):
                    mask = rate_group['Study'] == study
                    studies = plt.scatter(x[mask], y[mask], label=study, 
                                          color=colors[i], alpha=0.7, edgecolor='none')
                    studies_handles.append(study)

                #######NaOH#######
                x_NaOH = x1[0:len(naoh_columns)]
                #######Na2CO3#######
                x_Na2CO3 = x1[len(naoh_columns):len(naoh_columns)*2]

                # x ranges for line plotting and confidence intervals
                xpred_l = linspace(x.min(), x.max(), 100)
                xm = x.mean()
                xs = x.std()
                min_val = min([x.min(), x1.min()])
                max_val = max([x.max(), x1.max()])
                x1_smooth = np.linspace(min_val, max_val, 100)

                #Append # of datapoints, species names  and rate
                best_fits_studies.append(str(studies_handles).replace("'", "").strip("[]"))
                best_fits_n.append(len(x))
                best_fits_groups.append(rate_group['Group'].iloc[0])
                best_fits_species.append(species)
                
                match = re.search(r'\[.*?\]', rate)
                if match:
                    best_fits_rate.append(rate)
    
                #Select model
                #Choose model based on p value
                v = np.argmin([lp,qp,xp])
                smo = v
                # (Or uncomment to do lowest overall RMSE, p and R^2)
                # r = np.argmin([lrmse,qrmse,xrmse])
                # s = np.argmax([lrsqu, qrsqu, xqsqu])
                # smo = mode([r,v,s])

                # Linear = lowest p-value
                if smo==0:
                    best_fits.append('linear')
                    best_fits_p.append(lp) 
                    best_fits_r2.append(np.round(lrsqu, 2))
                    best_fits_rmse.append(np.round(lrmse, 4))

                    ypred_lin_l = fl(lparams, xpred_l)
                    plt.plot(xpred_l, ypred_lin_l, color='xkcd:golden rod', 
                             label='Experimental data species response')
                    ypred_smooth_l = fl(lparams, x1_smooth)
                    plt.plot(x1_smooth, ypred_smooth_l, 
                             label='Predicted response to TA addition', 
                             color='xkcd:golden rod', linestyle='--')
                    y_NaOH = fl(lparams, x_NaOH)
                    y_Na2CO3  = fl(lparams, x_Na2CO3)

                    # Confidence interval
                    x1_smooth_const = sm.add_constant(x1_smooth)
                    conf_int = model_lin.get_prediction(x1_smooth_const).conf_int(obs=True, alpha=0.1)
                    upper_lin = conf_int[:, 1]
                    lower_lin = conf_int[:, 0]
                    plt.fill_between(x1_smooth, upper_lin, lower_lin, 
                                     color='xkcd:golden rod', alpha=0.2, 
                                     label ='90% Prediction Interval', edgecolor='None')

                    # Pre-industrial rate
                    pre_industrial_rate = fl(lparams,1.16)
                    best_fits_pre_industrial.append(pre_industrial_rate)

                    # Current rate
                    current_rate = fl(lparams,df2_['Baseline'].values[0])
                    best_fits_current.append(current_rate)

                # Quadratic = lowest p-value                   
                elif smo==1:
                    best_fits.append('parabolic')
                    best_fits_p.append(qp)
                    best_fits_r2.append(np.round(qrsqu, 2))
                    best_fits_rmse.append(np.round(qrmse, 4))
                    
                    ypred_quad_l = fq(qparams, xpred_l)
                    plt.plot(xpred_l, ypred_quad_l, color='xkcd:golden rod', 
                             label='Experimental data species response')
                    ypred_smooth_q = fq(qparams, x1_smooth)
                    plt.plot(x1_smooth, ypred_smooth_q, 
                             label='Predicted response to TA addition', 
                             color='xkcd:golden rod', linestyle='--')
                    y_NaOH = fq(qparams, x_NaOH)
                    y_Na2CO3  = fq(qparams, x_Na2CO3)
                    
                    #Inflection point
                    # Find the index of the maximum value in ypred_quad_l
                    max_index = np.argmax(ypred_quad_l)
                    
                    # Find the corresponding value in xpred_l
                    x_at_max_y = xpred_l[max_index]

                    #Append inflection point info
                    inflection_species.append(species)
                    inflection_rate.append(rate)
                    inflection_point.append(np.round(x_at_max_y, 2))

                    xp1_smooth = polynomial_features.transform(x1_smooth.reshape(-1, 1))
                    conf_int = model_quad.get_prediction(xp1_smooth).conf_int(obs=True, alpha=0.1)
                    upper_lin = conf_int[:, 1]
                    lower_lin = conf_int[:, 0]
                    plt.fill_between(x1_smooth, upper_lin, lower_lin, color='xkcd:golden rod', alpha=0.2, label ='90% Prediction Interval', edgecolor='None')

                    # Pre-industrial rate
                    pre_industrial_rate = fq(qparams,1.16)
                    best_fits_pre_industrial.append(pre_industrial_rate)

                    # Current rate
                    current_rate = fq(qparams,df2_['Baseline'].values[0])
                    best_fits_current.append(current_rate)

                # Exponential = lowest p-value
                elif smo==2:
                    best_fits.append('exponential')
                    best_fits_p.append(xp)
                    best_fits_r2.append(np.round(xqsqu, 2))
                    best_fits_rmse.append(np.round(xrmse, 4))
                    
                    plt.plot(xpred_l, ypred_exp_l, color='xkcd:golden rod', 
                             label='Experimental data species response')
                    ypred_smooth_l=fe(xparam, x1_smooth)
                    plt.plot(x1_smooth, ypred_smooth_l, 
                             label='Predicted response to TA addition', 
                             color='xkcd:golden rod', linestyle='--')
                    y_NaOH = fe(xparam, x_NaOH)
                    y_Na2CO3  = fe(xparam, x_Na2CO3)

                    ypred_l, ciLo, ciUp, p, rsqr = prediction_interval(f,xparam,x1_smooth,x.values,y.values,.1,True)
                    plt.fill_between(x1_smooth, ciLo, ciUp, color="xkcd:golden rod", 
                                     label=r'90% prediction interval', 
                                     alpha=0.2, edgecolor='None')

                    # Pre-industrial rate
                    pre_industrial_rate = fe(xparam,1.16)
                    best_fits_pre_industrial.append(pre_industrial_rate)

                    # Current rate
                    current_rate = fe(xparam,df2_['Baseline'].values[0])
                    best_fits_current.append(current_rate)
          

                # Half current rate
                yhalf = current_rate - current_rate*0.5
                a,b=x_NaOH.shape
                if b>0:
                    val1 = np.min(np.abs(y_NaOH-yhalf))
                    # if val1<=yeps:
                    ind1 =  np.argmin(np.abs(y_NaOH-yhalf))
                    NaOH_threshold_species.append(species)
                    NaOH_threshold_groups.append(rate_group['Group'].iloc[0])
                    match = re.search(r'\[.*?\]', rate)
                    if match:
                        NaOH_threshold_rate.append(rate)
                    NaOH_thresholds.append(np.round(x_NaOH[ind1],3))
                    NaOH_threshold_amount.append(labels[ind1])

                a,b=x_Na2CO3.shape
                if b>0:
                    ind2 =  np.argmin(np.abs(y_Na2CO3-yhalf))
                    Na2CO3_threshold_species.append(species)
                    Na2CO3_threshold_groups.append(rate_group['Group'].iloc[0])
                    match = re.search(r'\[.*?\]', rate)
                    if match:
                        Na2CO3_threshold_rate.append(rate)
                    Na2CO3_thresholds.append(np.round(x_Na2CO3[ind2],3))
                    Na2CO3_threshold_amount.append(labels[ind2])
 
                for i in range(len(x_NaOH)):
                    plt.scatter(x_NaOH[i], y_NaOH[i],
                                    color=cmap_NaOH[i],
                                    label = f'{labels[i]} µmol/kg NaOH addition',
                                    marker='*', edgecolor='grey', linewidth=0.5, 
                                    s=75, zorder=10, alpha=0.8)

                for i in range(len(x_Na2CO3)):
                    plt.scatter(x_Na2CO3[i], y_Na2CO3[i],
                                color=cmap_Na2CO3[i],
                                label = f'{labels[i]} µmol/kg Na₂CO₃ addition',
                                marker='d', edgecolor='grey', linewidth=0.5, 
                                s=30, zorder=10, alpha=0.8)

                # Plot calc rate pre-industrial
                plt.axhline(y=pre_industrial_rate, color='xkcd:grey', 
                            label='Pre-industrial calc rate', alpha=0.4, 
                            zorder=0, linestyle = '-.')
 
                # Plot calc rate current
                plt.axhline(y=current_rate, color='xkcd:blue', 
                            label='Current calc rate', alpha=0.4, 
                            zorder=0, linestyle = '--')

            # If none of the regressions work = Neutal response               
            else:
                studies_handles = []
                for i, study in enumerate(rate_group['Study'][~nan_mask].unique()):
                    mask = rate_group['Study'] == study
                    studies = plt.scatter(x[mask], y[mask], label=study, 
                                          color=colors[i], alpha=0.7, edgecolor='none')
                    studies_handles.append(study)
                best_fits_studies.append(str(studies_handles).replace("'", "").strip("[]"))
                best_fits_n.append(len(x))
                best_fits_groups.append(rate_group['Group'].iloc[0])
                best_fits_species.append(species)
                match = re.search(r'\[.*?\]', rate)
                if match:
                    best_fits_rate.append(rate)

                best_fits.append('neutral')
                best_fits_p.append(np.nan)
                best_fits_r2.append(np.nan)
                best_fits_rmse.append(np.nan)
                best_fits_pre_industrial.append(np.nan)
                best_fits_current.append(np.nan)
                best_fits_pH_9.append(np.nan)
                
            # #Add line at baseline TA:DIC
            plt.axvline(x=df2_['Baseline'].values[0], color='xkcd:grey', 
                        label='Baseline TA:DIC', alpha=0.4, zorder=0, 
                        linestyle = ':')
            
            #Add line at y=0
            plt.axhline(y=0, color='xkcd:red', label='Calc rate = 0', 
                        alpha=0.5, zorder=0, linestyle='-')

            #Formatting
            species_formatted = species.replace(" ", "\u00A0")

            #Axes labels        
            plt.xlabel('TA:DIC', fontsize='medium')
            
            if rate == rate1:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[{mmol}$ ${g}^{-1}$ ${hr}^{-1}$]', fontsize='medium')
                rate_label = 'Rate 1'
            if rate == rate2:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[{mmol}$ ${m}^{-2}$ ${hr}^{-1}$]', fontsize='medium')
                rate_label = 'Rate 2'
            if rate == rate3:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[{mmol}$ ${m}^{-3}$ ${hr}^{-1}$]', fontsize='medium')
                rate_label = 'Rate 3'
            if rate == rate4:
                plt.ylabel('Calcification rate $CaCO_3$ $[{mmol}$ $\#^{-1}$ $hr^{-1}]$', fontsize='medium')
                rate_label = 'Rate 4'
            if rate == rate5:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[{mmol}$ ${hr}^{-1}$]', fontsize='medium')
                rate_label = 'Rate 5'
            if rate == rate6:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[\\mathrm{mmol}\\,\\#^{-1}]$', fontsize='medium')
                rate_label = 'Rate 6'
            if rate == rate7:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[{mmol}$ ${cm}^{-2}$]', fontsize='medium')
                rate_label = 'Rate 7'
            if rate == rate8:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[{\%}$]', fontsize='medium')
                rate_label = 'Rate 8'
            if rate == rate9:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[{mg}$]', fontsize='medium')
                rate_label = 'Rate 9'  
            if rate == rate10:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[{\%}$ ${hr}^{-1}$]', fontsize='medium')
                rate_label = 'Rate 10'
            if rate == rate11:
                plt.ylabel('Calcification rate ${CaCO}_3$ $[{\#}$ ${hr}^{-1}$]', fontsize='medium')
                rate_label = 'Rate 11'
            if rate == rate12:
                plt.ylabel('Calcification rate ${CaCO}_3$ ${hr}^{-1}$]', fontsize='medium')
                rate_label = 'Rate 12'
            if rate == rate13:
                plt.ylabel('Calcification rate (shell length) ${CaCO}_3$ $[{mm}$]', fontsize='medium')
                rate_label = 'Rate 13'
            if rate == rate14:
                plt.ylabel('Calcification rate (shell thickness) ${CaCO}_3$ $[{µm}$]', fontsize='medium')
                rate_label = 'Rate 14'

            #Formatting
            species_formatted = species.replace(" ", "\u00A0")
            plt.title('{} (${}$)'.format(rate_group['Group'].iloc[0], 
                                         species_formatted),fontsize=8)

            #Complete legend
            plt.legend(loc='center left', facecolor='white', framealpha=0.7, 
                       bbox_to_anchor=(1.05, 0.5), fontsize=6, ncol=1)

            plt.tight_layout()
            plt.savefig(filepath+'/your/path/{} ({}) {}'.format(
                rate_group['Group'].iloc[0], species_formatted, rate_label)+' TA-DIC.png', bbox_inches='tight')
            plt.show()

#Make df with responses and stats
best_fits_df = pd.DataFrame({
    'Studies': best_fits_studies,
    '# datapoints': best_fits_n,
    'Group': best_fits_groups,
    'Species': best_fits_species,
    'Rate unit': best_fits_rate,
    'Response': best_fits,
    'p-value': best_fits_p,
    'R2': best_fits_r2,
    'RMSE': best_fits_rmse,
    'Pre-industrial rate': best_fits_pre_industrial,
    'Current rate': best_fits_current,
    'Rate at pH 9': best_fits_pH_9,
    })

# Make df with inflection points
inflection_points = pd.DataFrame({
    'Studies': inflection_species,
    'Rate unit': inflection_rate,
    'Inflection point': inflection_point})

# Make df with thresholds
NaOH_threshold_df = pd.DataFrame({
    'Studies': NaOH_threshold_species,
    'Group': NaOH_threshold_groups,
    'Rate unit': NaOH_threshold_rate,
    'Threshold': NaOH_thresholds,
    'TA addition': NaOH_threshold_amount
    })

Na2CO3_threshold_df = pd.DataFrame({
    'Studies': Na2CO3_threshold_species,
    'Group': Na2CO3_threshold_groups,
    'Rate unit': Na2CO3_threshold_rate,
    'Threshold': Na2CO3_thresholds,
    'TA addition': Na2CO3_threshold_amount})