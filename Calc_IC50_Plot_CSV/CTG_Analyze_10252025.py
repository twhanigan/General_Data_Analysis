#Analyze CTG output
import os
import sys
import glob
import fnmatch
from itertools import cycle
import numpy as np
import pandas as pd
import math 
from pandas import ExcelWriter
from pandas import ExcelFile
import scipy
from scipy import stats
from scipy.stats import pearsonr
import warnings
from scipy.stats import ttest_ind
from sklearn.metrics.pairwise import pairwise_distances
from more_itertools import unique_everseen
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit.models import LinearModel, StepModel
from lmfit import Model
import lmfit
plt.rcParams['font.family']='Arial'
plt.rc('axes', linewidth=2)

my_path = "E:\\Clariostar\\B819_Combined\\10252025"

def get_filepaths(directory):
    """Get full filepaths of all files in a directory, including sub-directories"""
    file_paths = []
    to_scan = '.CSV'
    for root, directories, files in os.walk(directory,topdown=True):
        directories[:] = [d for d in directories]
        for filename in files:
            filepath = os.path.join(root, filename)
            if filepath.endswith(".CSV"):
                file_paths.append(filepath)
    return file_paths

files = get_filepaths(my_path)
print(files)
Compound_Names = ['Leflun_10','Leflun_5','B508']

def generate_concentrations(max_conc, number_dilutions, dilution_factor):
    """
    Generate a list of concentrations for serial dilutions
    max_conc : float, Starting (maximum) concentration.
    number_dilutions : int,Total number of concentrations (including the starting one).
    dilution_factor : float,The fold dilution between each step (e.g. 3 for 3-fold serial dilution).
    list of float, List of concentrations.
    """
    conc = [max_conc / (dilution_factor ** i) for i in range(number_dilutions)]
    conc_name = ['Negative Ctrl'] + conc[:-1]
    conc_list = [conc[-1]] + conc[:-1]
    log_conc = np.log10(conc_list)
    return conc_list, log_conc, conc_name

def ll4(x,b,c,d,e):
    '''This function is basically a copy of the LL.4 function from the R drc package with
     - b: hill slope
     - c: min response
     - d: max response
     - e: EC50'''
    return(c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))

def sigmoid(x,b,c,d,e):
    '''This function is basically a copy of the LL.4 functon from the R drc package with
     - b: hill slope
     - c: min response
     - d: max response
     - e: EC50'''
    return(c+(d-c)/(1+np.exp(b*(np.log(e)-np.log(x)))))

def IC50(x,a,b,c):
    # a: top
    # b: bottom
    # c: IC50
    return(b+(a-b)/1+(10**(x-np.log(c))))  

def pDose(x):
    return(np.log10(x))

def normalize(frame, frame_name):
    # Group by 'compound' and loop over each group (compound)
    compoundData = frame.groupby(['compound'], sort=False)
    
    for name, group in compoundData:
        # Fit the sigmoid model to the group (just an example, replace with your actual fitting logic)
        gmodel = Model(sigmoid)
        
        # Calculate the min and max Response values for the current compound
        val_max = group['Response'].loc[group['Concentration']==sorted(set(group['Concentration']))[0]].mean()
        val_min = negative_control

        # Print debugging info (you can remove these in production)
        print(f"Processing compound: {name}")
        print("val_min:", val_min, "val_max:", val_max)

        # Calculate normalized values for this group (avoid using a loop, vectorize it)
        group['Normalized'] = (group['Response'] - val_min) / (val_max - val_min)  # Example normalization

        # Update the main DataFrame with the normalized values (important to use loc to ensure correct indexing)
        frame.loc[group.index, 'Normalized'] = group['Normalized']

    # After processing all compounds, save the final DataFrame to Excel
    final = frame.copy()
    final.to_excel(frame_name + 'Normalized.xlsx', index=False)  # Save to Excel without index

    return final

def _initial_guesses_from_group(group):
    """Robust initial guesses for IC50 and Hill slope b."""
    g = group.copy()

    # Clean & sort
    g = g[['Concentration', 'Normalized']].dropna().sort_values('Concentration')
    if g.empty or g['Concentration'].nunique() < 3:
        # fallback in pathological cases
        xmid = np.median(g['Concentration']) if not g.empty else 1e-7
        return xmid, 1.0, 1.0, 0.0  # ic50, b, d(top), c(bottom)

    # Estimate top (d) at lowest conc; bottom (c) at highest conc
    # (use mean in case of replicates at same conc)
    d_top = g.loc[g['Concentration'] == g['Concentration'].min(), 'Normalized'].mean()
    c_bot = g.loc[g['Concentration'] == g['Concentration'].max(), 'Normalized'].mean()

    # If the curve is inverted (increasing with dose), swap roles
    # Decide by comparing overall trend
    trend = np.corrcoef(np.log10(g['Concentration'].values), g['Normalized'].values)[0,1]
    if trend > 0:  # increasing response with dose: bottom at min, top at max
        d_top, c_bot = (
            g.loc[g['Concentration'] == g['Concentration'].max(), 'Normalized'].mean(),
            g.loc[g['Concentration'] == g['Concentration'].min(), 'Normalized'].mean()
        )

    # Midpoint target response
    y_mid = (d_top + c_bot) / 2.0

    # Choose IC50 as conc whose response is closest to midpoint
    idx_mid = (g['Normalized'] - y_mid).abs().idxmin()
    init_ic50 = float(g.loc[idx_mid, 'Concentration'])

    # Estimate local slope dy/dlog10(x) around midpoint using a small window
    # Grab up to 5 points centered at idx_mid for a robust linear fit
    g = g.reset_index(drop=True)
    mid_pos = int(np.where(g.index == g.index[g['Concentration'] == init_ic50][0])[0][0]) \
              if (g['Concentration'] == init_ic50).any() else g['Concentration'].sub(init_ic50).abs().idxmin()
    lo = max(0, mid_pos - 2)
    hi = min(len(g) - 1, mid_pos + 2)
    window = g.iloc[lo:hi+1]

    X = np.log10(window['Concentration'].values)
    Y = window['Normalized'].values
    if len(window) >= 2 and np.isfinite(X).all():
        # simple linear fit Y = m*logC + b0
        m = np.polyfit(X, Y, 1)[0]  # dy/dlog10(x)
    else:
        m = 0.0  # fallback

    # Convert dy/dlog10(x) to Hill slope b using derivative at IC50 for a 4PL
    denom = (d_top - c_bot)
    if denom == 0 or not np.isfinite(denom):
        init_b = 1.0
    else:
        init_b = -4.0 * (m / denom) / np.log(10)  # b ≈ -4/(ln10) * (dy/dlog10(x))/(d-c)

    # Reasonable bounds/clip to avoid crazy starts
    if not np.isfinite(init_b):
        init_b = 1.0
    init_b = float(np.clip(init_b, -3.0, 3.0))
    if abs(init_b) < 1e-3:
        init_b = 0.5 if trend < 0 else -0.5  # keep direction

    return init_ic50, init_b, float(d_top), float(c_bot)

def get_curve_lm_b(frame):
    compoundData = frame.groupby(['compound','frame_name'], sort=False)
    fitData = []
    evaluated_frame = pd.DataFrame(
        columns=['dose','new','dely','upperbound','lowerbound','compound','group_num','frame_name'],
        dtype=float
    )
    nums = [0,1,2]
    licycle = cycle(nums)
    fig, ax = plt.subplots(figsize=(6, 4.5))

    for name, group in compoundData:
        gmodel = Model(sigmoid)

        # initial guesses
        init_ic50, init_slope, val_max, val_min = _initial_guesses_from_group(group)
        # decide parameter bounds based on direction
        decreasing = (val_max > val_min)  # typical inhibition curve
        params = gmodel.make_params()

        if decreasing:
            # typical: top ~1, bottom ~0-1
            params.add('d', min=0.7, max=1.5)   # top
            params.add('c', min=-0.1, max=0.7) # bottom
            params.add('e', min=frame['Concentration'].min()*1e-3, max=frame['Concentration'].max()*1e3)
            params.add('b', min=-3, max=3)
        else:
            # increasing curve
            params.add('d', min=-0.1, max=2.0)
            params.add('c', min=-0.1, max=2.0)
            params.add('e', min=frame['Concentration'].min()*1e-3, max=frame['Concentration'].max()*1e3)
            params.add('b', min=-3, max=3)

        f = next(licycle)
        xdata = group['Concentration'].values
        ydata = group['Normalized'].values
        refDose = np.linspace(np.min(xdata), np.max(xdata), 20000)
        dose = pDose(refDose)  # your plotting/transform function

        try:
            result = gmodel.fit(
                ydata, params, x=xdata,
                b=init_slope, c=val_min, d=val_max, e=init_ic50
            )
        except (ValueError, RuntimeError):
            # gentle fallback if the first try fails
            result = gmodel.fit(
                ydata, params, x=xdata,
                b=np.sign(init_slope) if init_slope != 0 else 1.0,
                c=np.clip(val_min, 0, 1.1), d=np.clip(val_max, 0, 1.5),
                e=np.clip(init_ic50, np.min(xdata), np.max(xdata))
            )

        new = result.eval(x=refDose)
        dely = 0.3 * max(1e-6, val_min)  # uncertainty ribbon
        upper_bound = new + dely
        lower_bound = new - dely

        evaluated_frame = pd.concat([
            evaluated_frame,
            pd.DataFrame({
                'dose': dose, 'new': new, 'dely': dely,
                'upperbound': upper_bound, 'lowerbound': lower_bound,
                'compound': name[0], 'group_num': f, 'frame_name': name[1]
            })
        ], ignore_index=True)

    return (frame, evaluated_frame)

Conc, LogConc, Concentrations = generate_concentrations(max_conc=1e-5, number_dilutions=12, dilution_factor=3)
print(Concentrations)
# Prepare an empty list to store DataFrames instead of using append (which is deprecated)
frame_list = []

# Loop through files and process each
for file in files:
    frame = pd.read_csv(file, skiprows=9, index_col=0, header=0)
    frame_name = str(os.path.basename(file))
    frame_cpd_1 = frame.loc['B':'C', :]
    frame_cpd_2 = frame.loc['D':'E', :]
    frame_cpd_3 = frame.loc['F':'G', :]
    negative_control = frame.loc['H', :].mean()
    frames = [frame_cpd_1, frame_cpd_2, frame_cpd_3]
    print(frame_name)
    lst = []
    d = -1

    if frame_name == "TRno1.CSV":
        for frame in frames:
            d += 1
            c = -1
            name = Compound_Names[d]
            for col in frame:
                c += 1
                average = frame[col].mean()
                error = np.std(frame[col])
                label = Concentrations[c]
                if name != 'Paclitaxel':
                    x_value = Conc[c]
                    log_x_value = LogConc[c]
                else:
                    x_value = alt_Conc[c]
                    log_x_value = LogConc_Alt[c]
                lst.append([frame_name, name, x_value, log_x_value, average, error])
    else:
        for frame in frames:
            d += 1
            c = -1
            name = Compound_Names[d]
            for col in frame:
                c += 1
                average = frame[col].mean()
                error = np.std(frame[col])
                label = Concentrations[c]
                if name != 'Paclitaxel':
                    x_value = Conc[c]
                    log_x_value = LogConc[c]
                else:
                    x_value = alt_Conc_b[c]
                    log_x_value = LogConc_Alt[c]
                lst.append([frame_name, name, x_value, log_x_value, average, error])

    # Create and normalize the DataFrame
    final = pd.DataFrame(lst, columns=['frame_name', 'compound', 'Concentration', 'Log Concentration', 'Response', 'Error'])
    normal = normalize(final, frame_name)

    # Append to list (no deprecated append used)
    frame_list.append(normal)

# Combine all normalized frames
total_frame = pd.concat(frame_list, ignore_index=True)

# Fit curve and export results
a, b = get_curve_lm_b(total_frame)
a.to_csv('Combined_Frame_10252025.csv')
b.to_csv('Combined_Evaluated_Frame_10252025.csv')
