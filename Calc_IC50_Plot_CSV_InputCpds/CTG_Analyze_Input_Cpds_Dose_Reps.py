# Analyze CTG output
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

# =========================
# USER INPUTS
# =========================
NUM_COMPOUNDS       = 2          # e.g., 2 compounds
REPS_PER_COMPOUND   = 3          # e.g., 3 replicates per compound
NUM_DOSES           = 12         # e.g., 12 columns of doses (1..12)
DILUTION_FACTOR     = 3.0        # e.g., 3-fold serial dilution
MAX_CONC_LIST       = [4e-5, 2e-5]  # per-compound starting conc; if len<NUM_COMPOUNDS it will repeat
START_ROW           = 'B'        # first data row
NEG_CTRL_ROW        = 'H'        # negative control row
Compound_Names      = ['Puromycin','Blasticidin']   # Compound names

# Directory or file
my_path = "C:\\Users\\thanigan\\Documents\\GitHub\\General_Data_Analysis\\General_Data_Analysis\\CTG_2Cpd\\"

def get_filepaths(directory):
    """Get full filepaths of all files in a directory, including sub-directories.
       If a single file is passed, return [that file].
    """
    if os.path.isfile(directory):
        return [directory]
    file_paths = []
    for root, directories, files in os.walk(directory, topdown=True):
        for filename in files:
            filepath = os.path.join(root, filename)
            if filepath.lower().endswith(".csv"):
                file_paths.append(filepath)
    return file_paths

files = get_filepaths(my_path)
print(files)

# ---------- helper functions ----------
def generate_concentrations(max_conc, number_dilutions, dilution_factor):
    """Left->right decreasing concentration (col 1 = highest)."""
    conc = [max_conc / (dilution_factor ** i) for i in range(number_dilutions)]
    conc_list = conc[:]
    conc_name = [f"C{i+1}" for i in range(number_dilutions)]
    log_conc = np.log10(conc_list)
    return conc_list, log_conc, conc_name

def ll4(x,b,c,d,e):
    return (c + (d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))

def sigmoid(x,b,c,d,e):
    return (c + (d-c)/(1+np.exp(b*(np.log(e)-np.log(x)))))

def IC50(x,a,b,c):
    return (b+(a-b)/1+(10**(x-np.log(c))))

def pDose(x):
    return (np.log10(x))

def normalize(frame, frame_name):
    frame = frame.copy()
    # dtype safety
    for col in ['Concentration','Response']:
        frame[col] = pd.to_numeric(frame[col], errors='coerce')
    compoundData = frame.groupby(['compound'], sort=False)
    for name, group in compoundData:
        gmodel = Model(sigmoid)
        val_max = group.loc[group['Concentration']==group['Concentration'].min(), 'Response'].mean()
        val_min = negative_control
        frame.loc[group.index, 'Normalized'] = (frame.loc[group.index, 'Response'] - val_min) / (val_max - val_min)
    frame['Normalized'] = pd.to_numeric(frame['Normalized'], errors='coerce')
    final = frame.copy()
    final.to_excel(frame_name + 'Normalized.xlsx', index=False)
    return final

def _initial_guesses_from_group(group):
    g = group[['Concentration','Normalized']].dropna().sort_values('Concentration').copy()
    if g.empty or g['Concentration'].nunique() < 3:
        xmid = np.median(g['Concentration']) if not g.empty else 1e-7
        return xmid, 1.0, 1.0, 0.0
    d_top = g.loc[g['Concentration'].min()==g['Concentration'],'Normalized'].mean()
    c_bot = g.loc[g['Concentration'].max()==g['Concentration'],'Normalized'].mean()
    trend = np.corrcoef(np.log10(g['Concentration'].values), g['Normalized'].values)[0,1]
    if trend > 0:
        d_top, c_bot = (
            g.loc[g['Concentration'].max()==g['Concentration'],'Normalized'].mean(),
            g.loc[g['Concentration'].min()==g['Concentration'],'Normalized'].mean()
        )
    y_mid = (d_top + c_bot)/2.0
    idx_mid = (g['Normalized'] - y_mid).abs().idxmin()
    init_ic50 = float(g.loc[idx_mid, 'Concentration'])
    g = g.reset_index(drop=True)
    mid_pos = int(g['Concentration'].sub(init_ic50).abs().idxmin())
    lo, hi = max(0, mid_pos-2), min(len(g)-1, mid_pos+2)
    window = g.iloc[lo:hi+1]
    X = np.log10(window['Concentration'].values)
    Y = window['Normalized'].values
    if len(window) >= 2 and np.isfinite(X).all():
        m = np.polyfit(X, Y, 1)[0]
    else:
        m = 0.0
    denom = (d_top - c_bot)
    if denom == 0 or not np.isfinite(denom):
        init_b = 1.0
    else:
        init_b = -4.0 * (m / denom) / np.log(10)
    if not np.isfinite(init_b):
        init_b = 1.0
    init_b = float(np.clip(init_b, -3.0, 3.0))
    if abs(init_b) < 1e-3:
        init_b = 0.5 if trend < 0 else -0.5
    return float(init_ic50), float(init_b), float(d_top), float(c_bot)

def get_curve_lm_b(frame):
    # dtype safety
    frame = frame.copy()
    frame['Concentration'] = pd.to_numeric(frame['Concentration'], errors='coerce')
    frame['Normalized']    = pd.to_numeric(frame['Normalized'], errors='coerce')
    frame = frame.dropna(subset=['Concentration','Normalized'])

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
        init_ic50, init_slope, val_max, val_min = _initial_guesses_from_group(group)
        decreasing = (val_max > val_min)
        params = gmodel.make_params()
        params.add('d', min=0.7 if decreasing else -0.1, max=1.5 if decreasing else 2.0)
        params.add('c', min=-0.1, max=0.7 if decreasing else 2.0)
        params.add('e', min=frame['Concentration'].min()*1e-3, max=frame['Concentration'].max()*1e3)
        params.add('b', min=-3, max=3)

        f = next(licycle)
        xdata = group['Concentration'].astype(float).values
        ydata = group['Normalized'].astype(float).values
        refDose = np.linspace(np.nanmin(xdata), np.nanmax(xdata), 20000)
        dose = pDose(refDose)

        try:
            result = gmodel.fit(ydata, params, x=xdata, b=init_slope, c=val_min, d=val_max, e=init_ic50)
        except (ValueError, RuntimeError):
            result = gmodel.fit(
                ydata, params, x=xdata,
                b=np.sign(init_slope) if init_slope != 0 else 1.0,
                c=np.clip(val_min, 0, 1.1), d=np.clip(val_max, 0, 1.5),
                e=np.clip(init_ic50, np.nanmin(xdata), np.nanmax(xdata))
            )

        new = result.eval(x=refDose)
        dely = 0.3 * max(1e-6, val_min)
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

# =========================
# generalized plate mapping
# =========================
row_labels = list("ABCDEFGH")
start_idx = row_labels.index(START_ROW)
ctrl_idx  = row_labels.index(NEG_CTRL_ROW)

# sanity: make enough names and max concs
if len(Compound_Names) < NUM_COMPOUNDS:
    Compound_Names = Compound_Names + [f"Compound_{i+1}" for i in range(len(Compound_Names), NUM_COMPOUNDS)]
if len(MAX_CONC_LIST) < NUM_COMPOUNDS:
    MAX_CONC_LIST = MAX_CONC_LIST + [MAX_CONC_LIST[0]]*(NUM_COMPOUNDS - len(MAX_CONC_LIST))

# precompute concentration lists per compound
conc_info = []
for i in range(NUM_COMPOUNDS):
    conc_i, logc_i, names_i = generate_concentrations(MAX_CONC_LIST[i], NUM_DOSES, DILUTION_FACTOR)
    conc_info.append((conc_i, logc_i, names_i))

# ================
# main file loop
# ================
frame_list = []

for file in files:
    # CLARIOstar export: header row begins at line 9 (0-based index 8)
    frame = pd.read_csv(file, skiprows=8, index_col=0, header=0)
    # keep only the first NUM_DOSES columns (1..NUM_DOSES)
    frame = frame[[str(i) for i in range(1, NUM_DOSES+1)]]
    frame = frame.apply(pd.to_numeric, errors='coerce')
    frame_name = str(os.path.basename(file))

    # per-plate negative control
    negative_control = frame.loc[NEG_CTRL_ROW, :].mean()

    # build blocks for each compound
    lst = []
    for comp_idx in range(NUM_COMPOUNDS):
        comp_name = Compound_Names[comp_idx]
        rows_for_comp = row_labels[start_idx + comp_idx*REPS_PER_COMPOUND : start_idx + (comp_idx+1)*REPS_PER_COMPOUND]
        # guard against overflow
        rows_for_comp = [r for r in rows_for_comp if r in row_labels and row_labels.index(r) != ctrl_idx]
        block = frame.loc[rows_for_comp].astype(float)
        avg = block.mean(axis=0)
        err = block.std(axis=0)

        conc_list, log_conc, _ = conc_info[comp_idx]
        for c_idx, col in enumerate(block.columns):
            x_value = conc_list[c_idx]
            log_x   = log_conc[c_idx]
            average = float(avg[col])
            error   = float(err[col])
            lst.append([frame_name, comp_name, x_value, log_x, average, error])

    final = pd.DataFrame(lst, columns=['frame_name', 'compound', 'Concentration', 'Log Concentration', 'Response', 'Error'])
    # dtype safety
    for col in ['Concentration','Log Concentration','Response','Error']:
        final[col] = pd.to_numeric(final[col], errors='coerce')

    normal = normalize(final, frame_name)
    frame_list.append(normal)

# Combine all normalized frames
total_frame = pd.concat(frame_list, ignore_index=True)
a, b = get_curve_lm_b(total_frame)
a.to_csv('Combined_Frame_10252025.csv', index=False)
b.to_csv('Combined_Evaluated_Frame_10252025.csv', index=False)
