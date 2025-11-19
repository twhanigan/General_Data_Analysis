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
plt.rcParams['font.family'] = 'Arial'
plt.rc('axes', linewidth=2)

# =========================
# USER INPUTS
# =========================
NUM_COMPOUNDS       = 2          # e.g., 2 compounds
REPS_PER_COMPOUND   = 3          # e.g., 3 replicates per compound
NUM_DOSES           = 11         # total columns including the no-drug column
DILUTION_FACTOR     = 2.0        # e.g., 2-fold serial dilution
MAX_CONC_LIST       = [4e-5, 2e-5]  # per-compound starting conc; if len<NUM_COMPOUNDS it will repeat
START_ROW           = 'B'        # first data row
NEG_CTRL_ROW        = 'H'        # negative control row
Compound_Names      = ['Puromycin', 'Blasticidin']   # Compound names

# Directory or file
my_path = "C:\\Users\\thanigan\\Downloads\\Win_CTG_Data"

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
print("Files found:", files)

# ---------- helper functions ----------
def generate_concentrations(max_conc, number_dilutions, dilution_factor):
    """ 
    Concentration list needs to be [0, max, max/DF, max/DF^2, ..., lowest]
    """

    # number of real (non-zero) concentrations
    num_nonzero = number_dilutions - 1

    # highest → lowest
    conc_desc = [max_conc / (dilution_factor ** i) for i in range(num_nonzero)]

    # prepend 0 for no drug
    conc_list = [0.0] + conc_desc

    # log10 values, using NaN for 0
    log_conc = [np.nan] + [np.log10(c) for c in conc_desc]

    conc_name = [f"C{i+1}" for i in range(number_dilutions)]
    return conc_list, np.array(log_conc, dtype=float), conc_name


def ll4(x, b, c, d, e):
    return (c + (d - c) / (1 + np.exp(b * (np.log(x) - np.log(e)))))

def sigmoid(x, b, c, d, e):
    return (c + (d - c) / (1 + np.exp(b * (np.log(e) - np.log(x)))))

def IC50(x, a, b, c):
    return (b + (a - b) / (1 + (10 ** (x - np.log(c)))))

def pDose(x):
    return (np.log10(x))

def normalize(frame, frame_name):
    frame = frame.copy()
    # dtype safety
    for col in ['Concentration', 'Response']:
        frame[col] = pd.to_numeric(frame[col], errors='coerce')
    compoundData = frame.groupby(['compound'], sort=False)
    for name, group in compoundData:
        gmodel = Model(sigmoid)
        val_max = group.loc[group['Concentration'] == group['Concentration'].min(), 'Response'].mean()
        val_min = negative_control
        frame.loc[group.index, 'Normalized'] = (frame.loc[group.index, 'Response'] - val_min) / (val_max - val_min)
    frame['Normalized'] = pd.to_numeric(frame['Normalized'], errors='coerce')
    final = frame.copy()
    final.to_excel(frame_name + 'Normalized.xlsx', index=False)
    return final

def _initial_guesses_from_group(group):
    """Heuristic initial guesses for the 4-param sigmoid, robust to a 0 (no-drug) column."""
    g = group[['Concentration', 'Normalized']].dropna().copy()
    # ensure proper ordering by concentration
    g = g.sort_values('Concentration')

    if g.empty or g['Concentration'].nunique() < 3:
        xmid = np.median(g['Concentration']) if not g.empty else 1e-7
        return float(xmid), 1.0, 1.0, 0.0

    # Work on a copy where non-positive concentrations are replaced
    g_pos = g.copy()
    # smallest positive concentration, fallback to 1e-7
    if (g_pos['Concentration'] > 0).any():
        min_pos = g_pos.loc[g_pos['Concentration'] > 0, 'Concentration'].min()
    else:
        min_pos = 1e-7
    g_pos.loc[g_pos['Concentration'] <= 0, 'Concentration'] = min_pos / 10.0

    # top = high response, bottom = low response
    d_top = g['Normalized'].loc[g['Concentration'] == g['Concentration'].min()].mean()
    c_bot = g['Normalized'].loc[g['Concentration'] == g['Concentration'].max()].mean()

    # estimate overall trend using only positive concentrations (after replacement)
    x_for_trend = np.log10(g_pos['Concentration'].values)
    y_for_trend = g_pos['Normalized'].values
    try:
        trend = np.corrcoef(x_for_trend, y_for_trend)[0, 1]
    except Exception:
        trend = -1.0  # default to decreasing

    # if the curve is increasing, swap top/bottom
    if np.isfinite(trend) and trend > 0:
        d_top, c_bot = (
            g['Normalized'].loc[g['Concentration'] == g['Concentration'].max()].mean(),
            g['Normalized'].loc[g['Concentration'] == g['Concentration'].min()].mean()
        )

    # pick an x around the mid-response as IC50 seed
    y_mid = (d_top + c_bot) / 2.0
    idx_mid = (g_pos['Normalized'] - y_mid).abs().idxmin()
    init_ic50 = float(g_pos.loc[idx_mid, 'Concentration'])

    # local window around IC50 for slope estimate
    g_pos = g_pos.reset_index(drop=True)
    mid_pos = int(g_pos['Concentration'].sub(init_ic50).abs().idxmin())
    lo, hi = max(0, mid_pos - 2), min(len(g_pos) - 1, mid_pos + 2)
    window = g_pos.iloc[lo:hi + 1]

    if len(window) >= 2 and (window['Concentration'] > 0).all():
        X = np.log10(window['Concentration'].values)
        Y = window['Normalized'].values
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

    # keep slopes in a reasonable range
    init_b = float(np.clip(init_b, -3.0, 3.0))
    if abs(init_b) < 1e-3:
        init_b = 0.5 if (np.isfinite(trend) and trend < 0) else -0.5

    return float(init_ic50), float(init_b), float(d_top), float(c_bot)

def get_curve_lm_b(frame):
    # dtype safety
    frame = frame.copy()
    frame['Concentration'] = pd.to_numeric(frame['Concentration'], errors='coerce')
    frame['Normalized']    = pd.to_numeric(frame['Normalized'], errors='coerce')
    frame = frame.dropna(subset=['Concentration', 'Normalized'])

    compoundData = frame.groupby(['compound', 'frame_name'], sort=False)
    evaluated_frame = pd.DataFrame(
        columns=['dose', 'new', 'dely', 'upperbound', 'lowerbound',
                 'compound', 'group_num', 'frame_name'],
        dtype=float
    )
    nums = [0, 1, 2]
    licycle = cycle(nums)
    fig, ax = plt.subplots(figsize=(6, 4.5))

    for name, group in compoundData:
        gmodel = Model(sigmoid)
        init_ic50, init_slope, val_max, val_min = _initial_guesses_from_group(group)
        decreasing = (val_max > val_min)
        params = gmodel.make_params()
        params.add('d', min=0.7 if decreasing else -0.1, max=1.5 if decreasing else 2.0)
        params.add('c', min=-0.1, max=0.7 if decreasing else 2.0)
        params.add('e', min=frame['Concentration'].min() * 1e-3,
                   max=frame['Concentration'].max() * 1e3)
        params.add('b', min=-3, max=3)

        f = next(licycle)
        xdata = group['Concentration'].astype(float).values
        ydata = group['Normalized'].astype(float).values
        # only positive doses for the fit evaluation grid
        refDose_min = np.nanmin(xdata[xdata > 0]) if (xdata > 0).any() else np.nanmin(xdata[xdata >= 0] + 1e-12)
        refDose_max = np.nanmax(xdata)
        refDose = np.linspace(refDose_min, refDose_max, 20000)
        dose = pDose(refDose)

        try:
            result = gmodel.fit(
                ydata, params, x=xdata,
                b=init_slope, c=val_min, d=val_max, e=init_ic50
            )
        except (ValueError, RuntimeError):
            result = gmodel.fit(
                ydata, params, x=xdata,
                b=np.sign(init_slope) if init_slope != 0 else 1.0,
                c=np.clip(val_min, 0, 1.1),
                d=np.clip(val_max, 0, 1.5),
                e=np.clip(init_ic50, refDose_min, refDose_max)
            )

        new = result.eval(x=refDose)
        dely = 0.3 * max(1e-6, val_min)
        upper_bound = new + dely
        lower_bound = new - dely

        evaluated_frame = pd.concat([
            evaluated_frame,
            pd.DataFrame({
                'dose': dose,
                'new': new,
                'dely': dely,
                'upperbound': upper_bound,
                'lowerbound': lower_bound,
                'compound': name[0],
                'group_num': f,
                'frame_name': name[1]
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
    Compound_Names = Compound_Names + [
        f"Compound_{i+1}" for i in range(len(Compound_Names), NUM_COMPOUNDS)
    ]
if len(MAX_CONC_LIST) < NUM_COMPOUNDS:
    MAX_CONC_LIST = MAX_CONC_LIST + [MAX_CONC_LIST[0]] * (NUM_COMPOUNDS - len(MAX_CONC_LIST))

# precompute concentration lists per compound
conc_info = []
for i in range(NUM_COMPOUNDS):
    conc_i, logc_i, names_i = generate_concentrations(MAX_CONC_LIST[i], NUM_DOSES, DILUTION_FACTOR)
    conc_info.append((conc_i, logc_i, names_i))
print("Example conc list:", conc_i, logc_i, names_i)

# ================
# main file loop
# ================
frame_list = []

for file in files:
    # For the current export format the data start on the first row after the index column.
    frame = pd.read_csv(file, skiprows=0, index_col=0, header=0)

    # keep only the first NUM_DOSES columns (1..NUM_DOSES), including the no-drug column 1
    frame = frame[[str(i) for i in range(1, NUM_DOSES + 1)]]
    frame = frame.apply(pd.to_numeric, errors='coerce')
    frame_name = str(os.path.basename(file))

    # per-plate negative control (row NEG_CTRL_ROW, across all dose columns)
    negative_control = frame.loc[NEG_CTRL_ROW, :].mean()

    # build blocks for each compound
    lst = []
    for comp_idx in range(NUM_COMPOUNDS):
        comp_name = Compound_Names[comp_idx]
        rows_for_comp = row_labels[start_idx + comp_idx * REPS_PER_COMPOUND:
                                   start_idx + (comp_idx + 1) * REPS_PER_COMPOUND]
        # guard against overflow and skip the control row if it falls in this range
        rows_for_comp = [r for r in rows_for_comp
                         if r in row_labels and row_labels.index(r) != ctrl_idx]

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

    final = pd.DataFrame(
        lst,
        columns=['frame_name', 'compound',
                 'Concentration', 'Log Concentration',
                 'Response', 'Error']
    )

    # dtype safety
    for col in ['Concentration', 'Log Concentration', 'Response', 'Error']:
        final[col] = pd.to_numeric(final[col], errors='coerce')

    normal = normalize(final, frame_name)
    frame_list.append(normal)

# Combine all normalized frames
total_frame = pd.concat(frame_list, ignore_index=True)
a, b = get_curve_lm_b(total_frame)
a.to_csv('Combined_Frame_10252025.csv', index=False)
b.to_csv('Combined_Evaluated_Frame_10252025.csv', index=False)
