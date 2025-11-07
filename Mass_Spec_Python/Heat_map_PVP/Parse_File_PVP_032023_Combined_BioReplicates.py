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
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
# Matplotlib defaults for fonts/lines (Arial, slightly thicker axes)
rc = {"fontname": "Arial", 'fontsize': 24, 'fontweight': 'bold', 'lines.linewidth': 2}
plt.rcParams.update({
    "font.family": 'sans-serif',
    "font.sans-serif": "Arial",
    "font.size": 12,
    "axes.linewidth": 2,
    "axes.labelsize": 14,
    # "axes.labelweight": "bold",  # you had this commented
})
plt.rc('axes', linewidth=2)

# --- Load input spreadsheets ---------------------------------------------------
# Raw_Frame: main data matrix (likely protein/peptide intensities by sample)
Raw_Frame     = pd.read_excel('052121_PvP_S1.xlsx')

# hits_pbs_b: an external hit list (with at least a "Name" column and p-values)
hits_pbs_b  = pd.read_excel('hits_S2.xlsx')
hits_pbs_b.index = hits_pbs_b.Name  # set gene/protein names as index for join/groupby later

# --- Log2 transform selected sample columns -----------------------------------
# Select a contiguous slice of columns from DMSO_R1 through BMT174819_R2 (inclusive)
Log2_Frame = Raw_Frame.loc[:, "DMSO_R1":"BMT174819_R2"].copy()

# Apply log2 transform to those intensity-like columns
# (Copy again is redundant but harmless)
Log2_Frame = np.log2(Raw_Frame.loc[:, "DMSO_R1":"BMT174819_R2"]).copy()

# --- Derive a compact "Name" column from description --------------------------
# Ensure 'description' is string (avoid .str errors if there are NaNs)
Raw_Frame['description'] = Raw_Frame['description'].astype(str)

# Extract the first token (before the first space) as the protein/gene "Name"
# e.g., "OXA1L isoform..." -> "OXA1L"
Log2_Frame['Name'] = Raw_Frame['description'].str.split(" ", n=1).str.get(0)

# Median across samples (not used directly below, but kept as in original)
Median_Frame = Log2_Frame.loc[:, Log2_Frame.columns != 'Name'].median()

# Create a working matrix limited to the raw (pre-log2) sample columns
# (Comment suggests subtraction of median, but it's currently disabled)
Probe_Targets = Log2_Frame.loc[:, "DMSO_R1":"BMT174819_R2"]  # - Median_Frame

# Attach the Name column for indexing, then clean up
Probe_Targets['Name'] = Log2_Frame.Name
Competed = []  # unused, but kept for parity with original

# Clean the matrix: fill NaNs, force numeric, set index to Name, drop infs/NaNs/dupes
Probe_Targets.fillna(0, inplace=True)
Probe_Targets.loc[:, "DMSO_R1":"BMT174819_R2"] = Probe_Targets.loc[:, "DMSO_R1":"BMT174819_R2"].apply(pd.to_numeric)
Probe_Targets = Probe_Targets.set_index('Name')
Probe_Targets.replace(np.inf,  np.nan, inplace=True)
Probe_Targets.replace(-np.inf, np.nan, inplace=True)
Probe_Targets = Probe_Targets.dropna()
Probe_Targets = Probe_Targets.drop_duplicates()

# --- Correlation matrix across samples (quality check / optional) -------------
corrMatrix = Probe_Targets.corr()

# --- Create per-compound replicate averages -----------------------------------
# Averages across R1/R2 for each treatment compound (and DMSO)
Probe_Targets['DMSO Avg']        = Probe_Targets[["DMSO_R1",      "DMSO_R2"]].mean(axis=1)
Probe_Targets['BMT181525 Avg']   = Probe_Targets[["BMT181525_R1", 'BMT181525_R2']].mean(axis=1)
Probe_Targets['BMT182526 Avg']   = Probe_Targets[["BMT182526_R1", 'BMT182526_R2']].mean(axis=1)
Probe_Targets['BMT179888 Avg']   = Probe_Targets[["BMT179888_R1", 'BMT179888_R2']].mean(axis=1)
Probe_Targets['BMT174819 Avg']   = Probe_Targets[["BMT174819_R1", 'BMT174819_R2']].mean(axis=1)

Probe_Targets['DMSO Avg']        = Probe_Targets[["DMSO_R1","DMSO_R2"]].mean(axis=1)

# Differences vs. BMT174819 (target compound)
Probe_Targets['diff_DMSO'] = Probe_Targets['BMT174819 Avg'] - Probe_Targets['DMSO Avg']
Probe_Targets['diff_525']  = Probe_Targets['BMT174819 Avg'] - Probe_Targets['BMT181525 Avg']
Probe_Targets['diff_526']  = Probe_Targets['BMT174819 Avg'] - Probe_Targets['BMT182526 Avg']
Probe_Targets['diff_888']  = Probe_Targets['BMT174819 Avg'] - Probe_Targets['BMT179888 Avg']

# Columns to compare
cols_b819 = ['BMT174819_R1', 'BMT174819_R2']
# Control = all other compounds + DMSO (change to just DMSO if you prefer)
cols_ctrl = [
    'BMT181525_R1', 'BMT181525_R2',
    'BMT182526_R1', 'BMT182526_R2',
    'BMT179888_R1', 'BMT179888_R2',
    'DMSO_R1',      'DMSO_R2'
]

def ttest_row(row):
    # Extract values, ensure float arrays, drop NaNs
    a = row[cols_b819].astype(float).to_numpy()
    b = row[cols_ctrl].astype(float).to_numpy()
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    # Need at least 2 replicates on each side for a stable t-test
    if a.size < 2 or b.size < 2:
        return pd.Series({'BMT174819_ttest': np.nan, 'BMT174819_pval': np.nan})
    t, p = ttest_ind(a, b, equal_var=False, nan_policy='omit')  # Welch's t-test
    return pd.Series({'BMT174819_ttest': t, 'BMT174819_pval': p})

Probe_Targets[['BMT174819_ttest','BMT174819_pval']] = Probe_Targets.apply(ttest_row, axis=1)

# Sort so that large BMT174819 responders come first, and low responses in others are up front
Probe_Targets = Probe_Targets.sort_values(
    by=['BMT174819 Avg', 'BMT182526 Avg', 'BMT181525 Avg', 'BMT179888 Avg', 'DMSO Avg'],
    ascending=[False,            True,              True,              True,         True]
)

# --- Hit selection -------------------------------------------------------------
# Keep only rows with strong nominal significance vs. “others” (per the column created above)
hits_pbs   = Probe_Targets[Probe_Targets['BMT174819_pval'] < 0.005]
# Filter the external hit list similarly
hits_pbs_b = hits_pbs_b[hits_pbs_b['BMT174819_pval'] < 0.005]

# Prepare a matrix for clustermap (uses 'DMSO_R2' as the average; see below)
Probe_Targetsb = Probe_Targets.copy()
Probe_Targetsb['DMSO_R2'] = Probe_Targets['DMSO Avg']  # reuses the 'Avg' column into R2 slot for plotting

# Merge in-lab hits and external hits by index (gene/protein Name) and average duplicates
s1 = pd.concat((hits_pbs, hits_pbs_b), join="inner")
output = (
    s1.groupby(s1.index)  # group by Name index
      .mean()
      .sort_values(
          by=['BMT174819 Avg','BMT182526 Avg','BMT181525 Avg','BMT179888 Avg','DMSO Avg'],
          ascending=[False, True, True, True, True]
      )
      .drop_duplicates()
)

# Remove common structural proteins that often dominate heatmaps
output = output.loc[~output.index.str.contains("TUBB|NEF")]
output.to_excel('output.xlsx')
# --- Clustermap visualization --------------------------------------------------
# Heatmap of raw replicate columns (plus the overwritten DMSO_R2 above)
g = sns.clustermap(
    data=Probe_Targetsb.loc[:, "DMSO_R1":"BMT174819_R2"],
    cmap="Reds",
    cbar_kws={'label': r'Log$_{2}$ Reporter Ion Intensity'},
    figsize=(6, 6.5),
    cbar_pos=(1.3, 1.1, 0.04, 0.2),
    col_cluster=True,
    row_cluster=False,
    method='centroid',
    vmin=13, vmax=18,
    dendrogram_ratio=0.15
)
# (Alternative palette commented in your original)

# --- Axis tick customization for the heatmap ----------------------------------
pos     = g.ax_heatmap.get_yticks()
pos_x   = g.ax_heatmap.get_xticks()
labels  = g.ax_heatmap.get_yticklabels()

# Keep only the first two y-tick positions (compact y-labels)
newpos    = pos[0], pos[1]
newpos_x  = pos_x[0::1]     # every x tick
newpos_xm = pos_x[0::2]     # every other x tick for minor ticks
newpos_xm = newpos_xm + 0.5 # shift minor ticks to between columns

# Custom x labels (order should match the columns shown in the clustermap)
xlabel   = ['DMSO', '182526', '181525', '179888', 'BMT-819']

# Replace the first two y-labels explicitly, then keep the rest as-is
labelsnew = ['VDAC2', 'OXA1L'] + labels[2:]

# Apply y-ticks and styling
plt.setp(g.ax_heatmap.set_yticks(newpos))
g.ax_heatmap.tick_params(length=12, width=1, which='major', axis='y')

# Set major and minor x-ticks
plt.setp(g.ax_heatmap.set_xticks(newpos_x))
plt.setp(g.ax_heatmap.set_xticks(newpos_x, minor=True))
g.ax_heatmap.tick_params(length=6,  width=1, which='minor', axis='x')
plt.setp(g.ax_heatmap.set_xticks(newpos_xm))
g.ax_heatmap.tick_params(length=10, width=2, which='major', axis='x', pad=0)

# Hide y tick labels (clean look); you could set them to labelsnew if desired:
#   g.ax_heatmap.set_yticklabels(labelsnew)
plt.setp(g.ax_heatmap.set_yticklabels(''))

# Draw a border around the heatmap axes
for _, spine in g.ax_heatmap.spines.items():
    spine.set_visible(True)

# Apply x labels, rotated, bold Arial
plt.setp(
    g.ax_heatmap.set_xticklabels(xlabel),
    rotation=45, fontsize=18, fontname='Arial', fontweight='bold', ha='right'
)

# Remove y-axis label
plt.setp(g.ax_heatmap.set_ylabel(''), fontsize=18, fontweight='bold', fontname='Arial', ha='center')

# --- Save and display ----------------------------------------------------------
plt.savefig(
    'HeatMap_PVP_S1_replicates_032023_E.png',
    dpi=600, transparent=True, bbox_inches='tight', pad_inches=2
)
plt.show()



