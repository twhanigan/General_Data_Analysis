# Plot Analyzed CTG
from itertools import cycle
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'
plt.rc('axes', linewidth=2)

# --- Load combined outputs from the first program ---
frame = pd.read_csv('Combined_Frame_10252025.csv')
evaluated = pd.read_csv('Combined_Evaluated_Frame_10252025.csv')

# Ensure numeric dtypes (robust against CSV typing issues)
for tbl, cols in [
    (frame, ['Log Concentration', 'Normalized']),
    (evaluated, ['dose', 'new', 'lowerbound', 'upperbound'])
]:
    for c in cols:
        tbl[c] = pd.to_numeric(tbl[c], errors='coerce')

# Color palette: stable color per compound across plates
compounds = frame['compound'].dropna().unique().tolist()
palette = sns.color_palette(["steelblue", "firebrick", "red", "k", "grey"], n_colors=max(5, len(compounds)))
comp2color = {cmpd: palette[i % len(palette)] for i, cmpd in enumerate(compounds)}

# Group data
compoundData = frame.groupby(['compound', 'frame_name'], sort=False)

fig, ax = plt.subplots(figsize=(4, 3.5))

# Track legend entries (one per compound)
seen = set()

for (cmpd, fname), group in compoundData:
    color = comp2color.get(cmpd, 'grey')

    # Slice evaluated rows for this compound/file
    mask_eval = (evaluated['compound'] == cmpd) & (evaluated['frame_name'] == fname)
    x_fit = evaluated.loc[mask_eval, 'dose'].values
    y_fit = evaluated.loc[mask_eval, 'new'].values
    y_lo  = evaluated.loc[mask_eval, 'lowerbound'].values
    y_hi  = evaluated.loc[mask_eval, 'upperbound'].values

    # Slice observed points
    mask_obs = (frame['compound'] == cmpd) & (frame['frame_name'] == fname)
    x_obs = frame.loc[mask_obs, 'Log Concentration'].values
    y_obs = frame.loc[mask_obs, 'Normalized'].values

    label = cmpd if cmpd not in seen else None
    if label is not None:
        seen.add(cmpd)

    ax.plot(x_fit, y_fit, linewidth=2, alpha=0.9, color=color, label=label)
    ax.scatter(x_obs, y_obs, s=18, alpha=0.85, color=color, edgecolors='none')
    ax.fill_between(x_fit, y_lo, y_hi, alpha=0.08, color=color)

# Dynamic axes based on actual evaluated dose range (from fitter)
if not evaluated['dose'].dropna().empty:
    x_min = np.nanmin(evaluated['dose'].values)
    x_max = np.nanmax(evaluated['dose'].values)
    # Keep left-to-right from low conc (more negative) to high conc (less negative), like your original
    ax.set_xlim(x_min, x_max)

# Y-limits with a bit of padding
all_y = np.concatenate([
    frame['Normalized'].dropna().values,
    evaluated['new'].dropna().values,
    evaluated['lowerbound'].dropna().values,
    evaluated['upperbound'].dropna().values
]) if not frame.empty and not evaluated.empty else np.array([0,1])
y_min = max(0.0, np.nanmin(all_y) - 0.05)
y_max = min(2.0, np.nanmax(all_y) + 0.05)
ax.set_ylim(y_min, max(1.0, y_max))

# Legend (compounds)
ax.legend(title='', bbox_to_anchor=(1.30, 1.05), fontsize=10, frameon=False)

# Labels & styling
plt.xlabel('Log[I], M', fontsize=18, fontname="Arial", fontweight='bold')
plt.ylabel('Relative Viability', fontsize=18, fontname="Arial", fontweight='bold')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.savefig('CTG_Output_10252025.png', dpi=600, bbox_inches='tight', pad_inches=2, transparent=True)
plt.show()