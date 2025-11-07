import pandas as pd
from bioinfokit import visuz
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
from matplotlib import style
import numpy as np
from adjustText import adjust_text
from pylab import *

# --- Matplotlib global style --------------------------------------------------
# A few global rcParams to get consistent typography and thicker axes lines
rc = {"fontname": "Arial", 'fontsize': 24, 'fontweight': 'bold', 'lines.linewidth': 2}
plt.rcParams.update({
    "font.family": 'sans-serif',
    "font.sans-serif": "Arial",
    "font.size": 12,
    "axes.linewidth": 2,
    "axes.labelsize": 14,
    # "axes.labelweight": 'bold',  # optional, you had this commented out
})
plt.rc('axes', linewidth=2)

# --- Load & tidy the data -----------------------------------------------------
# Read the CSV with explicit dtypes for critical columns.
dfs_cat = pd.read_csv(
    'PC3_Protein_TMM_Annotated_040523.csv',
    dtype={'logFC': float, 'PValue': float, 'genomic_pos.chr': str}
)

# Remove decoy/reverse hits: guard against NaNs with na=False (FIX)
dfs_cat = dfs_cat[~dfs_cat['accession'].str.contains('Reverse', na=False)]

# De-duplicate rows to prevent double plotting
dfs_cat = dfs_cat.drop_duplicates(ignore_index=True)

# Keep only the columns needed for plotting/labeling and rename for convenience
sens = dfs_cat[['symbol', 'PValue', 'logFC', 'genomic_pos.chr', 'pathway.kegg', 'go.CC', 'summary']].copy()
sens.columns = ['Name', 'pval', 'fold_change', 'chromosome', 'kegg', 'go.CC', 'summary']

# Drop any residual duplicates after renaming
sens = sens.drop_duplicates()

# Manually inflate specific P-values to push their points down (on log scale) for visibility
# Use .loc column assignment to avoid SettingWithCopyWarning (FIX)
sens.loc[sens['Name'] == 'TIMM13', 'pval'] = sens.loc[sens['Name'] == 'TIMM13', 'pval'] * 10_000_000
sens.loc[sens['Name'] == 'TIMM9',  'pval'] = sens.loc[sens['Name'] == 'TIMM9',  'pval'] * 400_000

# Precompute log10(P) (note: P should be > 0; NaNs/zeros would become -inf)
sens['logPValue'] = np.log10(sens['pval'])

# --- Biologically-motivated subsetting ---------------------------------------
# Focus on mitochondrial features:
#   - chromosome == "MT" OR GO Cellular Component contains "mitoch"
# Use case-insensitive search and ignore NaNs (na=False) (FIX)
subsetsens = sens.loc[
    (sens['chromosome'] == 'MT') |
    (sens['go.CC'].str.contains('mitoch', case=False, na=False))
].copy()

# Add a label for this cohort
subsetsens['Type'] = "sensitive"

# Ensure numeric dtype (robustness; coerce problematic strings to NaN)
subsetsens[['pval', 'fold_change']] = subsetsens[['pval', 'fold_change']].apply(
    pd.to_numeric, errors='coerce'
)

# Mark significance: p < 0.05 and |fold_change| >= ~1 (asymmetrically defined here to mirror your original)
subsetsens['Significant'] = (
    (subsetsens['pval'] < 0.05) &
    ((subsetsens['fold_change'] > 1) | (subsetsens['fold_change'] < -0.99))
)

# Extract only significant mitochondrial points for highlighting/labels
subsetsens_sub = subsetsens.loc[subsetsens['Significant']].copy().reset_index(drop=True)

# Convenience variables (retain your original names)
concat = subsetsens
concat_sub = concat.loc[concat['Significant']].copy().reset_index(drop=True)

# Sort for potential plotting order (smallest p first, then fold change)
sens = sens.sort_values(by=['pval', 'fold_change'], ascending=[True, True])

# --- Plot layout --------------------------------------------------------------
# Define a compact figure: one panel, tight spacing
colors = ["skyblue", "red"]
color_pal = sns.color_palette(colors, n_colors=2)
heights = [2]
spec = dict(height_ratios=heights, hspace=0.1, wspace=0)

fig, ax1 = plt.subplots(
    nrows=1, ncols=1,
    sharex=True, sharey=True,
    constrained_layout=True,
    figsize=(3.5, 4),
    gridspec_kw=spec
)

# --- Scatter: background (all genes) ------------------------------------------
# Light gray cloud for context; point size scales with fold_change magnitude (as you specified)
g = sns.scatterplot(
    x='fold_change', y='pval',
    data=sens,
    size='fold_change', sizes=(30, 10),      # note: smaller size for smaller |fold_change|
    legend='brief',
    linewidth=0.5, edgecolor='black',
    color='gray', alpha=.1,
    ax=ax1
    # size_norm=(-5, 5),  # you had this commented; include if you want symmetric scaling
)

# --- Scatter: highlighted cohort (significant mitochondrial) ------------------
g = sns.scatterplot(
    x='fold_change', y='pval',
    data=subsetsens_sub,
    size='fold_change', sizes=(100, 200),    # larger points to stand out
    alpha=.95, color='red',
    legend='brief',
    linewidth=0.5, edgecolor='black',
    ax=ax1
)

# --- Axes scaling, limits, and legend -----------------------------------------
# P-values span orders of magnitude; use log y-scale (smaller is better -> lower on plot)
plt.yscale('log')

# Explicit y-limits: top=1, bottom=1e-10 (be sure there are no zeros)
plt.ylim(1, 1e-10)

# Build a compact legend.
# WARNING: The handle indexing below is a bit brittle (depends on seaborn’s legend order).
# It’s left as-is to match your original behavior.
handles, labels = plt.gca().get_legend_handles_labels()
# Make the "red" handle actually appear red if seaborn didn't already
if len(handles) >= 5:
    handles[4].set_color('red')
    handles_new = (handles[4], handles[0])   # ('Mitochondria', 'Other') per your original
    labels_new = ['Mitochondria', 'Other']
    plt.legend(
        handles_new, labels_new, title='',
        bbox_to_anchor=(1.33, 1.12), fontsize=10
    )
else:
    # Fallback: if the automatic legend structure changes in a future seaborn/mpl version
    plt.legend([], [], frameon=False)

# Axis labels and ticks
plt.xlabel(r'Log$_{2}$(B508/Ctrl)', fontsize=18, fontweight='bold', fontname="Arial")
plt.xlim(-4, 4)
plt.xticks([-4, -2, 0, 2, 4], fontsize=12)
plt.ylabel('PValue', fontsize=18, fontweight='bold', fontname="Arial")
plt.yticks(fontsize=12)

# --- Label a subset of the most significant points ----------------------------
# Sort by logPValue (more negative == more significant) and pick top 12
texts = subsetsens_sub.sort_values(by=['logPValue'], ascending=[True]).reset_index(drop=True).head(12)

# Place text labels at their (fold_change, pval) positions
texts1 = [
    plt.text(
        texts.fold_change[i], texts.pval[i], texts.Name[i],
        color='black', fontweight='bold', fontsize=14
    )
    for i in range(len(texts))
]

# Use adjustText to nudge labels apart and draw arrows where needed
adjust_text(
    texts1,
    arrowprops=dict(arrowstyle='->', color='black', alpha=1),
    expand_text=(1, 1)
)

# --- Spine + ticks styling (robust across axes types) -------------------------
all_axes = fig.get_axes()

for ax in all_axes:
    # Hide all spines by default (only if present)
    for side in ('left', 'right', 'top', 'bottom'):
        if side in ax.spines:
            ax.spines[side].set_visible(False)

    # Determine if this axes is at the first row/col of a gridspec in a safe way
    is_first_row = False
    is_first_col = False

    try:
        sps = ax.get_subplotspec()  # works for subplot/grid-spec axes
        if sps is not None:
            is_first_row = sps.is_first_row()
            is_first_col = sps.is_first_col()
    except Exception:
        # Fallback for non-subplot axes (colorbars/insets) or old/new mpl variants
        is_first_row = getattr(ax, "is_first_row", lambda: False)()
        is_first_col = getattr(ax, "is_first_col", lambda: False)()

    # If first row: show bottom spine & simplify ticks on the y-axis
    if is_first_row:
        if 'bottom' in ax.spines:
            ax.spines['bottom'].set_visible(True)
        ax.tick_params(left=False)

    # If first col: show left spine & enable y ticks
    if is_first_col:
        if 'left' in ax.spines:
            ax.spines['left'].set_visible(True)
        ax.tick_params(left=True)

# --- Output -------------------------------------------------------------------
# Save a transparent, high-DPI PNG with tight bounding box and a bit of padding
plt.savefig('PC3_TMM_Mitochondria_xlim4_051624.png', dpi=600, bbox_inches='tight', pad_inches=1, transparent=True)

# Display the figure (useful in notebooks/interactive runs)
plt.show()
