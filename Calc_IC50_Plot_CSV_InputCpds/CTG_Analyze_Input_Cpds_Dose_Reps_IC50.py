# Analyze CTG output (UPDATED: IC50 extraction + averaging + Excel export + IC50-annotated plots)
import os
import numpy as np
import pandas as pd
from itertools import cycle
import matplotlib.pyplot as plt
from lmfit import Model

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
my_path = r"C:\\Users\\thanigan\\Downloads\\Win_CTG_Data_11182025"

# Outputs
OUTPUT_EXCEL = "CTG_IC50_Output.xlsx"
PLOT_DIR     = "CTG_Plots"


def get_filepaths(directory):
    """Get full filepaths of all files in a directory, including sub-directories.
       If a single file is passed, return [that file].
    """
    if os.path.isfile(directory):
        return [directory]
    file_paths = []
    for root, _, files in os.walk(directory, topdown=True):
        for filename in files:
            filepath = os.path.join(root, filename)
            if filepath.lower().endswith(".csv"):
                file_paths.append(filepath)
    return file_paths


# ---------- helper functions ----------
def generate_concentrations(max_conc, number_dilutions, dilution_factor):
    """
    Concentration list needs to be [0, max, max/DF, max/DF^2, ..., lowest]
    """
    num_nonzero = number_dilutions - 1
    conc_desc = [max_conc / (dilution_factor ** i) for i in range(num_nonzero)]
    conc_list = [0.0] + conc_desc
    log_conc = [np.nan] + [np.log10(c) for c in conc_desc]
    conc_name = [f"C{i+1}" for i in range(number_dilutions)]
    return conc_list, np.array(log_conc, dtype=float), conc_name


def sigmoid(x, b, c, d, e):
    # e is IC50 in the model parameterization
    return (c + (d - c) / (1 + np.exp(b * (np.log(e) - np.log(x)))))


def pDose(x):
    return np.log10(x)


def normalize(frame, frame_name, negative_control):
    """
    Normalized = (Response - val_min) / (val_max - val_min)
      where val_max is mean response at lowest concentration (no-drug)
      and val_min is plate negative control (NEG_CTRL_ROW average across dose columns)
    """
    frame = frame.copy()
    for col in ['Concentration', 'Response']:
        frame[col] = pd.to_numeric(frame[col], errors='coerce')

    compoundData = frame.groupby(['compound'], sort=False)
    for _, group in compoundData:
        val_max = group.loc[group['Concentration'] == group['Concentration'].min(), 'Response'].mean()
        val_min = negative_control
        frame.loc[group.index, 'Normalized'] = (frame.loc[group.index, 'Response'] - val_min) / (val_max - val_min)

    frame['Normalized'] = pd.to_numeric(frame['Normalized'], errors='coerce')
    return frame


def _initial_guesses_from_group(group):
    """Heuristic initial guesses for the 4-param sigmoid, robust to a 0 (no-drug) column."""
    g = group[['Concentration', 'Normalized']].dropna().copy()
    g = g.sort_values('Concentration')

    if g.empty or g['Concentration'].nunique() < 3:
        xmid = np.median(g['Concentration']) if not g.empty else 1e-7
        return float(xmid), 1.0, 1.0, 0.0

    g_pos = g.copy()
    if (g_pos['Concentration'] > 0).any():
        min_pos = g_pos.loc[g_pos['Concentration'] > 0, 'Concentration'].min()
    else:
        min_pos = 1e-7
    g_pos.loc[g_pos['Concentration'] <= 0, 'Concentration'] = min_pos / 10.0

    d_top = g['Normalized'].loc[g['Concentration'] == g['Concentration'].min()].mean()
    c_bot = g['Normalized'].loc[g['Concentration'] == g['Concentration'].max()].mean()

    x_for_trend = np.log10(g_pos['Concentration'].values)
    y_for_trend = g_pos['Normalized'].values
    try:
        trend = np.corrcoef(x_for_trend, y_for_trend)[0, 1]
    except Exception:
        trend = -1.0

    if np.isfinite(trend) and trend > 0:
        d_top, c_bot = (
            g['Normalized'].loc[g['Concentration'] == g['Concentration'].max()].mean(),
            g['Normalized'].loc[g['Concentration'] == g['Concentration'].min()].mean()
        )

    y_mid = (d_top + c_bot) / 2.0
    idx_mid = (g_pos['Normalized'] - y_mid).abs().idxmin()
    init_ic50 = float(g_pos.loc[idx_mid, 'Concentration'])

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

    init_b = float(np.clip(init_b, -3.0, 3.0))
    if abs(init_b) < 1e-3:
        init_b = 0.5 if (np.isfinite(trend) and trend < 0) else -0.5

    return float(init_ic50), float(init_b), float(d_top), float(c_bot)


def fit_curves_and_ic50(frame):
    """
    Fit sigmoid per (compound, frame_name), return:
      - fitted_params_df: includes IC50 and fit params
      - evaluated_frame: dense curve grid for plotting/export
    """
    frame = frame.copy()
    frame['Concentration'] = pd.to_numeric(frame['Concentration'], errors='coerce')
    frame['Normalized']    = pd.to_numeric(frame['Normalized'], errors='coerce')
    frame = frame.dropna(subset=['Concentration', 'Normalized'])

    compoundData = frame.groupby(['compound', 'frame_name'], sort=False)

    evaluated_frame = pd.DataFrame(
        columns=['dose', 'new', 'dely', 'upperbound', 'lowerbound', 'compound', 'group_num', 'frame_name'],
        dtype=float
    )

    fitted_rows = []
    nums = [0, 1, 2]
    licycle = cycle(nums)

    gmodel = Model(sigmoid)

    for (compound, frame_name), group in compoundData:
        init_ic50, init_slope, val_max, val_min = _initial_guesses_from_group(group)
        decreasing = (val_max > val_min)

        params = gmodel.make_params()
        params.add('d', min=0.7 if decreasing else -0.1, max=1.5 if decreasing else 2.0)
        params.add('c', min=-0.1, max=0.7 if decreasing else 2.0)
        params.add('e', min=max(frame['Concentration'].min(), 1e-18) * 1e-3,
                   max=max(frame['Concentration'].max(), 1e-18) * 1e3)
        params.add('b', min=-3, max=3)

        f = next(licycle)
        xdata = group['Concentration'].astype(float).values
        ydata = group['Normalized'].astype(float).values

        # evaluation grid uses positive doses only
        if (xdata > 0).any():
            refDose_min = np.nanmin(xdata[xdata > 0])
        else:
            refDose_min = 1e-12
        refDose_max = np.nanmax(xdata) if np.isfinite(np.nanmax(xdata)) else refDose_min * 10.0
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

        # keep your existing band logic
        dely = 0.3 * max(1e-6, float(val_min))
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
                'compound': compound,
                'group_num': f,
                'frame_name': frame_name
            })
        ], ignore_index=True)

        # ---- extract IC50 (parameter e) + stderr if available ----
        ic50_val = float(result.params['e'].value)
        ic50_se  = result.params['e'].stderr
        ic50_se  = float(ic50_se) if ic50_se is not None and np.isfinite(ic50_se) else np.nan

        fitted_rows.append({
            'frame_name': frame_name,
            'compound': compound,
            'IC50': ic50_val,
            'IC50_SE': ic50_se,
            'b_slope': float(result.params['b'].value),
            'c_bottom': float(result.params['c'].value),
            'd_top': float(result.params['d'].value),
            'fit_success': bool(getattr(result, "success", True)),
            'n_points': int(np.isfinite(xdata).sum())
        })

    fitted_params_df = pd.DataFrame(fitted_rows)
    return fitted_params_df, evaluated_frame


def summarize_ic50_across_replicates(ic50_per_file_df):
    """
    Average IC50 across biological replicates (here: across files / frame_name) for each compound.
    Outputs mean / SD / SEM / N, plus pIC50 (using mean IC50) for convenience.
    """
    df = ic50_per_file_df.copy()
    df['IC50'] = pd.to_numeric(df['IC50'], errors='coerce')
    df = df.dropna(subset=['IC50'])

    summary = (df.groupby('compound', as_index=False)
                 .agg(
                     N=('IC50', 'count'),
                     IC50_mean=('IC50', 'mean'),
                     IC50_sd=('IC50', 'std')
                 ))
    summary['IC50_sem'] = summary['IC50_sd'] / np.sqrt(summary['N'].clip(lower=1))
    summary['pIC50_from_meanIC50'] = -np.log10(summary['IC50_mean'])
    return summary


def plot_curves_with_ic50(raw_norm_df, evaluated_df, ic50_df, outdir):
    """
    One plot per frame_name. Plots raw normalized points and fitted curves.
    Annotates IC50 with a vertical line and text.
    """
    os.makedirs(outdir, exist_ok=True)

    raw = raw_norm_df.copy()
    raw['Concentration'] = pd.to_numeric(raw['Concentration'], errors='coerce')
    raw['Normalized'] = pd.to_numeric(raw['Normalized'], errors='coerce')

    for frame_name in sorted(raw['frame_name'].unique()):
        sub_raw = raw[raw['frame_name'] == frame_name].copy()
        sub_eval = evaluated_df[evaluated_df['frame_name'] == frame_name].copy()
        sub_ic50 = ic50_df[ic50_df['frame_name'] == frame_name].copy()

        if sub_raw.empty:
            continue

        fig, ax = plt.subplots(figsize=(6.5, 4.8))

        # determine a plotting x for the "0" dose point (log-safe)
        pos_concs = sub_raw.loc[sub_raw['Concentration'] > 0, 'Concentration']
        min_pos = float(pos_concs.min()) if not pos_concs.empty else 1e-12
        zero_plot_x = min_pos / 10.0

        for compound in sub_raw['compound'].unique():
            r = sub_raw[sub_raw['compound'] == compound].copy()
            r = r.dropna(subset=['Concentration', 'Normalized'])

            # plot raw points (replace 0 with zero_plot_x for visibility on log axis)
            x_plot = r['Concentration'].values.astype(float)
            x_plot = np.where(x_plot <= 0, zero_plot_x, x_plot)
            y_plot = r['Normalized'].values.astype(float)
            ax.semilogx(x_plot, y_plot, marker='o', linestyle='None', label=f"{compound} data")

            # plot fit curve
            e = sub_eval[sub_eval['compound'] == compound].copy()
            if not e.empty:
                x_fit = (10 ** e['dose'].astype(float).values)
                y_fit = e['new'].astype(float).values
                ax.semilogx(x_fit, y_fit, label=f"{compound} fit")

            # annotate IC50
            ic = sub_ic50[sub_ic50['compound'] == compound]
            if not ic.empty and np.isfinite(ic['IC50'].iloc[0]):
                ic50_val = float(ic['IC50'].iloc[0])
                ax.axvline(ic50_val, linewidth=1)
                ax.text(ic50_val, 1, f"IC50={ic50_val:.3g}", rotation=00, va='center', ha='right')

        ax.set_title(frame_name)
        ax.set_xlabel("Concentration")
        ax.set_ylabel("Normalized Response")
        ax.set_ylim(-0.1, 1.2)
        ax.legend(fontsize=8)

        outpath = os.path.join(outdir, f"{os.path.splitext(frame_name)[0]}_IC50.png")
        fig.tight_layout()
        fig.savefig(outpath, dpi=300)
        plt.close(fig)


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

# =========================
# main
# =========================
files = get_filepaths(my_path)
print("Files found:", files)

frame_list = []

for file in files:
    frame = pd.read_csv(file, skiprows=0, index_col=0, header=0)

    # keep only first NUM_DOSES columns (including no-drug col 1)
    frame = frame[[str(i) for i in range(1, NUM_DOSES + 1)]]
    frame = frame.apply(pd.to_numeric, errors='coerce')

    frame_name = str(os.path.basename(file))

    # per-plate negative control (row NEG_CTRL_ROW, across all dose columns)
    negative_control = float(frame.loc[NEG_CTRL_ROW, :].mean())

    # build blocks for each compound
    lst = []
    for comp_idx in range(NUM_COMPOUNDS):
        comp_name = Compound_Names[comp_idx]
        rows_for_comp = row_labels[start_idx + comp_idx * REPS_PER_COMPOUND:
                                   start_idx + (comp_idx + 1) * REPS_PER_COMPOUND]
        rows_for_comp = [r for r in rows_for_comp
                         if r in row_labels and row_labels.index(r) != ctrl_idx]

        block = frame.loc[rows_for_comp].astype(float)
        avg = block.mean(axis=0)
        err = block.std(axis=0)

        conc_list, log_conc, _ = conc_info[comp_idx]
        for c_idx, col in enumerate(block.columns):
            x_value = float(conc_list[c_idx])
            log_x   = float(log_conc[c_idx]) if np.isfinite(log_conc[c_idx]) else np.nan
            average = float(avg[col])
            error   = float(err[col])
            lst.append([frame_name, comp_name, x_value, log_x, average, error])

    final = pd.DataFrame(
        lst,
        columns=['frame_name', 'compound', 'Concentration', 'Log Concentration', 'Response', 'Error']
    )

    for col in ['Concentration', 'Log Concentration', 'Response', 'Error']:
        final[col] = pd.to_numeric(final[col], errors='coerce')

    normal = normalize(final, frame_name=frame_name, negative_control=negative_control)
    frame_list.append(normal)

# Combine all normalized frames
total_frame = pd.concat(frame_list, ignore_index=True)

# ---- fit + IC50 extraction ----
ic50_per_file, evaluated = fit_curves_and_ic50(total_frame)

# ---- average IC50 across biological replicates (files) ----
ic50_summary = summarize_ic50_across_replicates(ic50_per_file)

# ---- plots with IC50 annotation ----
plot_curves_with_ic50(
    raw_norm_df=total_frame,
    evaluated_df=evaluated,
    ic50_df=ic50_per_file,
    outdir=PLOT_DIR
)

# ---- Excel export with separate sheets ----
with pd.ExcelWriter(OUTPUT_EXCEL, engine='openpyxl') as writer:
    total_frame.to_excel(writer, sheet_name="Raw_Normalized_Data", index=False)
    evaluated.to_excel(writer, sheet_name="Evaluated_Curves", index=False)
    ic50_per_file.to_excel(writer, sheet_name="IC50_Per_File", index=False)
    ic50_summary.to_excel(writer, sheet_name="IC50_Summary", index=False)

print(f"Wrote Excel: {OUTPUT_EXCEL}")
print(f"Wrote plots to: {PLOT_DIR}")