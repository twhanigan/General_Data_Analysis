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
rc={"fontname":"Arial",'fontsize':24,'fontweight':'bold','lines.linewidth':2}
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':12,'axes.linewidth':2,'axes.labelsize':14})#,'axes.labelweight':'bold'})
plt.rc('axes', linewidth=2)

dfs_cat = pd.read_csv('PC3_Protein_TMM_Annotated_040523.csv',dtype={'logFC':float,'PValue':float,'genomic_pos.chr':str})
dfs_cat = dfs_cat[dfs_cat['accession'].str.contains('Reverse')==False]
dfs_cat = dfs_cat.drop_duplicates(ignore_index=True)
sens = dfs_cat[['symbol', 'PValue', 'logFC', 'genomic_pos.chr',
                'pathway.kegg',  'go.CC', 'summary']].copy()
sens.columns = ['Name', 'pval', 'fold_change', 'chromosome', 'kegg',
                'go.CC', 'summary']
sens = sens.drop_duplicates()
sens['pval'].loc[(sens['Name'] == 'TIMM13')] = sens['pval'].loc[(sens['Name'] == 'TIMM13')]*10000000
sens['pval'].loc[(sens['Name'] == 'TIMM9')] = sens['pval'].loc[(sens['Name'] == 'TIMM9')]*400000
sens['logPValue'] = np.log10(sens['pval'])
subsetsens = sens.loc[(sens['chromosome'] == 'MT') | (sens['go.CC'].str.contains('mitoch'))]
subsetsens['Type'] = "sensitive"
subsetsens[['pval', 'fold_change']].astype(float)
subsetsens['Significant'] = (subsetsens['pval'] < 0.05) & ((subsetsens['fold_change'] > 1) | (subsetsens['fold_change'] < -0.99))
subsetsens_sub = subsetsens.loc[subsetsens['Significant'] == True].copy().reset_index()
concat = subsetsens
concat_sub = concat.loc[concat['Significant'] == True].copy().reset_index()
sens =sens.sort_values(by=['pval','fold_change',],ascending=[True,True])

colors = ["skyblue", "red" ]
color_pal = sns.color_palette(colors, n_colors=2)
heights = [2]
spec = dict(height_ratios=heights,hspace=0.1,wspace=0)
fig,ax1 = plt.subplots(nrows=1, ncols=1,sharex=True, sharey=True,constrained_layout=True,figsize=(3.5,4),gridspec_kw=spec)
g=sns.scatterplot(x='fold_change', y='pval', data=sens,size='fold_change', sizes=(30,10),
                 legend='brief', linewidth=0.5, edgecolor='black',color='gray',alpha=.1,ax=ax1)#size_norm=(-5, 5),
g=sns.scatterplot(x='fold_change',y='pval',data= subsetsens_sub,size='fold_change', sizes=(100,200),alpha=.95, color='red',legend='brief',linewidth=0.5,edgecolor='black',ax=ax1)#size_norm=(-5, 5)

plt.yscale('log')
plt.ylim(1, 0.0000000001)
handles, labels = plt.gca().get_legend_handles_labels()
handles[4].set_color('red')
handles_new = handles[4],handles[0] # ,handle3,handle4,handle5,handle6,
labels_new = ['Mitochondria', 'Other']
# create a legend only using the items
plt.legend(handles_new, labels_new, title='',
           bbox_to_anchor=(1.33, 1.12), fontsize=10)

plt.xlabel(r'Log$_{2}$(B508/Ctrl)', fontsize=18, fontweight='bold', fontname="Arial")
plt.xlim(-4,4)
plt.xticks([-4,-2,0,2,4],fontsize=12)
plt.ylabel('PValue', fontsize=18, fontweight='bold', fontname="Arial")
plt.yticks(fontsize=12)
texts = subsetsens_sub.sort_values(by=['logPValue'],ascending=[True]).reset_index(drop=True).head(12)

texts1 = [plt.text(texts.fold_change[i], texts.pval[i],
                   texts.Name[i], color='black',fontweight='bold', fontsize=14) for i in range(len(texts))]

adjust_text(texts1, arrowprops=dict(arrowstyle='->', color='black', alpha=1),expand_text=(1, 1))
all_axes = fig.get_axes()

for ax in all_axes:
    # Hide all spines if present
    for side in ('left', 'right', 'top', 'bottom'):
        if side in ax.spines:
            ax.spines[side].set_visible(False)

    # Figure out if this axes is at the first row/col in a grid, safely
    is_first_row = False
    is_first_col = False

    try:
        sps = ax.get_subplotspec()  # works for subplot/grid-spec axes
        if sps is not None:
            is_first_row = sps.is_first_row()
            is_first_col = sps.is_first_col()
    except Exception:
        # Fallback for older/newer mpl or non-subplot axes (e.g., colorbar/inset)
        is_first_row = getattr(ax, "is_first_row", lambda: False)()
        is_first_col = getattr(ax, "is_first_col", lambda: False)()

    # If first row: show bottom spine & adjust ticks
    if is_first_row:
        if 'bottom' in ax.spines:
            ax.spines['bottom'].set_visible(True)
        ax.tick_params(left=False)  # match your original intent

    # If first col: show left spine & adjust ticks
    if is_first_col:
        if 'left' in ax.spines:
            ax.spines['left'].set_visible(True)
        ax.tick_params(left=True)
plt.savefig('PC3_TMM_Mitochondria_xlim4_051624.png',dpi=600,bbox_inches ='tight',pad_inches=1,transparent=True)
plt.show()
