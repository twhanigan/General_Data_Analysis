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
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
rc={"fontname":"Arial",'fontsize':24,'fontweight':'bold','lines.linewidth':2}
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':12,'axes.linewidth':2,'axes.labelsize':14})#,'axes.labelweight':'bold'})
plt.rc('axes', linewidth=2)

# Import Spreadsheets
Raw_Frame = pd.read_csv('abundance_protein_MD_Annotated.csv')

# Get Protein Names from Description
Median_Frame = Raw_Frame.loc[:, "sample-01":"sample-10"].median()
Raw_Frame.loc[:, "sample-01":"sample-10"] = Raw_Frame.loc[:, "sample-01":"sample-10"] - Median_Frame

#check the data integrity using a box and whisker plot
g1 = sns.boxplot(data=Raw_Frame.loc[:, "sample-01":"sample-10"], fill=False, gap=.1)
plt.xticks(rotation=45)

# Add labels and title
plt.title('Normalized Reporter Ion Intensities per Channel')
plt.xlabel('Sample')
plt.ylabel('Log2 Intensity')
plt.tight_layout()
plt.show()

Raw_Frame.fillna(0,inplace=True)
Raw_Frame.loc[:, "sample-01":"sample-10"] = Raw_Frame.loc[:, "sample-01":"sample-10"].apply(pd.to_numeric)
Raw_Frame = Raw_Frame.set_index('Gene')
Raw_Frame.replace(np.inf, np.nan, inplace =True)
Raw_Frame.replace(-np.inf, np.nan, inplace =True)
Raw_Frame = Raw_Frame.dropna()
Raw_Frame = Raw_Frame.drop_duplicates()
Raw_Frame['mito_ribosome'] = 0
Raw_Frame.loc[Raw_Frame['go.CC'].str.contains('mitochondrion', na=False), 'mito_ribosome'] = 1
Raw_Frame.loc[Raw_Frame.index.str.contains('MRPL', na=False), 'mito_ribosome'] = 2
Raw_Frame.loc[Raw_Frame.index.str.contains('MRPS', na=False), 'mito_ribosome'] = 3
print(Raw_Frame['mito_ribosome'].loc[Raw_Frame['mito_ribosome']==1])
corrMatrix = Raw_Frame.loc[:, "sample-01":"sample-10"].corr()

# Dynamically calculate difference columns relative to 'sample-01'
for i in range(1, 11):
    col = f'sample-{i:02d}'
    diff_col = f'diff_S{i}'
    Raw_Frame[diff_col] = Raw_Frame[col] - Raw_Frame['sample-01']

#Find ratio of B819 to inactive control protein and take only rows with B508/B143 greater than 2
one =  Raw_Frame["sample-01"],Raw_Frame['sample-02'],Raw_Frame["sample-03"],Raw_Frame['sample-04'],Raw_Frame["sample-05"]
two = Raw_Frame["sample-06"],Raw_Frame['sample-07'],Raw_Frame["sample-08"],Raw_Frame['sample-09'],Raw_Frame["sample-10"]
for row in one:
    t, p = ttest_ind(one, two)
    Raw_Frame['4S_5S_ttest'] = t
    Raw_Frame['4S_5S_pval'] = p

#Raw_Frame = Raw_Frame.sort_values(by=['diff_S3',"diff_S8"],ascending =[False,False])#ascending =[False,True,True,True,True]
Raw_Frame = Raw_Frame.sort_values(by=['mito_ribosome'], ascending =[False])#ascending =[False,True,True,True,True]
row_colors = Raw_Frame['mito_ribosome'].map({0: 'lightgrey', 1: 'red',2: 'darkred',3:'crimson'})
hits_pbs = Raw_Frame[Raw_Frame['4S_5S_pval']<0.1]#&(Raw_Frame['diff_525']>2.9)&(Raw_Frame['diff_526']>1)&(Raw_Frame['diff_888']>1)]
s1 = (hits_pbs)
print(row_colors.unique())
print(Raw_Frame['mito_ribosome'].value_counts(dropna=False))
#Plot the figure
g = sns.clustermap(
    data=Raw_Frame.loc[:, "diff_S1":"diff_S10"],
    cmap="vlag",
    cbar_kws={'label': r'Input Normalized Reporter Ion Intensity'},
    figsize=(6, 6.5),
    col_cluster=False,
    row_cluster=False,
    z_score=0,
    method='centroid',
    dendrogram_ratio=0.15,
    cbar_pos=(-0.2, .2, .03, .4),
    row_colors=row_colors 
)
xlabel = ['4S Lysate','4S_Mito','4S_F1','4S_F6','4S_F12','5S Lysate','5S_Mito','5S_F1','5S_F6','5S_F12']
g.ax_heatmap.tick_params(length=12,width=1,which='major',axis='y')
g.ax_heatmap.tick_params(length=6,width=1,which='minor',axis='x')
g.ax_heatmap.tick_params(length=10,width=2,which='major',axis='x',pad=0)
for _, spine in g.ax_heatmap.spines.items():
    spine.set_visible(True)
plt.savefig('HeatMap_mitoribo_fractions_A.png',dpi=600,transparent=True,bbox_inches='tight',pad_inches=2)
plt.show()




