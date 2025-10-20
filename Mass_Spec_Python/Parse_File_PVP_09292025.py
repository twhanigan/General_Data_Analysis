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
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as pyplot
from pylab import *

rc={"fontname":"Arial",'fontsize':24,'fontweight':'bold','lines.linewidth':2}
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':12,'axes.linewidth':2,'axes.labelsize':14})#,'axes.labelweight':'bold'})
plt.rc('axes', linewidth=2)

# Import Spreadsheets
Raw_Frame = pd.read_excel('052121_PvP_S1.xlsx')
Raw_Frame.median('sample-01')
print(Raw_Frame.median('sample-01'))
print('hello world')
hits_pbs_b = pd.read_excel('hits_S2.xlsx')
hits_pbs_b.index = hits_pbs_b.Name
Mitorib = ['MCU','AGK','PHB','GADD45GIP1','ERAL1','TBRG4','MRPL11','MRPL51','ERAL1','MRPL4','MRPL21','MRPS18A','SSBP1','IMMT','MRPL15','MRPL30','MRPS30','MRPL3','MRPL39','MRPL19','MRPL32','MRPL27','OXA1L','MRPL45','MRPL24','ELAC2','MRPL2','FASTKD2','MRPL37','MRP63','DDX5','MRPL28']
OXA1L = ['OXA1L']
Cyto = ['NQO1','SLC16A3','KDELR1']
ComplexIV = ['MT-CO1','MT-CO2','MT-CO3','COX5B','COX5A','COX6C','COX7C','MT-CYB','COX6B1','COX7A2','PET100'] #'OXA1L',
ComplexI = ['NDUFC1','NDUFB8','NDUFB6','NDUFB3','NDUFA1','NDUFA12','NDUFS5','NDUFA2','NDUFA8','NDUFA6','NDUFAB1','NDUFV2','NDUFA11','NDUFS4','NDUFS6','NDUFS5','NDUFA5','NDUFB10','NDUFA10','ACAD9',]
ComplexII = ['SDHD','SDHA','SDHC','SDHAF2','FH','NDUFV1','UQCRFS1','NDUFS2','SUCLG1','SUCLA2','SDF2']
ComplexIII = ['UQCRB','MT-CYB','UQCRC1','CYC1','UQCRH','NDUFV1','UQCRC2','UQCR10','UQCRFS1','UQCRQ','UQCR11']
ComplexV = ['SLC25A24','ATP13A1','SLC25A5', 'SLC25A6','ATP5G1','ATP5G3','ATP5F1','ATP5I','ATP5D','ATP5J2','ATP5O','ATP5L','ATP5B','ATP5H','ATP5C1','ATP5A1','ATP6V1H','MT-ATP6','ATP5I','ATP5H','ATP2A2']
UreaCycle = ['ASS1','PYCR2','PYCR1','ALDH18A1','ATP1B3','CPS1','SLC25A15','OAT']
Thymidine = ['DTYMK','TK2','TK1','CDA','NTPCR','CAPN2','LDHA']

# log2 transform
Log2_Frame = Raw_Frame.loc[:, "DMSO_R1":"BMT174819_R2"].copy()
Log2_Frame = np.log2(
    Raw_Frame.loc[:, "DMSO_R1":"BMT174819_R2"]).copy()

Raw_Frame['description'] = Raw_Frame['description'].astype(str)
Log2_Frame['Name'] = Raw_Frame['description'].str.split(" ", n=1).str.get(0)
Median_Frame = Log2_Frame.loc[:, Log2_Frame.columns != 'Name'].median()
Probe_Targets = Log2_Frame.loc[:,"DMSO_R1":"BMT174819_R2"] #- Median_Frame

Probe_Targets['Name']=Log2_Frame.Name
Competed = []
Probe_Targets.fillna(0,inplace=True)
Probe_Targets.loc[:,"DMSO_R1":"BMT174819_R2"] = Probe_Targets.loc[:,"DMSO_R1":"BMT174819_R2"].apply(pd.to_numeric)
Probe_Targets = Probe_Targets.set_index('Name')
Probe_Targets.replace(np.inf, np.nan, inplace =True)
Probe_Targets.replace(-np.inf, np.nan, inplace =True)
Probe_Targets = Probe_Targets.dropna()
Probe_Targets = Probe_Targets.drop_duplicates()

#Average each replicate
Probe_Targets['DMSO Avg'] = Probe_Targets[["DMSO_R1", "DMSO_R2"]].mean(axis = 1)
Probe_Targets['BMT181525 Avg'] = Probe_Targets[["BMT181525_R1", 'BMT181525_R2']].mean(axis = 1)
Probe_Targets['BMT182526 Avg'] = Probe_Targets[["BMT182526_R1", 'BMT182526_R2']].mean(axis = 1)
Probe_Targets['BMT179888 Avg'] = Probe_Targets[["BMT179888_R1", 'BMT179888_R2']].mean(axis = 1)
Probe_Targets['BMT174819 Avg'] = Probe_Targets[["BMT174819_R1", 'BMT174819_R2']].mean(axis = 1)
Probe_Targets['DMSO Avg'] = Probe_Targets[["DMSO_R1", "DMSO_R1"]].mean(axis = 1)
Probe_Targets['diff_DMSO'] = Probe_Targets['BMT174819 Avg']-Probe_Targets['DMSO Avg']
Probe_Targets['diff_525'] = Probe_Targets['BMT174819 Avg']-Probe_Targets['BMT181525 Avg']
Probe_Targets['diff_526'] = Probe_Targets['BMT174819 Avg']-Probe_Targets['BMT182526 Avg']
Probe_Targets['diff_888'] = Probe_Targets['BMT174819 Avg']-Probe_Targets['BMT179888 Avg']

#Find ratio of B819 to inactive control protein and take only rows with B508/B143 greater than 2
one =  Probe_Targets["BMT174819_R1"],Probe_Targets['BMT174819_R2']
two = Probe_Targets["BMT181525_R1"],Probe_Targets['BMT181525_R2'],Probe_Targets["BMT182526_R1"],Probe_Targets['BMT182526_R2'],Probe_Targets["BMT179888_R1"],Probe_Targets["BMT179888_R2"],Probe_Targets["DMSO_R1"],Probe_Targets['DMSO_R2']
for row in one:
    t, p = ttest_ind(one, two)
    Probe_Targets['BMT174819_ttest'] = t
    Probe_Targets['BMT174819_pval'] = p
Probe_Targets = Probe_Targets.sort_values(by=['BMT174819 Avg','BMT182526 Avg','BMT181525 Avg','BMT179888 Avg','DMSO Avg'],ascending =[False,True,True,True,True])#ascending =[False,True,True,True,True]

#hits_pbs = Probe_Targets.loc[(Probe_Targets['diff_DMSO']>4)&(Probe_Targets['diff_525']>2.9)&(Probe_Targets['diff_526']>1)&(Probe_Targets['diff_888']>1)]
hits_pbs = Probe_Targets[Probe_Targets['BMT174819_pval']<0.004]#&(Probe_Targets['diff_525']>2.9)&(Probe_Targets['diff_526']>1)&(Probe_Targets['diff_888']>1)]
hits_pbs_b = hits_pbs_b[hits_pbs_b['BMT174819_pval']<0.004]#&(Probe_Targets['diff_525']>2.9)&(Probe_Targets['diff_526']>1)&(Probe_Targets['diff_888']>1)]

Probe_Targetsb = Probe_Targets.copy()
Probe_Targetsb['DMSO_R2'] = Probe_Targets['DMSO Avg']
s1 = pd.concat((hits_pbs, hits_pbs_b), join="inner")
output = s1.groupby(s1.index).mean().sort_values(by=['BMT174819 Avg','BMT182526 Avg','BMT181525 Avg','BMT179888 Avg','DMSO Avg'],ascending =[False,True,True,True,True]).drop_duplicates()
output = output.loc[~output.index.str.contains("TUBB|NEF")]

Probe_Targets.replace(np.inf, np.nan, inplace =True)
Probe_Targets.replace(-np.inf, np.nan, inplace =True)
Probe_Targets = Probe_Targets.sort_values(by=['BMT174819 Avg','BMT182526 Avg','BMT181525 Avg','BMT179888 Avg','DMSO Avg'],ascending =[False,True,True,True,True])#ascending =[False,True,True,True,True]
Probe_Targets.replace(np.nan, np.nan, inplace=True)
Probe_Targets.dropna(inplace=True)
Probe_Targets.to_excel('Analyzed_B819_PVP_S2_09292025.xlsx')
output.to_excel('Combined_hitsc.xlsx')

#Make heatmap of subset
heatmap_set = Probe_Targets.iloc[:,1:10].copy()
heatmap_set.dropna(inplace=True)
xlabel = ['DMSO','BMT181525','BMT182526','BMT179888','BMT174819']
colors = ["skyblue", "red",'red' ]
color_pal = sns.color_palette(colors, n_colors=2)
#g=sns.clustermap(heatmap_set,cmap = "RdYlBu", col_cluster= True,figsize = (7,7),method='centroid') #col_cluster=False #standard_scale=1, metric = 'correlation' 
g=sns.clustermap(Probe_Targets[['DMSO Avg','BMT181525 Avg','BMT182526 Avg','BMT179888 Avg','BMT174819 Avg']],cmap = "RdYlBu",cbar_kws={'label':r'Log$_{2}$ Reporter Ion Intensity'},figsize=(5.35,6),cbar_pos=(1.3, 1.1,0.04,0.2), col_cluster= True,row_cluster=False,method = 'average') #vmin=-2, vmax=6,standard_scale=1,,linewidths=0.00001,linecolor='white', method='median'col_cluster=False #standard_scale=1, metric = 'correlation' 
#g=sns.clustermap(output[['BMT181525 Avg','DMSO Avg','BMT182526 Avg','BMT179888 Avg','BMT174819 Avg']],cmap = "Reds",cbar_kws={'label':r'Log$_{2}$ Reporter'+'\nIon Intensity',"orientation": "horizontal",'extend':'both','ticks':[10,15,20]},figsize=(2,6.4),cbar_pos=(1.3, 1.1,0.5,0.03), col_cluster= False,row_cluster=False,method = 'average',vmin=10, vmax=20) #vmin=-2, vmax=6,standard_scale=1,,linewidths=0.00001,linecolor='white', method='median'col_cluster=False #standard_scale=1, metric = 'correlation' 

pos= g.ax_heatmap.get_yticks()
labels = g.ax_heatmap.get_yticklabels()
print(labels)
labelsnew =  labels[:-3]+['OXA1L']+labels[-3:-1]
print(len(labelsnew))

labelsb = ['SLC16A3','MT-CO2','NQO1',  'KDELR1','ATP5L', 'COX7A2', 'SLC25A3', 'ECH1', 'HSD17B12', 'ATP5F1','MTCH1', 'MTCH2', 'ATP13A1', 'ATP5G1', 'SLC25A10', 'NDUFB3', 'OXA1L', 'ACAD9',  'SLC25A24']
colors = ['lightcoral','tomato','red','crimson','darkred','k',]
newpos = pos[0::1]
labelsnew = ['VDAC2','OXA1L'] + labels[2:]
plt.setp(g.ax_heatmap.set_yticks(newpos))
plt.setp(g.ax_heatmap.set_yticklabels(labelsnew),fontweight='bold')
for _, spine in g.ax_heatmap.spines.items():
    spine.set_visible(True)

plt.setp(g.ax_heatmap.set_ylabel(''),fontsize=18,fontweight='bold',fontname='Arial',ha='center')
ticks = g.ax_heatmap.get_yticklabels()
print(ticks)
for tick in ticks:
    print(tick)
    if tick.get_text() in ComplexI:
        tick.set_color('lightcoral')
    elif tick.get_text() in ComplexII:
        tick.set_color('tomato')
    elif tick.get_text() in ComplexIII:
        tick.set_color('red')
    elif tick.get_text() in ComplexIV:
        tick.set_color('crimson')
    elif tick.get_text() in ComplexV:
        tick.set_color('darkred')
    elif tick.get_text() in Cyto:
        tick.set_color('gray')
    elif tick.get_text() in OXA1L:
        tick.set_color('steelblue')
    else:
        tick.set_color('black')

handle =  [Line2D([0], [0], color='tomato',markerfacecolor='tomato',alpha=0.5, markersize=7),Line2D([0], [0], color='red',markerfacecolor='red',alpha=0.7, markersize=7),
    Line2D([0], [0], color='crimson',markerfacecolor='crimson',alpha=0.25, markersize=7),Line2D([0], [0], color='darkred',markerfacecolor='darkred',alpha=0.7, markersize=7),Line2D([0], [0], color='gray',markerfacecolor='gray',alpha=0.25, markersize=7),Line2D([0], [0],color='blue',markerfacecolor='blue',alpha=0.7, markersize=7),Line2D([0], [0], color='black',markerfacecolor='black',alpha=0.25, markersize=7),Line2D([0], [0], color='gray',markerfacecolor='gray',alpha=0.25, markersize=7)]#Line2D([0], [0],color='gray',markerfacecolor='gray',alpha=0.7, markersize=7)]
label = ['Complex I','Complex II','Complex III','Complex IV','Complex V','Mitochondrial']

leg = plt.legend(handle,label, handlelength=0, handletextpad=0,bbox_to_anchor=(5, 1, 1, 1), loc='upper right',borderaxespad=0.,ncol=3,columnspacing=0.5,)

i=-1
for text in leg.get_texts():
    i+=1
    text.set_color(colors[i])
plt.savefig('HeatMap_PVP_HitsOnly_Legend_09292025.png',dpi=600,transparent=True,bbox_inches='tight',pad_inches=2)
plt.show()