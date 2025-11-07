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
rc={"fontname":"Arial",'fontsize':24,'fontweight':'bold','lines.linewidth':2}
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':12,'axes.linewidth':2,'axes.labelsize':14})#,'axes.labelweight':'bold'})
plt.rc('axes', linewidth=2)

# Import Spreadsheets
Raw_Frame = pd.read_excel('052121_PvP_S1.xlsx')
hits_pbs_b = pd.read_excel('hits_S2.xlsx')
hits_pbs_b.index = hits_pbs_b.Name

print(Raw_Frame.head())
# log2 transform
Log2_Frame = Raw_Frame.loc[:, "DMSO_R1":"BMT174819_R2"].copy()
Log2_Frame = np.log2(
    Raw_Frame.loc[:, "DMSO_R1":"BMT174819_R2"]).copy()
#Log2_Frame.loc[:, "DMSO_R1":"BMT174819_R2"] = Log2_Frame.loc[:, "DMSO_R1":"BMT174819_R2"].div(Log2_Frame.sum(axis=1), axis=0)
# Get Protein Names from Description
Raw_Frame['description'] = Raw_Frame['description'].astype(str)
Log2_Frame['Name'] = Raw_Frame['description'].str.split(" ", n=1).str.get(0)
Median_Frame = Log2_Frame.loc[:, Log2_Frame.columns != 'Name'].median()
Probe_Targets = Log2_Frame.loc[:,"DMSO_R1":"BMT174819_R2"] #- Median_Frame

Probe_Targets['Name']=Log2_Frame.Name
Competed = []
Probe_Targets.fillna(0,inplace=True)
#Probe_Targets = Probe_Targets.replace({'-':''}, regex=True)
Probe_Targets.loc[:,"DMSO_R1":"BMT174819_R2"] = Probe_Targets.loc[:,"DMSO_R1":"BMT174819_R2"].apply(pd.to_numeric)
Probe_Targets = Probe_Targets.set_index('Name')
Probe_Targets.replace(np.inf, np.nan, inplace =True)
Probe_Targets.replace(-np.inf, np.nan, inplace =True)
Probe_Targets = Probe_Targets.dropna()
Probe_Targets = Probe_Targets.drop_duplicates()



corrMatrix = Probe_Targets.corr()
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
hits_pbs = Probe_Targets[Probe_Targets['BMT174819_pval']<0.005]#&(Probe_Targets['diff_525']>2.9)&(Probe_Targets['diff_526']>1)&(Probe_Targets['diff_888']>1)]
hits_pbs_b = hits_pbs_b[hits_pbs_b['BMT174819_pval']<0.005]#&(Probe_Targets['diff_525']>2.9)&(Probe_Targets['diff_526']>1)&(Probe_Targets['diff_888']>1)]

#hits_pbs = hits_pbs[~hits_pbs['Name'].str.contains('TUBB')]
Probe_Targetsb = Probe_Targets.copy()
Probe_Targetsb['DMSO_R2'] = Probe_Targets['DMSO Avg']
#Probe_Targets = Probe_Targets.sort_values(by=['BMT174819_R1','BMT174819_R2','BMT182526_R1','BMT181525_R1','BMT179888_R1','DMSO_R1'],ascending =[False,False,True,True,True,True])#ascending =[False,True,True,True,True]
s1 = pd.concat((hits_pbs, hits_pbs_b), join="inner")
output = s1.groupby(s1.index).mean().sort_values(by=['BMT174819 Avg','BMT182526 Avg','BMT181525 Avg','BMT179888 Avg','DMSO Avg'],ascending =[False,True,True,True,True]).drop_duplicates()
output = output.loc[~output.index.str.contains("TUBB|NEF")]
#color_pal = sns.color_palette(colors, n_colors=2)
#g=sns.clustermap(heatmap_set,cmap = "RdYlBu", col_cluster= True,figsize = (7,7),method='centroid') #col_cluster=False #standard_scale=1, metric = 'correlation' 
g=sns.clustermap(data=Probe_Targetsb.loc[:,"DMSO_R1":"BMT174819_R2"],cmap = "Reds",cbar_kws={'label':r'Log$_{2}$ Reporter Ion Intensity'},figsize=(6,6.5),cbar_pos=(1.3, 1.1,0.04,0.2), col_cluster= True,row_cluster=False,method='centroid',vmin=13,vmax=18,dendrogram_ratio=0.15) #g=sns.clustermap(data=Probe_Targets.loc[:,"DMSO_R1":"BMT174819_R2"],cmap = "RdYlBu",cbar_kws={'label':r'Log$_{2}$ Reporter Ion Intensity'},figsize=(5.5,6),cbar_pos=(1.3, 1.1,0.04,0.2), col_cluster= True,row_cluster=False,method = 'median') 

pos= g.ax_heatmap.get_yticks()
pos_x = g.ax_heatmap.get_xticks()
labels = g.ax_heatmap.get_yticklabels()
newpos = pos[0],pos[1]
newpos_x = pos_x[0::1]
newpos_xm = pos_x[0::2]
newpos_xm = newpos_xm+0.5
#xlabel = ['BMT181525','BMT181525','DMSO','DMSO','BMT182526','BMT182526','BMT179888','BMT179888','BMT174819','BMT174819',]
xlabel = ['DMSO','182526','181525','179888','BMT-819']
#xlabel = ['BMT181525','','DMSO','','BMT182526','','BMT179888','','BMT174819']

labelsnew = ['VDAC2','OXA1L'] + labels[2:]
plt.setp(g.ax_heatmap.set_yticks(newpos))
g.ax_heatmap.tick_params(length=12,width=1,which='major',axis='y')
plt.setp(g.ax_heatmap.set_xticks(newpos_x))
plt.setp(g.ax_heatmap.set_xticks(newpos_x,minor=True))
g.ax_heatmap.tick_params(length=6,width=1,which='minor',axis='x')
plt.setp(g.ax_heatmap.set_xticks(newpos_xm))
g.ax_heatmap.tick_params(length=10,width=2,which='major',axis='x',pad=0)
plt.setp(g.ax_heatmap.set_yticklabels(''))
#plt.setp(g.ax_heatmap.set_yticklabels(labelsnew))
for _, spine in g.ax_heatmap.spines.items():
    spine.set_visible(True)
plt.setp(g.ax_heatmap.set_xticklabels(xlabel), rotation=45,fontsize=18,fontname='Arial',fontweight='bold',ha='right')#fontweight='bold',f=sns.clustermap(Total_Frame, cmap='RdYlBu',col_cluster=False, yticklabels=Compounds,xticklabels=xlabel, method = 'centroid')
#plt.setp(g.ax_heatmap.set_yticklabels(['OXA1L','MT-CO2']),fontsize=18,fontweight='bold',fontname='Arial',)
plt.setp(g.ax_heatmap.set_ylabel(''),fontsize=18,fontweight='bold',fontname='Arial',ha='center')
plt.savefig('HeatMap_PVP_S1_replicates_032023_E.png',dpi=600,transparent=True,bbox_inches='tight',pad_inches=2)
plt.show()
Probe_Targets.replace(np.inf, np.nan, inplace =True)
Probe_Targets.replace(-np.inf, np.nan, inplace =True)
Probe_Targets = Probe_Targets.sort_values(by=['BMT174819 Avg','BMT182526 Avg','BMT181525 Avg','BMT179888 Avg','DMSO Avg'],ascending =[False,True,True,True,True])#ascending =[False,True,True,True,True]
#print(table)
Probe_Targets.replace(np.nan, np.nan, inplace=True)
Probe_Targets.dropna(inplace=True)
#Probe_Targets.to_excel('Analyzed_B819_PVP_S2.xlsx')
#output.to_excel('Combined_hitsc.xlsx')

### Make the heatmap based on the table

heatmap_set = Probe_Targets.iloc[:,1:10].copy()
heatmap_set.dropna(inplace=True)
xlabel = ['DMSO','BMT181525','BMT182526','BMT179888','BMT174819']
colors = ["skyblue", "red",'red' ]
color_pal = sns.color_palette(colors, n_colors=2)
#g=sns.clustermap(heatmap_set,cmap = "RdYlBu", col_cluster= True,figsize = (7,7),method='centroid') #col_cluster=False #standard_scale=1, metric = 'correlation' 
#g=sns.clustermap(Probe_Targets[['DMSO Avg','BMT181525 Avg','BMT182526 Avg','BMT179888 Avg','BMT174819 Avg']],cmap = "RdYlBu",cbar_kws={'label':r'Log$_{2}$ Reporter Ion Intensity'},figsize=(5.35,6),cbar_pos=(1.3, 1.1,0.04,0.2), col_cluster= True,row_cluster=False,method = 'average') #vmin=-2, vmax=6,standard_scale=1,,linewidths=0.00001,linecolor='white', method='median'col_cluster=False #standard_scale=1, metric = 'correlation' 
