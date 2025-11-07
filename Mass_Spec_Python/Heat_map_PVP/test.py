import numpy as np
import pandas as pd
import math 
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
#from sklearn.metrics import jaccard_similarity_score
from more_itertools import unique_everseen
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
Probe_Targets['diff_DMSO'] = Probe_Targets['BMT174819 Avg']-Probe_Targets['DMSO Avg']
Probe_Targets['diff_525'] = Probe_Targets['BMT174819 Avg']-Probe_Targets['BMT181525 Avg']
Probe_Targets['diff_526'] = Probe_Targets['BMT174819 Avg']-Probe_Targets['BMT182526 Avg']
Probe_Targets['diff_888'] = Probe_Targets['BMT174819 Avg']-Probe_Targets['BMT179888 Avg']

hits_pbs = Probe_Targets[Probe_Targets['DMSO Avg']>0]
print(hits_pbs)
