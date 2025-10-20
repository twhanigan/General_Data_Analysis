# Plot Analyzed CTG
from itertools import cycle
import fnmatch
import glob
import sys
import os
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

frame = pd.read_csv('Combined_Frame_10252025.csv')
evaluated = pd.read_csv('Combined_Evaluated_Frame_10252025.csv')
print(evaluated.head())
colors = ["steelblue", "firebrick", "red",'k','grey']
cmap = matplotlib.cm.get_cmap('RdYlBu')
#color_pal = sns.color_palette('RdYlBu', n_colors=3)
color_pal = sns.color_palette(colors, n_colors=5)
hex_colors = list(color_pal.as_hex())
#subset = frame.loc[((frame['frame_name'] == 'H460_NH4Cl_Final.CSV') & (frame['compound'] == 'Leflun_20')) | ((frame['frame_name'] == 'PC3_NH4Cl_Final.CSV') & (frame['compound'] == 'Leflun_20'))  | (
#    (frame['frame_name'] == 'PC3_Leflunamide.CSV') & (frame['compound'] == 'Leflun_10')) | ((frame['frame_name'] == 'PC3_NH4Cl_Final.CSV') & (frame['compound'] == 'Leflun_10')) | ((frame['frame_name'] == 'PC3_NH4Cl_Final.CSV') & (frame['compound'] == 'Leflun_5'))].copy().reset_index(drop=True)
#print(subset)
compoundData = frame.groupby(['compound', 'frame_name'], sort=False)
#compoundData = frame.groupby(['compound','frame_name'],sort=False)
nums = [0, 1, 2]
licycle = cycle(nums)

fig, ax = plt.subplots(figsize=(4, 3.5))
for name, group in compoundData:
    h = next(licycle)
    if name[1] == 'TRno1.CSV':
        #ax.scatter(x='Log Concentration',y='Normalized',data=group,color=hex_colors[h],cmap='RdYlBu',s=1)
        ax.plot(evaluated['dose'].loc[(evaluated['compound'] == name[0]) & (evaluated['frame_name'] == name[1])], evaluated['new'].loc[(
            evaluated['compound'] == name[0]) & (evaluated['frame_name'] == name[1])], linewidth=2, alpha=0.75, color=hex_colors[h])
        ax.scatter(frame['Log Concentration'].loc[(frame['compound'] == name[0]) & (frame['frame_name'] == name[1])], frame['Normalized'].loc[(
            frame['compound'] == name[0]) & (frame['frame_name'] == name[1])], linewidth=2, alpha=0.75, color=hex_colors[h])
        # sns.lineplot(x='dose',y='new',data=evaluated_frame,hue='compound',palette=color_pal)
        # print(evaluated['dose'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])])
        plt.fill_between(evaluated['dose'].loc[(evaluated['compound'] == name[0]) & (evaluated['frame_name'] == name[1])], evaluated['lowerbound'].loc[(evaluated['compound'] == name[0]) & (
            evaluated['frame_name'] == name[1])], evaluated['upperbound'].loc[(evaluated['compound'] == name[0]) & (evaluated['frame_name'] == name[1])], alpha=0.01, color=hex_colors[h], cmap='RdYlBu')
    else:
        #ax.scatter(x='Log Concentration',y='Normalized',data=group,color=hex_colors[h],cmap='RdYlBu',s=1)
        ax.plot(evaluated['dose'].loc[(evaluated['compound'] == name[0]) & (evaluated['frame_name'] == name[1])], evaluated['new'].loc[(
            evaluated['compound'] == name[0]) & (evaluated['frame_name'] == name[1])], linewidth=2, alpha=0.75, color=hex_colors[h])
        ax.scatter(frame['Log Concentration'].loc[(frame['compound'] == name[0]) & (frame['frame_name'] == name[1])], frame['Normalized'].loc[(
            frame['compound'] == name[0]) & (frame['frame_name'] == name[1])], linewidth=2, alpha=0.75, color=hex_colors[h])
        # sns.lineplot(x='dose',y='new',data=evaluated_frame,hue='compound',palette=color_pal)
        # print(evaluated['dose'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])])
        plt.fill_between(evaluated['dose'].loc[(evaluated['compound'] == name[0]) & (evaluated['frame_name'] == name[1])], evaluated['lowerbound'].loc[(evaluated['compound'] == name[0]) & (
            evaluated['frame_name'] == name[1])], evaluated['upperbound'].loc[(evaluated['compound'] == name[0]) & (evaluated['frame_name'] == name[1])], alpha=0.01, color=hex_colors[h], cmap='RdYlBu')

handles, labels = plt.gca().get_legend_handles_labels()
print(handles)
#handle3 = handles[3].set_color('red')
#handle4 = handles[4].set_color('red')
#handle5 = handles[5].set_color('red')
#handle6 = handles[6].set_color('red')
#handle7 = handles[7].set_color('red')
#lst_handle = handles[8].set_color('red')

#handles_new = [handles[0], handles[1], handles[2], handles[3],
#               handles[5], handles[4]]



labels_new =labels
ax.legend(labels_new, title='',
           bbox_to_anchor=(1.33, 1.12), fontsize=10)

plt.xlabel('Log[I], M', fontsize=18, fontname="Arial", fontweight='bold')
plt.xticks(fontsize=12)
plt.ylabel('Relative Viability', fontsize=18,
           fontname="Arial", fontweight='bold')
plt.yticks(fontsize=12)
plt.ylim(0, 2)
plt.xlim(-9.5, -5)

plt.savefig('CTG_Output_10252025.png', dpi=600,bbox_inches ='tight',pad_inches=2,transparent=True)
plt.show()
