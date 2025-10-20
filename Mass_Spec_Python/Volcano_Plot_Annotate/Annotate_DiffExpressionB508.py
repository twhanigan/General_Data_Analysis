import pandas as pd
from bioinfokit import visuz
import numpy as np
import pandas as pd
import math
import numpy as d
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
from sklearn.cluster import KMeans
from more_itertools import unique_everseen
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
from matplotlib.animation import FFMpegWriter
from matplotlib import style
import statistics
import mygene

df= pd.read_excel('H460_ExperimentalvsControl_DEResults_TMMf.xlsx')
pep = pd.read_excel('proList_2pep.xlsx')
peps = pep[['description','peptide num','accession']]
peps = peps.loc[~(peps['accession'].str.contains('Reverse'))]

peps['description'] = peps['description'].str.split(" ", 1).str.get(0)
print(peps)
Mutations = df.merge(peps,how='left', on='description')
print(Mutations.head())

mg = mygene.MyGeneInfo()
gene_name_list = Mutations['accession']#.str.split("] ", 1).str.get(1).str.split(';',1).str.get(0)
Mutations.index = gene_name_list

protein_name = mg.querymany(gene_name_list, scopes='uniprot', fields=[
                            'symbol','ensembl.gene','genomic_pos.chr','genomic_pos','summary','exons','kegg','go','pathway','reactome','type_of_gene'], species = 'human', as_dataframe = True, df_index = True)
protein_name = protein_name.drop_duplicates(subset='symbol')
merged_mutations = Mutations.merge(
    protein_name, how='left', left_index=True, right_index=True)
merged_mutations.to_csv('H460_Protein_TMM_Annotated_102.csv')