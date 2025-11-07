import pandas as pd
import mygene

# Load the data
df = pd.read_csv('abundance_protein_MD.tsv', sep='\t')
print("Columns in input file:", df.columns)

# Initialize MyGeneInfo
mg = mygene.MyGeneInfo()

# Clean the gene name list
gene_name_list = df['Gene'].dropna().unique().tolist()  # Drop NaNs and ensure uniqueness

# Query MyGene.info
protein_name = mg.querymany(
    gene_name_list,
    scopes='symbol',
    fields=[
        'uniprot', 'ensembl.gene', 'genomic_pos.chr', 'genomic_pos',
        'summary', 'exons', 'kegg', 'go', 'pathway', 'reactome', 'type_of_gene'
    ],
    species='human',
    as_dataframe=True,
    df_index=False
)

print(protein_name.head())

# Save MyGene results
protein_name.to_csv('test.csv')

# Merge with original dataframe
merged_df = df.merge(
    protein_name, how='left', left_on='Gene', right_on='query'
)

# Save annotated dataframe
merged_df.to_csv('abundance_protein_MD_Annotated.csv', index=False)