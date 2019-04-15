import pandas as pd
import numpy as np

fdr_filename = ""
protein_filename = ""
threshold = 0.01


fdr = pd.read_csv("")
peptides = pd.read_csv("")
tryp = "sp|P00761|TRYP_PIG"
cluster = {}
current = 0
for i in range(7, 85, 1):
    if (i-7)%3 == 0:
        current += 1
    cluster[i] = current

for i in range(85, 103, 1):
    if i not in cluster:
        current += 1
        for i2 in range(i, i+6, 2):
            cluster[i2] = current

proteins = {}

for ind, row in peptides.iterrows():
    for i in range(7, len(peptides.columns), 1):
        if peptides["Protein"] not in proteins:
            proteins[peptides["Protein"]] = {}
        if peptides.columns[i] not in proteins[peptides["Protein"]]:
            proteins[peptides["Protein"]] = 0
        if row[peptides.columns[i]] >= threshold:
            peptides.at[ind, peptides.columns[i]] = np.nan
        else:
            proteins[peptides["Protein"]] += peptides.at[ind]

for ind, row in peptides.iterrows():
