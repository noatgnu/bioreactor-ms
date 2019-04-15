import pandas as pd

fdr_filename = ""
protein_filename = ""

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

for ind, row in peptides.iterrows():
    for i in range(7, len(peptides.columns), 1):


for ind, row in peptides.iterrows():
