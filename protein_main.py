import pandas as pd
from scipy.stats import sem
import numpy as np

protein_filename = "20190204_Schulz_Luci_all_DDA_9999P6T6min - Proteins.txt"
normalized_filename = "norm.20190204_Schulz_Luci_all_DDA_9999P6T6min - Proteins.txt"
average_filename = "ave.20190204_Schulz_Luci_all_DDA_9999P6T6min - Proteins.txt"
sem_filename = "sem.20190204_Schulz_Luci_all_DDA_9999P6T6min - Proteins.txt"

proteins = pd.read_csv(protein_filename, sep="\t", index_col=0)


cluster = {}
current = 0
for i in range(0, 78, 1):
    if i%3 == 0:
        current += 1
    cluster[i] = current

for i in range(78, 96, 1):
    if i not in cluster:
        current += 1
        for i2 in range(i, i+6, 2):
            cluster[i2] = current

tryp = "sp|P00761|TRYP_PIG"

trypsin = proteins.loc[tryp].copy(deep=True)
data = []

for ind, row in proteins.iterrows():
    stats = {}
    for i in range(0, len(proteins.columns), 1):
        # if pd.notnull(trypsin[proteins.columns[i]]) and pd.notnull(proteins.at[ind, proteins.columns[i]]):
            proteins.at[ind, proteins.columns[i]] = proteins.at[ind, proteins.columns[i]]/trypsin[proteins.columns[i]]
            if cluster[i] not in stats:
                stats[cluster[i]] = {"name": [], "values": [], "nona": 0}
            stats[cluster[i]]["name"].append(proteins.columns[i])
            stats[cluster[i]]["values"].append(proteins.at[ind, proteins.columns[i]])
            if pd.notnull(proteins.at[ind, proteins.columns[i]]):
                stats[cluster[i]]["nona"] += 1

    for i in stats:
        se = sem(stats[i]["values"], nan_policy="omit")
        ave = np.average(stats[i]["values"])
        data.append([ind, "+".join(stats[i]["name"]), ave, se, stats[i]["nona"]])

data = pd.DataFrame(data, columns=["Protein", "Sample", "Average", "SEM", "Non_Nan_Count"])
proteins.reset_index().to_csv(normalized_filename, sep="\t")
data[["Protein", "Sample", "Average"]].pivot(index="Protein", columns="Sample").reset_index().to_csv(average_filename, sep="\t", index=False)
data[["Protein", "Sample", "SEM"]].pivot(index="Protein", columns="Sample").reset_index().to_csv(sem_filename, sep="\t", index=False)
