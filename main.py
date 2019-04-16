import pandas as pd
from scipy.stats import sem
import numpy as np

fdr_filename = "20190204_Schulz_Luci_all_DDA_9999P6T6min - FDR.txt"
peptide_filename = "20190204_Schulz_Luci_all_DDA_9999P6T6min - Peptides.txt"
threshold = 0.01
check_decoy = False

class FDRRow:
    def __init__(self, protein, peptide, charge, fdr):
        self.protein = protein
        self.peptide = peptide
        self.charge = charge
        self.fdr = fdr
        self.fdr_pass = False

    def check_threshold(self, threshold, reverse=None):
        check = self.fdr < threshold
        if not reverse:
            self.fdr_pass = check
        else:
            self.fdr_pass = check == True and reverse >= threshold
        return self.fdr_pass


fdr = pd.read_csv(fdr_filename, sep="\t")
peptides = pd.read_csv(peptide_filename, sep="\t", index_col=[0,1,3])
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

current_row = ""
for ind, row in fdr.iterrows():
    for i in range(7, len(fdr.columns), 1):

        if not row["Decoy"]:
            current_row = FDRRow(row["Protein"], row["Peptide"], row["Precursor Charge"], row[fdr.columns[i]])
            if not check_decoy:
                if current_row.check_threshold(threshold):
                    if peptides.columns[i-5] not in proteins:
                        proteins[peptides.columns[i-5]] = {}

                    if row["Protein"] not in proteins[peptides.columns[i-5]]:
                        proteins[peptides.columns[i-5]][row["Protein"]] = 0

                    proteins[peptides.columns[i-5]][row["Protein"]] += peptides.at[(current_row.protein, current_row.peptide, current_row.charge), peptides.columns[i-5]]
        else:
            if check_decoy:
                if current_row.check_threshold(threshold, row[fdr.columns[i]]):
                    if peptides.columns[i-5] not in proteins:
                        proteins[peptides.columns[i-5]] = {}

                    if row["Protein"] not in proteins[peptides.columns[i-5]]:
                        proteins[peptides.columns[i-5]][row["Protein"]] = 0

                    proteins[peptides.columns[i-5]][row["Protein"]] += peptides.at[(current_row.protein, current_row.peptide, current_row.charge), peptides.columns[i-5]]

proteins = pd.DataFrame(proteins)

trypsin = proteins.loc[tryp]

for ind, row in proteins.iterrows():
    stats = {}
    for i in range(len(proteins.columns)):
        if pd.notnull(trypsin[proteins.columns[i]]) and pd.notnull(proteins.at[ind, proteins.columns[i]]):
            proteins.at[ind, proteins.columns[i]] = proteins.at[ind, proteins.columns[i]]/trypsin[proteins.columns[i]]
            if cluster[i+7] not in stats:
                stats[cluster[i+7]] = {"name": "", "values": []}
            stats[cluster[i+7]]["name"] += proteins.columns[i]
            stats[cluster[i+7]]["values"].append(proteins.at[ind, proteins.columns[i]])

    for i in stats:
        se = sem(stats[i]["values"], nan_policy="omit")
        ave = np.average(stats[i]["values"])



