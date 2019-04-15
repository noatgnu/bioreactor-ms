import pandas as pd
import numpy as np

fdr_filename = "20190204_Schulz_Luci_all_DDA_9999P6T6min - FDR.txt"
peptide_filename = "20190204_Schulz_Luci_all_DDA_9999P6T6min - Peptides.txt"
threshold = 0.01
check_decoy = True
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
            self.fdr_pass = check == (reverse >= threshold)
        return self.fdr_pass


fdr = pd.read_csv(fdr_filename, sep="\t")
peptides = pd.read_csv(peptide_filename, sep="\t")
tryp = "sp|P00761|TRYP_PIG"
cluster = {}
current = 0
for i in range(4, 85, 1):
    if (i-4)%3 == 0:
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
    for i in range(4, len(fdr.columns), 1):

        if not row["Decoy"]:
            current_row = FDRRow(row["Protein"], row["Peptide"], row["Precursor Charge"], row[fdr.columns[i]])
        else:
            if check_decoy:
                current_row.check_threshold(threshold, row[fdr.columns[i]])

        if row[fdr.columns[i]] < threshold:
            if fdr.columns[i] not in proteins:
                proteins[peptides.columns[i]] = {}

            if row["Protein"] not in proteins[peptides.columns[i]]:
                proteins[peptides.columns[i]][row["Protein"]] = 0

            if fdr.at[ind, fdr.columns[i+2]] < threshold:
                proteins[peptides.columns[i]][row["Protein"]]  += row[peptides.columns[i]]

proteins = pd.DataFrame(proteins)

for ind, row in peptides.iterrows():
