import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import pandas as pd
import math


# Lire le fichier CSV
df = pd.read_csv("../data/Centromere_Positions/centromeres.csv")

# Créer un dictionnaire avec le nom du chromosome et ses coordonnées de centromère
centromeres = {
    row['Chromosome']: (row['Centromere_Start'], row['Centromere_End'])
    for _, row in df.iterrows()
}

# Afficher les positions du centromère
for chromosome, (cent_start, cent_end) in centromeres.items():
    print(f"{chromosome}: taille = {cent_end-cent_start}")
