from sequana import FastA, Repeats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Charger le fichier FASTA
fasta = FastA("../data/Fasta/TriTrypDB-68_LmajorFriedlin_Genome.fasta")

df = pd.read_csv("../data/Centromere_Positions/centromeres.csv")
centromeres = {
    row['Chromosome']: (row['Centromere_Start'], row['Centromere_End'])
    for _, row in df.iterrows()
}
# Pour chaque séquence dans le fichier
for maseq in fasta:
    print(f"Analyse de : {maseq.name}")

    # Calculer les répétitions
    r = Repeats("../data/Fasta/TriTrypDB-68_LmajorFriedlin_Genome.fasta", name=maseq.name)
    r.threshold = 10  # Seuil minimum pour une répétition

    # Fenêtre glissante
    step = 1000
    N = len(maseq.sequence)
    positions = []
    max_repeats = []

    for i in np.arange(step / 2, N, step):
        start = int(i - step / 2)
        stop = int(i + step / 2)

        subset = r.df_shustring.query("position >= @start and position <= @stop")
        max_len = subset.shustring_length.max() if not subset.empty else 0

        positions.append(i)
        max_repeats.append(max_len)

    # Affichage du graphique
    plt.figure(figsize=(12, 5))
    plt.plot(positions, max_repeats, label=f"{maseq.name}")
    plt.xlabel("Position dans le génome")
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]

        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    plt.ylabel("Longueur max de répétition")
    plt.title(f"Profil des répétitions dans {maseq.name}")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"../output/dataExploration/repeat/{maseq.name}.png")
