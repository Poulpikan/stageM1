from sequana import FastA
import numpy as np
import matplotlib.pyplot as plt
import re
import pandas as pd

# Définir une fonction pour compter les homopolymères
def count_homopolymers(seq, min_length=5):
    # Ex : trouve AAAAA ou TTTTT, etc.
    pattern = re.compile(rf"(A{{{min_length},}}|T{{{min_length},}}|C{{{min_length},}}|G{{{min_length},}})")
    return len(pattern.findall(seq.upper()))

# Paramètres
fasta = FastA("../data/Fasta/TriTrypDB-68_LmajorFriedlin_Genome.fasta")
window_size = 500
min_homopolymer_length = 5
df = pd.read_csv("../data/Centromere_Positions/centromeres.csv")
centromeres = {
    row['Chromosome']: (row['Centromere_Start'], row['Centromere_End'])
    for _, row in df.iterrows()
}
# Parcourir les séquences
for maseq in fasta:
    sequence = maseq.sequence.upper()
    N = len(sequence)
    x = []
    y = []

    for i in range(0, N, window_size):
        window = sequence[i:i + window_size]
        nb = count_homopolymers(window, min_length=min_homopolymer_length)
        x.append(i)
        y.append(nb)

    # Affichage du graphique
    plt.figure(figsize=(12, 5))
    plt.plot(x, y, label=f"{maseq.name}")
    plt.title(f"Nombre de homopolymères ≥ {min_homopolymer_length} nt dans {maseq.name}")
    plt.xlabel("Position dans la séquence (bp)")
    plt.ylabel("Nombre de homopolymères")
    plt.grid(True)
    plt.legend()
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]

        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    plt.savefig(f"../output/dataExploration/homopolymer/{maseq.name}.png")
