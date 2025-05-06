import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import pandas as pd
import math

def calculate_entropy(seq_window):
    length = len(seq_window)
    if length == 0:
        return 0.0

    freq = {}
    for base in 'ACGT':
        freq[base] = seq_window.count(base) / length

    entropy = 0.0
    for p in freq.values():
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy


df = pd.read_csv("../data/Centromere_Positions/centromeres.csv")
centromeres = {
    row['Chromosome']: (row['Centromere_Start'], row['Centromere_End'])
    for _, row in df.iterrows()
}


mesSeq = SeqIO.parse("../data/Fasta/TriTrypDB-68_LmajorFriedlin_Genome.fasta","fasta")
window_size = 1000
print(f"Pourcentage GC%")
for maseq in mesSeq:
    gc_content = round((maseq.seq.count("G")+maseq.seq.count("C"))/len(maseq.seq)*100)
    at_content = round((maseq.seq.count("A")+maseq.seq.count("T"))/len(maseq.seq)*100)

    print(f"{maseq.name} GC {str(gc_content)} %  AT {str(at_content)} %")

    seq = maseq.seq
    positions = []
    gc_skews = []
    at_skews = []
    entropies = []

    # Calcul du GC skew
    for i in range(0, len(seq) - window_size + 1, window_size):
        window = seq[i:i + window_size]
        g = window.count("G")
        c = window.count("C")
        a = window.count("A")
        t = window.count("T")


        if g + c == 0:
            skew = 0
        else:
            skew = (g - c) / (g + c)

        gc_skews.append(skew)

        if  a + t == 0:
            skew = 0
        else:
            skew = (a - t) / (a + t)
        at_skews.append(skew)

        entropy = calculate_entropy(window)
        entropies.append(entropy)
        positions.append(i + window_size // 2)  



    # Tracé
    plt.figure(figsize=(12, 4))
    plt.plot(positions, gc_skews, label=f'GC skew - {maseq.name}', color='teal')
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    plt.title(f'GC skew plot for {maseq.name}')
    plt.xlabel('Position in sequence')
    plt.ylabel('GC skew')
    plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../output/dataExploration/skewGC/{maseq.name}.png")

    plt.figure(figsize=(12, 4))
    plt.plot(positions, at_skews, label=f'AT skew - {maseq.name}', color='teal')
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    plt.title(f'GC skew plot for {maseq.name}')
    plt.xlabel('Position in sequence')
    plt.ylabel('GC skew')
    plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../output/dataExploration/skewAT/{maseq.name}.png")


    # === Entropy Plot ===
    plt.figure(figsize=(12, 4))
    plt.plot(positions, entropies, label='Entropy', color='darkgreen')
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    plt.title(f'Entropy plot for {maseq.name}')
    plt.xlabel('Position in sequence')
    plt.ylabel('Shannon Entropy')
    plt.axhline(y=2.0, color='gray', linestyle='--', label='Max entropy')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../output/dataExploration/entropy/{maseq.name}.png")
    plt.close()

    # === Entropy Plot - Zoom centromere ===
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        zoom_margin = 10000  # tu peux ajuster cette marge

        # Filtrer les données dans cette zone
        zoom_positions = []
        zoom_entropies = []
        for pos, entropy in zip(positions, entropies):
            if cent_start - zoom_margin <= pos <= cent_end + zoom_margin:
                zoom_positions.append(pos)
                zoom_entropies.append(entropy)

        plt.figure(figsize=(12, 4))
        plt.plot(zoom_positions, zoom_entropies, label='Entropy (zoom)', color='darkgreen')
        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
        plt.title(f'Entropy plot (zoom) for {maseq.name}')
        plt.xlabel('Position in sequence')
        plt.ylabel('Shannon Entropy')
        plt.axhline(y=2.0, color='gray', linestyle='--')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"../output/dataExploration/entropy_zoom/{maseq.name}.png")
        plt.close()

    # === GC Skew Plot - Zoom centromere ===
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        zoom_margin = 10000  # ajustable

        # Filtrer les données autour du centromère
        zoom_positions = []
        zoom_gc_skews = []
        for pos, skew in zip(positions, gc_skews):
            if cent_start - zoom_margin <= pos <= cent_end + zoom_margin:
                zoom_positions.append(pos)
                zoom_gc_skews.append(skew)

        plt.figure(figsize=(12, 4))
        plt.plot(zoom_positions, zoom_gc_skews, label='GC skew (zoom)', color='teal')
        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
        plt.title(f'GC skew plot (zoom) for {maseq.name}')
        plt.xlabel('Position in sequence')
        plt.ylabel('GC skew')
        plt.axhline(0, color='gray', linestyle='--')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"../output/dataExploration/skewGC_zoom/{maseq.name}.png")
        plt.close()
