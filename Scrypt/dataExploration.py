import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import pandas as pd

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

