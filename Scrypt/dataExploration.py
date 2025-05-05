import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

window_size = 1000
mesSeq = SeqIO.parse("../Fasta/TriTrypDB-68_LmajorFriedlin_Genome.fasta","fasta")
print(f"Pourcentage GC%")
t = 1
for maseq in mesSeq:
    gc_content = round((maseq.seq.count("G")+maseq.seq.count("C"))/len(maseq.seq)*100)
    print(f"{maseq.name}  {str(gc_content)} %")

    seq = maseq.seq
    positions = []
    gc_skews = []

    # Calcul du GC skew
    for i in range(0, len(seq) - window_size + 1, window_size):
        window = seq[i:i + window_size]
        g = window.count("G")
        c = window.count("C")
        if g + c == 0:
            skew = 0
        else:
            skew = (g - c) / (g + c)

        positions.append(i + window_size // 2)  
        gc_skews.append(skew)

    # Tracé
    plt.figure(figsize=(12, 4))
    plt.plot(positions, gc_skews, label=f'GC skew - {maseq.name}', color='teal')
    plt.title(f'GC skew plot for {maseq.name}')
    plt.xlabel('Position in sequence')
    plt.ylabel('GC skew')
    plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../output/dataExploration/skewGC/{t}-{maseq.name}.png")
    t+=1

