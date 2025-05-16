import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import pandas as pd
import math
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objs as go
import numpy as np
from scipy.signal import find_peaks

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
#mesSeq = SeqIO.parse("../data/Fasta/GCA_000002765.1_ASM276v1_genomic.fna","fasta")

window_size = 500
print(f"Pourcentage GC%")
for maseq in mesSeq:
    
    seq = maseq.seq.upper()
    gc_content = round((seq.count("G")+seq.count("C"))/len(seq)*100)
    at_content = round((seq.count("A")+seq.count("T"))/len(seq)*100)

    print(f"{maseq.name} GC {str(gc_content)} %  AT {str(at_content)} %")

    positions = []
    gc_skews = []
    at_skews = []
    entropies = []
    gc = []
    at = []

    # Calcul du GC skew
    for i in range(0, len(seq) - window_size + 1, window_size):
        window = seq[i:i + window_size]
        g = window.count("G")
        c = window.count("C")
        a = window.count("A")
        t = window.count("T")

        gc.append((g + c) / (a + t + g + c))
        at.append((a + t) / (a + t + g + c))


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

    positions = range(len(seq))
    y_vals = [0]
    for base in seq:
        y = y_vals[-1]
        if base == 'A':
            y += 1
        elif base == 'C':
            y += 1
        elif base == 'G':
            y -= 1
        elif base == 'T':
            y -= 1
        y_vals.append(y)
    
    # Supprimer la première valeur (0) pour aligner avec les positions
    y_vals = y_vals[1:]
    
    # Détection des pics
    # Ajuster ces paramètres selon vos besoins
    prominence = 500  # Ajuster cette valeur pour contrôler la sensibilité
    distance = 5000  # Distance minimale entre les pics en bases
    
    # Trouver les pics positifs
    peaks_pos, _ = find_peaks(y_vals, prominence=prominence, distance=distance)
    
    # Tracé
    plt.figure(figsize=(12, 4))
    plt.plot(positions, y_vals, label=f'Y skew - {maseq.name}', color='teal')
    
    # Ajouter des barres verticales pour les pics positifs
    for peak in peaks_pos:
        plt.axvline(x=peak, color='green', linestyle='--', linewidth=2)
    

    # Centromères
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        plt.axvline(x=cent_start, color='red', linestyle='--', linewidth=2, label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', linewidth=2, label='Centromere End')
    
    # Autres éléments
    plt.title(f'Y skew plot for {maseq.name} with peaks')
    plt.xlabel('Position in sequence')
    plt.ylabel('Y skew')
    plt.axhline(0, color='gray', linestyle='--')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../output/dataExploration/Y/{maseq.name}.png")
    plt.close()

    # Tracé
    plt.figure(figsize=(12, 4))
    plt.plot(positions, y_vals, label=f'Y skew - {maseq.name}', color='teal')
    
    # Ajouter des barres verticales pour les pics positifs
    for peak in peaks_pos:
        plt.axvline(x=peak, color='green', linestyle='--', linewidth=2)
    

    # Centromères
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        plt.axvline(x=cent_start, color='red', linestyle='--', linewidth=2, label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', linewidth=2, label='Centromere End')
    margin = 10000
    zoom_start = max(0, cent_start - margin)
    zoom_end = min(len(y_vals), cent_end + margin)
    plt.xlim(zoom_start, zoom_end)
    # Autres éléments
    plt.title(f'Y skew plot for {maseq.name} with peaks')
    plt.xlabel('Position in sequence')
    plt.ylabel('Y skew')
    plt.axhline(0, color='gray', linestyle='--')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../output/dataExploration/Y/{maseq.name}_zomm.png")
    plt.close()




    

    # for i in range(half_window, len(y_vals) - half_window,500):
    #     window = y_vals[i - half_window: i + half_window + 1]
    #     if y_vals[i] == max(window):
    #         peaks.append(i)
            
    # plt.figure(figsize=(12, 4))
    # x_vals = list(range(len(y_vals)))  # ou utilise une vraie position si nécessaire
    # plt.plot(x_vals, y_vals, label=f'GC - {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # #Sommets
    # for peak in peaks:
    #     plt.axvline(x=peak, color='green', linestyle='--', alpha=0.6)
    # plt.title(f'Y skew plot for {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('Y')
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/Y/{maseq.name}.png")
    # plt.close()

    # # Calcul de la dérivée discrète de y
    # y_deriv = np.diff(y_vals)
    # x_vals_deriv = x_vals[1:]  # aligné avec y_deriv

    # # Lissage de la dérivée sur une fenêtre de 150
    # window = 1000
    # kernel = np.ones(window) / window
    # smoothed_deriv = np.convolve(y_deriv, kernel, mode='valid')

    # # Centrage de l'axe x pour correspondre à la fenêtre
    # half_window = window // 2
    # x_smoothed = x_vals_deriv[window - 1:]
    # # Tracé du graphique
    # plt.figure(figsize=(12, 4))
    # plt.plot(x_smoothed, smoothed_deriv, label=f'Smoothed dY/dx ({window}bp)', color='orange')

    # # Ajout des centromères si présents
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')

    # # Mise en forme du graphique
    # plt.title(f'Smoothed Y derivative ({window}bp window) for {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('Smoothed dY/dx')
    # plt.axhline(0, color='gray', linestyle='--')
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()

    # # Sauvegarde du fichier
    # plt.savefig(f"../output/dataExploration/Y/{maseq.name}_deriv.png")
    # plt.close()
    #     # Convertir en tableau numpy pour faciliter le calcul de la dérivée
    # gc_skews_np = np.array(gc_skews)

    # # Calcul de la dérivée d(gc_skews) / d(position)
    # d_gc_skew = np.gradient(gc_skews_np) / np.gradient(positions)
    # # Tracer la dérivée de GC skew par rapport à la position
    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, d_gc_skew, label=f'Dérivée GC skew - {maseq.name}', color='darkorange')

    # # Si centromère est présent, ajouter une ligne pour marquer le début et la fin du centromère
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')

    # # Détails du graphique
    # plt.title(f'Dérivée de GC skew pour {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel("d(GC skew)/d(Position)")
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()

    # # Sauvegarder l'image
    # plt.savefig(f"../output/dataExploration/skewGC_derivative/{maseq.name}_derivative.png")
    # plt.close()


    # # Tracé
    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, gc, label=f'GC - {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'GC skew plot for {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('GC')
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/GC/{maseq.name}.png")
    # plt.close()

    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, at, label=f'AT - {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'AT plot for {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('AT')
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/AT/{maseq.name}.png")
    # plt.close()


    # # Tracé
    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, gc_skews, label=f'GC skew - {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'GC skew plot for {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('GC skew')
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/skewGC/{maseq.name}.png")
    # plt.close()

    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, at_skews, label=f'AT skew - {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'GC skew plot for {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('AT skew')
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/skewAT/{maseq.name}.png")
    # plt.close()

    # # === Entropy Plot ===
    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, entropies, label='Entropy', color='darkgreen')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'Entropy plot for {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('Shannon Entropy')
    # plt.axhline(y=2.0, color='gray', linestyle='--', label='Max entropy')
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/entropy/{maseq.name}.png")
    # plt.close()

    # # Calcul de la dérivée de l'entropie (dérivée numérique)
    # derivative = np.diff(entropies) / np.diff(positions)

    # # Affichage des résultats
    # plt.figure(figsize=(12, 4))

  

    # # Tracé de la dérivée de l'entropie
    # plt.plot(positions[:-1], derivative*100, label="Derivative of Entropy", color='darkgreen')
    # if maseq.name in centromeres:
    #         cent_start, cent_end = centromeres[maseq.name]
    #         plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #         plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title('Entropy and its Derivative with Respect to Position')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('Shannon Entropy')
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/entropyDerive/{maseq.name}.png")
    # plt.close()

    # # === Entropy Plot - Zoom centromere ===
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     zoom_margin = 10000  # tu peux ajuster cette marge

    #     # Filtrer les données dans cette zone
    #     zoom_positions = []
    #     zoom_entropies = []
    #     for pos, entropy in zip(positions, entropies):
    #         if cent_start - zoom_margin <= pos <= cent_end + zoom_margin:
    #             zoom_positions.append(pos)
    #             zoom_entropies.append(entropy)

    #     plt.figure(figsize=(12, 4))
    #     plt.plot(zoom_positions, zoom_entropies, label='Entropy (zoom)', color='darkgreen')
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    #     plt.title(f'Entropy plot (zoom) for {maseq.name}')
    #     plt.xlabel('Position in sequence')
    #     plt.ylabel('Shannon Entropy')
    #     plt.axhline(y=2.0, color='gray', linestyle='--')
    #     plt.grid(True)
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.savefig(f"../output/dataExploration/entropy_zoom/{maseq.name}.png")
    #     plt.close()

    # # === GC Skew Plot - Zoom centromere ===
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     zoom_margin = 10000  # ajustable

    #     # Filtrer les données autour du centromère
    #     zoom_positions = []
    #     zoom_gc_skews = []
    #     for pos, skew in zip(positions, gc_skews):
    #         if cent_start - zoom_margin <= pos <= cent_end + zoom_margin:
    #             zoom_positions.append(pos)
    #             zoom_gc_skews.append(skew)

    #     plt.figure(figsize=(12, 4))
    #     plt.plot(zoom_positions, zoom_gc_skews, label='GC skew (zoom)', color='teal')
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    #     plt.title(f'GC skew plot (zoom) for {maseq.name}')
    #     plt.xlabel('Position in sequence')
    #     plt.ylabel('GC skew')
    #     plt.axhline(0, color='gray', linestyle='--')
    #     plt.grid(True)
    #     plt.legend()
    #     plt.tight_layout()
    #     plt.savefig(f"../output/dataExploration/skewGC_zoom/{maseq.name}.png")
    #     plt.close()

    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     length_Centromere = cent_end-cent_start
    #     print(f"Taille du centromere = {length_Centromere} et le pourcentage sur du genome {length_Centromere/len(maseq.seq)*100}")

    # === Z-curve Plot ===
    # x_vals = [0]
    # y_vals = [0]
    # z_vals = [0]

    # for base in seq:
    #     x = x_vals[-1]
    #     y = y_vals[-1]
    #     z = z_vals[-1]
    #     if base == 'A':
    #         x += 1
    #         y += 1
    #         z += 1
    #     elif base == 'C':
    #         x -= 1
    #         y += 1
    #         z -= 1
    #     elif base == 'G':
    #         x += 1
    #         y -= 1
    #         z -= 1
    #     elif base == 'T':
    #         x -= 1
    #         y -= 1
    #         z += 1
    #     x_vals.append(x)
    #     y_vals.append(y)
    #     z_vals.append(z)


    # fig = plt.figure(figsize=(10, 8))
    # ax = fig.add_subplot(111, projection='3d')

    # # Centromère présent ?
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]

    #     # Tracé avant centromère
    #     ax.plot(
    #         x_vals[:cent_start],
    #         y_vals[:cent_start],
    #         z_vals[:cent_start],
    #         color='darkblue',
    #         label='Outside Centromere'
    #     )

    #     # Tracé centromère (en rouge)
    #     ax.plot(
    #         x_vals[cent_start:cent_end],
    #         y_vals[cent_start:cent_end],
    #         z_vals[cent_start:cent_end],
    #         color='red',
    #         label='Centromere'
    #     )

    #     # Tracé après centromère
    #     ax.plot(
    #         x_vals[cent_end:],
    #         y_vals[cent_end:],
    #         z_vals[cent_end:],
    #         color='darkblue'
    #     )
    # else:
    #     # Pas de centromère : tracer tout en bleu
    #     ax.plot(x_vals, y_vals, z_vals, color='darkblue', label='Z-curve')
    # ax.view_init(elev=0, azim=0) # Y et z
    # #ax.view_init(elev=90, azim=-90) #Y et X
    # ax.set_title(f'3D Z-curve for {maseq.name}')
    # ax.set_xlabel('X (Purine - Pyrimidine)')
    # ax.set_ylabel('Y (Amino - Keto)')
    # ax.set_zlabel('Z (Weak - Strong)')
    # ax.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/zcurve/A{maseq.name}.png")
    # plt.close()


# trace = go.Scatter3d(
#     x=x_vals,
#     y=y_vals,
#     z=z_vals,
#     mode='lines',
#     line=dict(color='blue', width=2),
#     name='Z-curve'
# )

# # Segment centromère en rouge
# if maseq.name in centromeres:
#     cent_start, cent_end = centromeres[maseq.name]
#     trace_centromere = go.Scatter3d(
#         x=x_vals[cent_start:cent_end],
#         y=y_vals[cent_start:cent_end],
#         z=z_vals[cent_start:cent_end],
#         mode='lines',
#         line=dict(color='red', width=4),
#         name='Centromere'
#     )
#     data = [trace, trace_centromere]
# else:
#     data = [trace]

# layout = go.Layout(
#     title=f"Interactive Z-curve for {maseq.name}",
#     scene=dict(
#         xaxis_title='X',
#         yaxis_title='Y',
#         zaxis_title='Z'
#     )
# )

# fig = go.Figure(data=data, layout=layout)
# fig.write_html(f"../output/dataExploration/zcurve_3D/{maseq.name}.html")