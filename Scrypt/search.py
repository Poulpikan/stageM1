from sequana import DNA
from sequana import FastA
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pywt
from scipy.signal import find_peaks
import re


f = FastA("../data/Fasta/GCA_000002765.1_ASM276v1_genomic.fna")
df = pd.read_csv("../data/Centromere_Positions/centromeres.csv")
centromeres = {
    row['Chromosome']: (row['Centromere_Start'], row['Centromere_End'])
    for _, row in df.iterrows()
}

for maseq in f:
    #Etape 5
    s = DNA(maseq.sequence)
    min_homopolymer_length = 5
    s.window = 100

    karlin = s.get_karlin_signature_difference()
    karlin = karlin[5000:-5000]
    positions = range(len(karlin))


    df_karlin = pd.DataFrame({'Karlin': karlin})
    df_karlin.to_csv(f"../data/karlin/{maseq.name}_karlin_100.csv", index=False)

    mean_karlin = sum(karlin) / len(karlin)
    min_karlin = min(karlin)
    max_karlin = max(karlin)
    print(f"{maseq.name} moy = {mean_karlin:.3f}  min = {min_karlin:.3f}  max = {max_karlin:.3f}")
    entropie = s.get_entropy(500)
    # Fonction pour calculer la moyenne autour d'une position
    def average_around(position, window_size=500):
        start = max(0, position - window_size)
        end = min(len(karlin), position + window_size)
        return sum(karlin[start:end]) / (end - start)

    def entropy_max_around(position,entropy, window_size=1000):
        start = max(0, position - window_size)
        end = min(len(entropy), position + window_size)
        return max(entropy[start:end])
    def count_homopolymers(seq, min_length=5):
        # Ex : trouve AAAAA ou TTTTT, etc.
        pattern = re.compile(rf"(A{{{min_length},}}|T{{{min_length},}}|C{{{min_length},}}|G{{{min_length},}})")
        return len(pattern.findall(seq.upper()))
    window_size = 500
    N = len(maseq.sequence)
    x_homopolymere = []
    y_homopolymere  = []

    for i in range(5000, N-5000, window_size):
        window = maseq.sequence[i:i + window_size]
        nb = count_homopolymers(window, min_length=min_homopolymer_length)
        x_homopolymere.append(i)
        y_homopolymere.append(nb)
    # Tracé
    plt.figure(figsize=(12, 4))
    plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    
    y_vals = [0]
    for base in maseq.sequence:
        y = y_vals[-1]
        if base == 'A' or base == 'C':
            y += 1
        elif base == 'G' or base == 'T':
            y -= 1
        y_vals.append(y)

    if len(y_vals) > 10000:  # Assure-toi qu'il y a assez d'éléments
        y_vals = y_vals[5000:-5000]

    prominence = 500  # Ajuster cette valeur pour contrôler la sensibilité
    distance = 5000  # Distance minimale entre les pics en bases
    
    # Trouver les pics positifs
    peaks_pos, _ = find_peaks(y_vals, prominence=prominence, distance=distance)
    best = []
    # Ajouter des lignes verticales aux pics
    for i in range(1, len(karlin) - 1,100):  # Ne pas inclure les bords
        # Vérification si c'est un pic (2 fois supérieur à la moyenne générale)
        if karlin[i] > 1.75 * mean_karlin:
            # Vérification de la condition de la moyenne autour de 2000
            avg_around = average_around(i)
            if avg_around > mean_karlin and  entropy_max_around(i,entropie) > 1.97:
                # Vérifie si ce pic est proche d'un pic Y_skew
                is_near_y_peak = False
                is_near_y_peak = any(abs(i - yp) <= 1500 for yp in peaks_pos)
                if is_near_y_peak:
                    best.append(i)
                    # Affiche la ligne et imprime l’info
    final = []
    max_avg = -1  # pour suivre la meilleure moyenne

    for i in best:
        # Convertir les positions en indices dans y_homopolymere
        start = max(0, (i - 1000) // window_size)
        end = min(len(y_homopolymere), (i + 1000) // window_size)

        if start < end:
            avg_homopolymers = np.mean(y_homopolymere[start:end])

            if avg_homopolymers > max_avg:
                max_avg = avg_homopolymers
                final = [(i, avg_homopolymers)]
            elif avg_homopolymers == max_avg:
                final.append((i, avg_homopolymers))
    # Ajouter des lignes verticales aux positions des meilleures moyennes
    for i, _ in final:
        plt.axvline(x=i, color='orange', linestyle='--', label='Meilleure moyenne')

    # Centromères
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        cent_start = cent_start-5000
        cent_end = cent_end -5000
        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    
    plt.title(f'Karlin {maseq.name}')
    plt.xlabel('Position in sequence')
    plt.ylabel('Karlin')
    plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../output/dataExploration/search/5/{maseq.name}.png")
    plt.close()

  # Tracé Zoom
    plt.figure(figsize=(12, 4))
    plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    best = []

    # Ajouter des lignes verticales aux pics
    for i in range(1, len(karlin) - 1,100):  # Ne pas inclure les bords
        # Vérification si c'est un pic (2 fois supérieur à la moyenne générale)
        if karlin[i] > 1.75 * mean_karlin:
            # Vérification de la condition de la moyenne autour de 2000
            avg_around = average_around(i)
            if avg_around > mean_karlin and  entropy_max_around(i,entropie) > 1.97:
                # Vérifie si ce pic est proche d'un pic Y_skew
                is_near_y_peak = False
                is_near_y_peak = any(abs(i - yp) <= 1500 for yp in peaks_pos)
                if is_near_y_peak:
                    best.append(i)

    final = []
    max_avg = -1  # pour suivre la meilleure moyenne

    for i in best:
        # Convertir les positions en indices dans y_homopolymere
        start = max(0, (i - 1000) // window_size)
        end = min(len(y_homopolymere), (i + 1000) // window_size)

        if start < end:
            avg_homopolymers = np.mean(y_homopolymere[start:end])

            if avg_homopolymers > max_avg:
                max_avg = avg_homopolymers
                final = [(i, avg_homopolymers)]
            elif avg_homopolymers == max_avg:
                final.append((i, avg_homopolymers))
    # Ajouter des lignes verticales aux positions des meilleures moyennes
    for i, _ in final:
        plt.axvline(x=i, color='orange', linestyle='--', label='Meilleure moyenne')

    # Centromères
    if maseq.name in centromeres:
        cent_start, cent_end = centromeres[maseq.name]
        cent_start = cent_start-5000
        cent_end = cent_end -5000
        plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
        plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
                # Zoom autour du centromère
        margin = 10000  # marge de 10kb autour du centromère
        zoom_start = max(0, cent_start - margin)
        zoom_end = min(len(karlin), cent_end + margin)
        plt.xlim(zoom_start, zoom_end)
    
    plt.title(f'Karlin {maseq.name}')
    plt.xlabel('Position in sequence')
    plt.ylabel('Karlin')
    plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../output/dataExploration/search/5/{maseq.name}_.png")
    plt.close()

#     #Etape 4
#     s = DNA(maseq.sequence)
    
#     s.window = 100
#     df_karlin = pd.read_csv(f"../data/karlin/{maseq.name}_karlin_100.csv")
#     karlin = df_karlin['Karlin'].tolist()
#     karlin = karlin[5000:-5000]
#     positions = range(len(karlin))

#     mean_karlin = sum(karlin) / len(karlin)
#     min_karlin = min(karlin)
#     max_karlin = max(karlin)
#     print(f"{maseq.name} moy = {mean_karlin:.3f}  min = {min_karlin:.3f}  max = {max_karlin:.3f}")
#     entropie = s.get_entropy(500)
#     # Fonction pour calculer la moyenne autour d'une position
#     def average_around(position, window_size=500):
#         start = max(0, position - window_size)
#         end = min(len(karlin), position + window_size)
#         return sum(karlin[start:end]) / (end - start)

#     def entropy_max_around(position,entropy, window_size=1000):
#         start = max(0, position - window_size)
#         end = min(len(entropy), position + window_size)
#         return max(entropy[start:end])
#     # Tracé
#     plt.figure(figsize=(12, 4))
#     plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    
#     y_vals = [0]
#     for base in maseq.sequence:
#         y = y_vals[-1]
#         if base == 'A' or base == 'C':
#             y += 1
#         elif base == 'G' or base == 'T':
#             y -= 1
#         y_vals.append(y)

#     if len(y_vals) > 10000:  # Assure-toi qu'il y a assez d'éléments
#         y_vals = y_vals[5000:-5000]

#     prominence = 500  # Ajuster cette valeur pour contrôler la sensibilité
#     distance = 5000  # Distance minimale entre les pics en bases
    
#     # Trouver les pics positifs
#     peaks_pos, _ = find_peaks(y_vals, prominence=prominence, distance=distance)
    
#     valid_peaks = []

#     # Étape 1 : Collecter les pics valides
#     for i in range(1, len(karlin) - 1, 100):
#         if karlin[i] > 1.75 * mean_karlin:
#             avg_around = average_around(i)
#             if avg_around > mean_karlin and entropy_max_around(i, entropie) > 1.97:
#                 is_near_y_peak = any(abs(i - yp) <= 1500 for yp in peaks_pos)
#                 if is_near_y_peak:
#                     valid_peaks.append(i)

#     #  Vérifier qu'on a des pics valides
#     if valid_peaks:
#         valid_peaks.sort()
#         grouped_peaks = []
#         group = [valid_peaks[0]]

#         for i in range(1, len(valid_peaks)):
#             if valid_peaks[i] - valid_peaks[i-1] <= 1500:
#                 group.append(valid_peaks[i])
#             else:
#                 if len(group) >= 3:
#                     grouped_peaks.append(group)
#                 group = [valid_peaks[i]]

#         if len(group) >= 3:
#             grouped_peaks.append(group)

#         # Trouver le groupe le plus dense
#         if grouped_peaks:
#             most_dense_group = max(grouped_peaks, key=len)

#             # Tracer uniquement ce groupe
#             for peak in most_dense_group:
#                 plt.axvline(x=peak, color='orange', linestyle='--', label='Zone la plus dense')
#     else:
#         print("Aucun pic valide trouvé.")

#     # Centromères
#     if maseq.name in centromeres:
#         cent_start, cent_end = centromeres[maseq.name]
#         cent_start = cent_start-5000
#         cent_end = cent_end -5000
#         plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
#         plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    
#     plt.title(f'Karlin {maseq.name}')
#     plt.xlabel('Position in sequence')
#     plt.ylabel('Karlin')
#     plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(f"../output/dataExploration/search/4/{maseq.name}.png")
#     plt.close()

#   # Tracé Zoom
#     plt.figure(figsize=(12, 4))
#     plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    
#     valid_peaks = []

#     # Étape 1 : Collecter les pics valides
#     for i in range(1, len(karlin) - 1, 100):
#         if karlin[i] > 1.75 * mean_karlin:
#             avg_around = average_around(i)
#             if avg_around > mean_karlin and entropy_max_around(i, entropie) > 1.97:
#                 is_near_y_peak = any(abs(i - yp) <= 1500 for yp in peaks_pos)
#                 if is_near_y_peak:
#                     valid_peaks.append(i)

#     # Vérifier qu'on a des pics valides
#     if valid_peaks:
#         valid_peaks.sort()
#         grouped_peaks = []
#         group = [valid_peaks[0]]

#         for i in range(1, len(valid_peaks)):
#             if valid_peaks[i] - valid_peaks[i-1] <= 1500:
#                 group.append(valid_peaks[i])
#             else:
#                 if len(group) >= 3:
#                     grouped_peaks.append(group)
#                 group = [valid_peaks[i]]

#         if len(group) >= 3:
#             grouped_peaks.append(group)

#         # Trouver le groupe le plus dense
#         if grouped_peaks:
#             most_dense_group = max(grouped_peaks, key=len)

#             # Tracer uniquement ce groupe
#             for peak in most_dense_group:
#                 plt.axvline(x=peak, color='orange', linestyle='--', label='Zone la plus dense')
#     else:
#         print("Aucun pic valide trouvé.")

#     # Centromères
#     if maseq.name in centromeres:
#         cent_start, cent_end = centromeres[maseq.name]
#         cent_start = cent_start-5000
#         cent_end = cent_end -5000
#         plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
#         plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
#                 # Zoom autour du centromère
#         margin = 10000  # marge de 10kb autour du centromère
#         zoom_start = max(0, cent_start - margin)
#         zoom_end = min(len(karlin), cent_end + margin)
#         plt.xlim(zoom_start, zoom_end)
    
#     plt.title(f'Karlin {maseq.name}')
#     plt.xlabel('Position in sequence')
#     plt.ylabel('Karlin')
#     plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(f"../output/dataExploration/search/4/{maseq.name}_.png")
#     plt.close()


    #Etape 3
#     s = DNA(maseq.sequence)
    
#     s.window = 100
#     df_karlin = pd.read_csv(f"../data/karlin/{maseq.name}_karlin_100.csv")
#     karlin = df_karlin['Karlin'].tolist()
#     karlin = karlin[5000:-5000]
#     positions = range(len(karlin))

#     mean_karlin = sum(karlin) / len(karlin)
#     min_karlin = min(karlin)
#     max_karlin = max(karlin)
#     print(f"{maseq.name} moy = {mean_karlin:.3f}  min = {min_karlin:.3f}  max = {max_karlin:.3f}")
#     entropie = s.get_entropy(500)
#     # Fonction pour calculer la moyenne autour d'une position
#     def average_around(position, window_size=500):
#         start = max(0, position - window_size)
#         end = min(len(karlin), position + window_size)
#         return sum(karlin[start:end]) / (end - start)

#     def entropy_max_around(position,entropy, window_size=1000):
#         start = max(0, position - window_size)
#         end = min(len(entropy), position + window_size)
#         return max(entropy[start:end])
#     # Tracé
#     plt.figure(figsize=(12, 4))
#     plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    
#     y_vals = [0]
#     for base in maseq.sequence:
#         y = y_vals[-1]
#         if base == 'A' or base == 'C':
#             y += 1
#         elif base == 'G' or base == 'T':
#             y -= 1
#         y_vals.append(y)

#     if len(y_vals) > 10000:  # Assure-toi qu'il y a assez d'éléments
#         y_vals = y_vals[5000:-5000]

#     prominence = 500  # Ajuster cette valeur pour contrôler la sensibilité
#     distance = 5000  # Distance minimale entre les pics en bases
    
#     # Trouver les pics positifs
#     peaks_pos, _ = find_peaks(y_vals, prominence=prominence, distance=distance)
    
#     # Ajouter des lignes verticales aux pics
#     for i in range(1, len(karlin) - 1,100):  # Ne pas inclure les bords
#         # Vérification si c'est un pic (2 fois supérieur à la moyenne générale)
#         if karlin[i] > 1.75 * mean_karlin:
#             # Vérification de la condition de la moyenne autour de 2000
#             avg_around = average_around(i)
#             if avg_around > mean_karlin and  entropy_max_around(i,entropie) > 1.97:
#                 # Vérifie si ce pic est proche d'un pic Y_skew
#                 is_near_y_peak = False
#                 is_near_y_peak = any(abs(i - yp) <= 1500 for yp in peaks_pos)
#                 if is_near_y_peak:

#                     # Affiche la ligne et imprime l’info
#                     plt.axvline(x=i, color='orange', linestyle='--', label='Pic local')

#     # Centromères
#     if maseq.name in centromeres:
#         cent_start, cent_end = centromeres[maseq.name]
#         cent_start = cent_start-5000
#         cent_end = cent_end -5000
#         plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
#         plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    
#     plt.title(f'Karlin {maseq.name}')
#     plt.xlabel('Position in sequence')
#     plt.ylabel('Karlin')
#     plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(f"../output/dataExploration/search/3/{maseq.name}.png")
#     plt.close()

#   # Tracé Zoom
#     plt.figure(figsize=(12, 4))
#     plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    
#     # Ajouter des lignes verticales aux pics
#     for i in range(1, len(karlin) - 1,100):  # Ne pas inclure les bords
#         # Vérification si c'est un pic (2 fois supérieur à la moyenne générale)
#         if karlin[i] > 1.75 * mean_karlin:
#             # Vérification de la condition de la moyenne autour de 2000
#             avg_around = average_around(i)
#             if avg_around > mean_karlin and  entropy_max_around(i,entropie) > 1.97:
#                 # Vérifie si ce pic est proche d'un pic Y_skew
#                 is_near_y_peak = False
#                 is_near_y_peak = any(abs(i - yp) <= 1500 for yp in peaks_pos)
#                 if is_near_y_peak:

#                     # Affiche la ligne et imprime l’info
#                     plt.axvline(x=i, color='orange', linestyle='--', label='Pic local')

#     # Centromères
#     if maseq.name in centromeres:
#         cent_start, cent_end = centromeres[maseq.name]
#         cent_start = cent_start-5000
#         cent_end = cent_end -5000
#         plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
#         plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
#                 # Zoom autour du centromère
#         margin = 10000  # marge de 10kb autour du centromère
#         zoom_start = max(0, cent_start - margin)
#         zoom_end = min(len(karlin), cent_end + margin)
#         plt.xlim(zoom_start, zoom_end)
    
#     plt.title(f'Karlin {maseq.name}')
#     plt.xlabel('Position in sequence')
#     plt.ylabel('Karlin')
#     plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(f"../output/dataExploration/search/3/{maseq.name}_.png")
#     plt.close()

   


#     #Etape 2
#     s = DNA(maseq.sequence)
    
#     s.window = 100
#     df_karlin = pd.read_csv(f"../data/karlin/{maseq.name}_karlin_100.csv")
#     karlin = df_karlin['Karlin'].tolist()
#     karlin = karlin[5000:-5000]
#     positions = range(len(karlin))

#     mean_karlin = sum(karlin) / len(karlin)
#     min_karlin = min(karlin)
#     max_karlin = max(karlin)
#     print(f"{maseq.name} moy = {mean_karlin:.3f}  min = {min_karlin:.3f}  max = {max_karlin:.3f}")
#     entropie = s.get_entropy(500)
#     # Fonction pour calculer la moyenne autour d'une position
#     def average_around(position, window_size=500):
#         start = max(0, position - window_size)
#         end = min(len(karlin), position + window_size)
#         return sum(karlin[start:end]) / (end - start)

#     def entropy_max_around(position,entropy, window_size=1000):
#         start = max(0, position - window_size)
#         end = min(len(entropy), position + window_size)
#         return max(entropy[start:end])
#     # Tracé
#     plt.figure(figsize=(12, 4))
#     plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    
#     # Ajouter des lignes verticales aux pics
#     for i in range(1, len(karlin) - 1,100):  # Ne pas inclure les bords
#         # Vérification si c'est un pic (2 fois supérieur à la moyenne générale)
#         if karlin[i] > 1.75 * mean_karlin:
#             # Vérification de la condition de la moyenne autour de 2000
#             avg_around = average_around(i)
#             if avg_around > mean_karlin and  entropy_max_around(i,entropie) > 1.97:
#                 plt.axvline(x=i, color='orange', linestyle='--', label='Pic local')

#     # Centromères
#     if maseq.name in centromeres:
#         cent_start, cent_end = centromeres[maseq.name]
#         cent_start = cent_start-5000
#         cent_end = cent_end -5000
#         plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
#         plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    
#     plt.title(f'Karlin {maseq.name}')
#     plt.xlabel('Position in sequence')
#     plt.ylabel('Karlin')
#     plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(f"../output/dataExploration/search/2/{maseq.name}.png")
#     plt.close()

#   # Tracé Zoom
#     plt.figure(figsize=(12, 4))
#     plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    
#     # Ajouter des lignes verticales aux pics
#     for i in range(1, len(karlin) - 1,100):  # Ne pas inclure les bords
#         # Vérification si c'est un pic (2 fois supérieur à la moyenne générale)
#         if karlin[i] > 1.75 * mean_karlin:
#             # Vérification de la condition de la moyenne autour de 2000
#             avg_around = average_around(i)
#             if avg_around > mean_karlin and entropy_max_around(i,entropie) > 1.97:
#                 plt.axvline(x=i, color='orange', linestyle='--', label='Pic local')

#     # Centromères
#     if maseq.name in centromeres:
#         cent_start, cent_end = centromeres[maseq.name]
#         cent_start = cent_start-5000
#         cent_end = cent_end -5000
#         plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
#         plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
#                 # Zoom autour du centromère
#         margin = 10000  # marge de 10kb autour du centromère
#         zoom_start = max(0, cent_start - margin)
#         zoom_end = min(len(karlin), cent_end + margin)
#         plt.xlim(zoom_start, zoom_end)
    
#     plt.title(f'Karlin {maseq.name}')
#     plt.xlabel('Position in sequence')
#     plt.ylabel('Karlin')
#     plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(f"../output/dataExploration/search/2/{maseq.name}_.png")
#     plt.close()

   
    # dnaFlex = s.get_dna_flexibility()
    # # Tracé
    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, dnaFlex, label=f'DNA Flex - {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'DNA Flex {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('DNA Flex')
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/dna_Flex/{maseq.name}.png")
    # plt.close()

    # information_entropy = s.get_informational_entropy()
    # # Tracé
    # plt.figure(figsize=(12, 4))
    # plt.plot(positions, information_entropy, label=f'informational entropy- {maseq.name}', color='teal')
    # if maseq.name in centromeres:
    #     cent_start, cent_end = centromeres[maseq.name]
    #     plt.axvline(x=cent_start, color='red', linestyle='--', label='Centromere Start')
    #     plt.axvline(x=cent_end, color='purple', linestyle='--', label='Centromere End')
    # plt.title(f'informational entropy {maseq.name}')
    # plt.xlabel('Position in sequence')
    # plt.ylabel('informational entropy')
    # plt.axhline(0, color='gray', linestyle='--')  # ligne à y=0 pour référence
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f"../output/dataExploration/information_entropy/{maseq.name}.png")
    # plt.close()
