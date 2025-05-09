from sequana import DNA
from sequana import FastA
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pywt
from scipy.signal import find_peaks

f = FastA("../data/Fasta/TriTrypDB-68_LmajorFriedlin_Genome.fasta")
df = pd.read_csv("../data/Centromere_Positions/centromeres.csv")
centromeres = {
    row['Chromosome']: (row['Centromere_Start'], row['Centromere_End'])
    for _, row in df.iterrows()
}

for maseq in f:
    #Etape 3
    s = DNA(maseq.sequence)
    
    s.window = 100
    df_karlin = pd.read_csv(f"../data/karlin/{maseq.name}_karlin_100.csv")
    karlin = df_karlin['Karlin'].tolist()
    karlin = karlin[5000:-5000]
    positions = range(len(karlin))

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

    # Supprimer les 5000 premières et dernières valeurs
    if len(y_vals) > 10000:  # Assure-toi qu'il y a assez d'éléments
        y_vals = y_vals[5000:-5000]
    peaks, _ = find_peaks(y_vals, distance=500, prominence=10)

    # Ajouter des lignes verticales aux pics
    for i in range(1, len(karlin) - 1,100):  # Ne pas inclure les bords
        # Vérification si c'est un pic (2 fois supérieur à la moyenne générale)
        if karlin[i] > 1.75 * mean_karlin:
            # Vérification de la condition de la moyenne autour de 2000
            avg_around = average_around(i)
            if avg_around > mean_karlin and  entropy_max_around(i,entropie) > 1.97:
                # Vérifie si ce pic est proche d'un pic Y_skew
                is_near_y_peak = False
                is_near_y_peak = any(abs(i - yp) <= 100 for yp in peaks)
                if is_near_y_peak:

                    # Affiche la ligne et imprime l’info
                    plt.axvline(x=i, color='orange', linestyle='--', label='Pic local')

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
    plt.savefig(f"../output/dataExploration/search/3/{maseq.name}.png")
    plt.close()

  # Tracé Zoom
    plt.figure(figsize=(12, 4))
    plt.plot(positions, karlin, label=f'Karlin - {maseq.name}', color='teal')
    
    # Ajouter des lignes verticales aux pics
    for i in range(1, len(karlin) - 1,100):  # Ne pas inclure les bords
        # Vérification si c'est un pic (2 fois supérieur à la moyenne générale)
        if karlin[i] > 1.75 * mean_karlin:
            # Vérification de la condition de la moyenne autour de 2000
            avg_around = average_around(i)
            if avg_around > mean_karlin and  entropy_max_around(i,entropie) > 1.97:
                # Vérifie si ce pic est proche d'un pic Y_skew
                is_near_y_peak = False
                is_near_y_peak = any(abs(i - yp) <= 100 for yp in peaks)
                if is_near_y_peak:

                    # Affiche la ligne et imprime l’info
                    plt.axvline(x=i, color='orange', linestyle='--', label='Pic local')

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
    plt.savefig(f"../output/dataExploration/search/3/{maseq.name}_.png")
    plt.close()

   


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
