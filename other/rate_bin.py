from lib import FastaReader
import json
import numpy as np
from matplotlib import pyplot as plt

if __name__ == '__main__':
    bin_filepath = "C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\Ergebnisse\\030622_GUI+Sequenzen\\Bin01_n176.fasta"
    labels_filepath = "C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\Datensatz_25s\\25s_mycc\\member.txt"

    bin_file = open(bin_filepath, 'r')
    labels_file = open(labels_filepath, 'r')

    reader = FastaReader(FastaReader.BINMINER)
    header_vec = reader.read_header_only(bin_file)

    labels_vec = labels_file.read().split('\n')
    print(labels_vec[0].split('\t'))
    organisms_dict = {}
    for l in labels_vec:
        info = l.split('\t')
        contig = info[0].strip()
        org = info[1].strip()

        if org not in organisms_dict:
            organisms_dict[org] = [contig]
        else:
            organisms_dict[org].append(contig)

    # json_file = open("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\Ergebnisse\\030622_GUI+Sequenzen\\Bin01_n176_AUSWERUNG_DICT.json", 'w')
    # json.dump(organisms_dict, json_file)

    organisms_vec = list(organisms_dict.keys())
    count_vec_abs = np.zeros(len(organisms_vec))
    count_vec_rel = np.zeros(len(organisms_vec))

    # check header of given bin
    for h in header_vec:
        o = h.organism
        o_idx = organisms_vec.index(o)
        count_vec_abs[o_idx] += 1

    sum_contigs = sum(count_vec_abs)
    i = 0
    while i < len(count_vec_abs):
        if count_vec_abs[i] != 0:
            count_vec_rel[i] = count_vec_abs[i] / len(organisms_dict[organisms_vec[i]])
            print(count_vec_rel[i])
        i += 1

    result = [None] * (len(count_vec_rel) * 2)
    result[::2] = organisms_vec
    result[1::2] = count_vec_rel
    print(result)

    fig, ax = plt.subplots()
    ax.bar(organisms_vec, count_vec_rel)
    plt.yticks(ticks=[])
    plt.xticks(rotation=90, ha='center')
    plt.tight_layout()
    plt.savefig("C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\Ergebnisse\\030622_GUI+Sequenzen\\Bin01_n176_plot.png")
    plt.show()
