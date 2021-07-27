import numpy as np
from sklearn.manifold import TSNE

fname = 'C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor Arbeit\\Data\\precomputed\\10s_data.npy'
x_raw_mat = np.load(fname)

contig_lengths = np.sum(x_raw_mat, 1)
x_mat = (x_raw_mat.T / contig_lengths).T

data = TSNE(n_components=2).fit_transform(x_mat)

print(x_raw_mat, "\n", x_mat, "\n", data)
