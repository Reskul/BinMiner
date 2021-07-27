import numpy as np
from KDEpy import FFTKDE
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from skimage.feature import peak_local_max

fname = 'C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor Arbeit\\Data\\precomputed\\10s_data.npy'

# TEST DATA LOAD IN AND TSNE
x_raw_mat = np.load(fname)

contig_lengths = np.sum(x_raw_mat, 1)
x_mat = (x_raw_mat.T / contig_lengths).T

data = TSNE(n_components=2).fit_transform(x_mat)

print(x_raw_mat, "\n", x_mat, "\n", data)


# TEST KDE
bw = 1

data_contour = np.zeros(len(data))
#data = np.random.randn(2 ** 4)
# (1) Automatic bw selection using Improved Sheather Jones (ISJ)
#x, y = FFTKDE(bw='ISJ').fit(data).evaluate(grid_points)
# (2) Explicit choice of kernel and bw (standard deviation of kernel)
#x, y = FFTKDE(kernel='triweight', bw=0.5).fit(data).evaluate()
weights = data + 10
# (3) Using a grid and weights for the data
#y = FFTKDE(kernel='epa', bw=0.5).fit(data, weights).evaluate(x)
# (4) If you supply your own grid, it must be equidistant
#y = FFTKDE().fit(data)(np.linspace(-10, 10, num=2 ** 12))

#print("After transformation")
#print("X:\n", x)
#print("Y:\n", y)
#print("Z:\n", z)

  # std value? i guess
fig, ax = plt.subplots()
# Contours
grid_points = 100
grid, points = FFTKDE(norm=1, bw=bw).fit(data).evaluate(grid_points)
x, y = np.unique(grid[:, 0]), np.unique(grid[:, 1])
z = points.reshape(grid_points, grid_points).T
contours = 26
cont_data = ax.contourf(x, y, z, contours, cmap='RdBu_r')

print("X:\n", grid)
print("Y:\n", points)
# Datapoints
scatter_data = ax.scatter(data[:, 0], data[:, 1], marker=".", s=2, c=data_contour)
# Peaks
peakdata = peak_local_max(z, threshold_rel=0.02)
ax.scatter(x[peakdata[:, 1]], y[peakdata[:, 0]], marker=10, c="orange")

ax.set_title("Data")
plt.show()
