import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.axes import Axes
from scipy.spatial import ConvexHull

fig, (ax1, ax2) = plt.subplots(2)

ax1: Axes = ax1
ax2: Axes = ax2

# x_vec = [10, 200, 250, 800]
# y_vec = [1, 100, 200, 600]
#
# # ax.plot(x_vec, y_vec)
#
# ellipse = Ellipse((300, 600), width=400, height=100, edgecolor='black', facecolor='green', linewidth=2)
# # ax.add_patch(ellipse)
#
# # Get the path
# path = ellipse.get_path()
# # Get the list of path codes
# codes = path.codes
# # Get the list of path vertices
# vertices = path.vertices.copy()
# # Transform the vertices so that they have the correct coordinates
# vertices = ellipse.get_patch_transform().transform(vertices)
#
# # Do your transforms here
# vertices[0] = vertices[0] - 100
#
# # Create a new Path
# modified_path = Path(vertices, codes)
# # Create a new PathPatch
# modified_ellipse = PathPatch(modified_path, edgecolor='black', facecolor='black', linewidth=2)
#
# # Add the modified ellipse
# # ax.add_patch(modified_ellipse)
#
# elli = Ellipse((300, 300), width=400, height=100, edgecolor='black', facecolor='green', linewidth=2)
#
# path = elli.get_path()
# codes = path.codes
# print(codes, codes.shape)
# vert = path.vertices.copy()
# # vert = elli.get_patch_transform().transform(vert)
#
# # print(vert, vert.shape)
# # print("x:", vert[:, 0])
# # x = vert[:, 0]
# # y = vert[:, 1]
# # print("y:", y)
# # z = np.meshgrid(x, y)
# # print("Z:", z)
# #
# # ax.plot(vert)
#
#
#
# cont = [[1, 2, 3, 4, 5, 4, 3, 2], [5, 7, 9, 7, 5, 3, 1, 3]]
#
# delta = 0.5
# y = np.arange(-3.0, 3.0, delta)
# print(y, "Dim:", y.shape)
# x = np.arange(0, 240)
#
# # ax.plot(x, y)
#
# x = np.arange(-3.0, 3.0, delta)
# y = np.arange(-2.0, 2.0, delta)
# X, Y = np.meshgrid(x, y)
#
# # ax.plot(X, Y)
# print(X.shape, Y.shape, vert.shape)
#
# ax.contour(cont, 30)
#
# Z1 = np.exp(-X ** 2 - Y ** 2)
# Z2 = np.exp(-(X - 1) ** 2 - (Y - 1) ** 2)
# Z = (Z1 - Z2) * 2
# print(Z1, Z1.shape)
# print(Z.shape)
#
# CS = ax.contour(X, Y, Z)


points = np.random.rand(30, 2)  # 30 random points in 2-D
hull = ConvexHull(points)  # Indices of hull points
print(hull.vertices)

ax1.plot(points[:, 0], points[:, 1], 'o')
ax2.plot(points[:, 0], points[:, 1], 'o')
# plot convex hull polygon
ax1.plot(points[hull.vertices, 0], points[hull.vertices, 1], 'r--', lw=2)
# plot convex full vertices
ax1.plot(points[hull.vertices[0], 0], points[hull.vertices[0], 1], 'ro')

hull_vertices = points[hull.vertices]
n = len(hull_vertices)
codes = [1]
i = 1
while i < n - 1:
    codes.append(4)
    i += 1
codes.append(79)
path = Path(hull_vertices, closed=True)
patch = PathPatch(path, facecolor='red', lw=1, edgecolor='black', fill=False)

ax2.add_patch(patch)

plt.show()
