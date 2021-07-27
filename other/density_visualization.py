import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from KDEpy import FFTKDE  # Density ?
from skimage.feature import peak_local_max
from sklearn.manifold import TSNE

################################
#
# Features:
#
# - Analyses and visualizes 2D point data based on their density using FFTKDE.
# - Shows interactive live changes of contour numbers and bandwidth.
# - Allows export of data points in a chosen cluster.
#
# Usage notes:
#
# It is possible to adjust the amount of visualized contours and the bandwidth of the FFTKDE.
# Afterwards you may process the points and then select contour regions to be exported to a file.
# Note that while multiple regions might be colored, only the connected region that is selected will be exported.
#
################################

################################
# Variables
################################

# 2D data points in variable "data"
# in shape (n,2) where n is the number of datapoints

# prepare data - example usage    
fname = 'C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor Arbeit\\Data\\precomputed\\10s_data.npy'
x_raw_mat = np.load(fname)
(n, d) = x_raw_mat.shape
contig_lengths = np.sum(x_raw_mat, 1)
x_mat = (x_raw_mat.T / contig_lengths).T  # relative frequencies (alternative normalization possible, e.g. centered log-ratio)
# labels, only use for control/benchmarking!
labs_ref = np.load('C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor Arbeit\\Data\\precomputed\\10s_labels.npy')
z_mat = TSNE(n_components=2).fit_transform(x_mat)  # t-SNE

# final feature vectors 
data = z_mat
# data labels in variable "data_labels" in shape (n,), can be anything informative 
data_labels = labs_ref  # contig_lengths # labs_ref # np.arange(len(data)) # simple indices

# resolution in one dimension
grid_points = 100

# filename to be created or overwritten
export_filename = "data_output_indices.npy"

show_scatterplot = True
show_labels = True
show_local_maxima = True

################################

min_contour = 3
max_contour = 50
min_bandwidth = 0.01
max_bandwidth = 100

################################
# End of variables
################################

contours = (min_contour + max_contour) // 2
data_contour = np.zeros(len(data))
data_contoured = False

point_col = np.zeros((data.shape[0], 3))

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.20)

prev_mouse_x = -1
prev_mouse_y = -1


def get_selected_path(event):
    if not data_contoured:
        return
    if not event.xdata or not event.ydata:
        return
    i = 0
    for paths in cont_data.collections:
        if not paths.contains(event)[0] or i == len(cont_data.collections) - 1:
            if not paths.contains(event)[0]:
                cur = cont_data.collections[i - 1]
            else:
                cur = cont_data.collections[i]
            if i > 1:
                return cur
            else:
                return None
        i += 1


def cursor_press(event):
    global prev_mouse_x
    global prev_mouse_y
    prev_mouse_x = event.xdata
    prev_mouse_y = event.ydata


def cursor_release(event):
    if prev_mouse_x != event.xdata or prev_mouse_y != event.ydata:
        return
    selected_path = get_selected_path(event)
    if not selected_path:
        return
    for path in selected_path.get_paths():
        if path.contains_point((event.xdata, event.ydata), radius=0.5):
            global data_output_indices
            data_output_indices = np.nonzero(path.contains_points(data))[0]
            np.save(export_filename, data_output_indices)
            print("saved " + export_filename)


def cursor_hover(event):
    if show_scatterplot and show_labels:
        contains, index = scatter_data.contains(event)
        if contains:
            current_index = index["ind"][0]
            label = data_labels[current_index]
            annot.set_text(label)
            annot.xy = scatter_data.get_offsets()[current_index]
            annot.set_visible(True)
        else:
            annot.set_visible(False)

    for i in range(len(cont_data.collections)):
        cont_data.collections[i].set_color(contour_colors[i])
    selected_path = get_selected_path(event)
    if selected_path:
        for path in selected_path.get_paths():
            if path.contains_point((event.xdata, event.ydata), radius=0.5):
                selected_path.set_edgecolor("green")
    fig.canvas.draw_idle()


fig.canvas.mpl_connect('button_press_event', cursor_press)
fig.canvas.mpl_connect('button_release_event', cursor_release)
fig.canvas.mpl_connect('motion_notify_event', cursor_hover)


def replot():
    global cont_data
    cont_data = ax.contourf(x, y, z, contours, cmap='RdBu_r')
    global contour_colors
    contour_colors = []
    for paths in cont_data.collections:
        contour_colors.append(paths.get_facecolor()[0])
    if show_scatterplot:
        global scatter_data
        scatter_data = ax.scatter(data[:, 0], data[:, 1], marker=".", s=2, c=data_contour)
    if show_local_maxima:
        peakdata = peak_local_max(z, threshold_rel=0.02)
        ax.scatter(x[peakdata[:, 1]], y[peakdata[:, 0]], marker=10, c="orange")
    global annot
    annot = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
                        bbox=dict(fc="w"))
    annot.set_visible(False)


def drawkde(bw):
    kde = FFTKDE(norm=1, bw=bw)
    grid, points = kde.fit(data).evaluate(grid_points)
    global x, y, z
    x, y = np.unique(grid[:, 0]), np.unique(grid[:, 1])
    z = points.reshape(grid_points, grid_points).T
    replot()


drawkde(1)
min_b = np.log10(min_bandwidth)
max_b = np.log10(max_bandwidth)
s_bandwidth = Slider(plt.axes([0.25, 0.02, 0.35, 0.03]), "Bandwidth", min_b, max_b, (min_b + max_b) / 2)
v_bandwidth = 10 ** (s_bandwidth.val)
s_bandwidth.valtext.set_text(f"{v_bandwidth:.3f}")

s_contours = Slider(plt.axes([0.25, 0.08, 0.35, 0.03]), "Countour Lines", min_contour, max_contour, contours, valfmt="%i")

b_process = Button(plt.axes([0.70, 0.04, 0.18, 0.075]), "Process Points")


def update_band(val):
    v_bandwidth = 10 ** val
    s_bandwidth.valtext.set_text(f"{v_bandwidth:.3f}")
    global data_contour
    global data_contoured
    data_contour = np.zeros(len(data))
    data_contoured = False
    ax.cla()
    drawkde(v_bandwidth)


s_bandwidth.on_changed(update_band)


def update_cont(val):
    global contours
    contours = int(val)
    global data_contour
    global data_contoured
    data_contour = np.zeros(len(data))
    data_contoured = False
    ax.cla()
    replot()


s_contours.on_changed(update_cont)


def button_process(event):
    global data_contour
    global data_contoured
    data_contour = np.zeros(len(data))
    for paths in cont_data.collections:
        for path in paths.get_paths():
            newadd = path.contains_points(data)
            data_contour += newadd
    ax.cla()
    replot()
    data_contoured = True


b_process.on_clicked(button_process)

plt.show()
