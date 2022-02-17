import colorsys
import matplotlib.pyplot as plt
import matplotlib.colors as colors

hue = [0.78, 0.4, 0.21, 0.98]

rgb_col = []

i = 0
while i < len(hue):
    rgb_col.append(colors.hsv_to_rgb([hue[i], 1.0, 1.0]))
    i += 1

x_data = [1, 1, 2, 2]
y_data = [1, 2, 2, 1]
s = 10
size = [s * 4 ** 1, s * 4 ** 2, s * 4 ** 3, s * 4 ** 4]
plt.scatter(x_data, y_data, s=size, c=rgb_col)
plt.show()
