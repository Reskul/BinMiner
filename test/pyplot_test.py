#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math

pi = math.pi

x_list = np.arange(-2*pi,2*pi,0.1)
y_list = [math.cos(x) for x in x_list]

fig = plt.figure(1)

plot = fig.add_subplot(111)

plt.plot(x_list,y_list)

plot.tick_params(axis='x', labelsize=30)

plt.grid()
plt.title('Change label axis font size in matplotlib')

plt.show()
