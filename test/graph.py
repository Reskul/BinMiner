import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    print(np.log2(40), np.log2(40) ** 2)

    n = 51
    x = [i for i in range(n)[1:]]
    x_mean = sum(x) / len(x)
    print(x_mean)
    print(x)
    y = [(i ** 2) / x_mean for i in x]
    y2 = [i ** 2 for i in x]

    plt.plot(x, y)
    plt.plot(x, y2)
    plt.show()
