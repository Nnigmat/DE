import matplotlib.pyplot as plt
from math import exp
from numpy import linspace

err = []
x =  linspace(-5, 0, 1000)
for i in x:
    err.append(1/(-4.5 - i) + exp(i))

plt.plot(x, err)
plt.show()
