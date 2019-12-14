import numpy as np
import scipy
import random
import matplotlib.pyplot as plt

def rando():
    return random.random()

X = np.zeros(100, dtype=int)
T = np.zeros(100, dtype=int)

for i in range(10000):
    r = int(rando()*100)
    X[r] += 1

plt.xlabel('N')
plt.ylabel('Height')
plt.plot(X, 'r-')
plt.savefig('plot.png')