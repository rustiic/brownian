import numpy as np
#import random
from random import SystemRandom, random, randrange
import matplotlib.pyplot as plt


size = 1000   #Specifies the size of the arrays
passes  = 5

cryptogen = SystemRandom()
randon, crypton = np.zeros(passes - 1, dtype=int), np.zeros(passes - 1, dtype=int)

def rando():
    return randrange(size)

def crypto():
    return cryptogen.randrange(size)


X = np.zeros((passes, size), dtype=int)
Y = np.zeros((passes, size), dtype=int)
t = np.zeros(size, dtype=int)


for i in range(size*1000):
    for j in range(passes):
        x = rando()
        y = crypto()
        X[j][x] += 1
        Y[j][y] += 1

#'''
for i in range(passes - 1):
    for j in range(size):
        X[i + 1][j] -= X[0][i]
        X[i + 1][j] = abs(X[i + 1][j])
        Y[i + 1][j] -= Y[0][i]
        Y[i + 1][j] = abs(Y[i + 1][j])
#'''

for i in range(size):
    X[0][i] = 0
    Y[0][i] = 0
    t[i] = 1
    for j in range(passes - 1):
        randon[j]  += X[j+1][i]
        crypton[j] += Y[j+1][i]

X = np.transpose(X)
Y = np.transpose(Y)

f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(X, 'b-')
axarr[0].set_title('random.random()')
axarr[1].plot(Y, 'r-')
axarr[1].set_title('SystemRandom')
plt.savefig('plot.png')

for i in range(passes - 1):
    print("\n {}) Rando {} \t Crypto: {}".format( i+1, randon[i], crypton[i]))