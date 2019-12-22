import numpy as np
#import random
from random import SystemRandom, random, randint
import matplotlib.pyplot as plt

def Sample(size, moment):
    # size specifies the size of the arrays
    # moment specifies the # of moments to be calculated

    cryptogen = SystemRandom()
    randon, crypton = np.zeros((moment), dtype=float), np.zeros((moment), dtype=float)
    #silly = np.zeros(moment, dtype=float)

    def rando():
        return randint(0, 1000)

    def crypto():
        return cryptogen.randint(0, 1000)


    X = np.zeros(size, dtype=float)
    Y = np.zeros(size, dtype=float)
    # rand_moments = np.zeros(size, dtype=float)
    # crypt_moments = np.zeros(size, dtype=float)
    #Z = np.zeros(size, dtype=float)

    for i in range(size):
        x = rando()
        y = crypto()
        #z = randint(500, 999)
        X[x-1] += 1/size
        Y[y-1] += 1/size
        #Z[z-1] += 1/size


    for i in range(moment):
        for j in range(size):
            randon[i] += j**(i+1)*X[j]
            crypton[i] += j**(i+1)*Y[j]

    
    # f, axarr = plt.subplots(2, sharex=True)
    # axarr[0].plot(X, 'b-')
    # axarr[0].set_title('random.random()')
    # axarr[1].plot(Y, 'r-')
    # axarr[1].set_title('SystemRandom')
    # plt.savefig('plot.png')
    
    file = open('out', 'w')

    for i in range(moment):
        file.write("\n {}th moment of \n Rando : {} \n Crypto: {}".format( i+1, randon[i], crypton[i]))
    file.close()

Sample(1000, 20)
