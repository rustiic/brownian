from numba import jit
import numpy as np
#import random
from random import SystemRandom, random, randrange
import matplotlib.pyplot as plt

#The two fucntions get rand# in two different ways
def random_num(a, b):
    return randrange(a, b)

sysrand = SystemRandom()        #old name cryptogen, needed for crypto() to work
def sys_rand_num(a, b):
    return sysrand.randrange(a, b)


#Sampling program for rngs random and SystemRandom
def Sample(size, moment):
    # size specifies the size of the arrays
    # moment specifies the # of moments to be calculated

    #The following line makes arrays for moment calculation
    # randon, crypton = np.zeros((moment), dtype=float), np.zeros((moment), dtype=float)



    X = np.zeros(size, dtype=float)
    Y = np.zeros(size, dtype=float)
    dist = np.zeros(size, dtype=float)
    # rand_moments = np.zeros(size, dtype=float)
    # crypt_moments = np.zeros(size, dtype=float)
    #Z = np.zeros(size, dtype=float)

    for i in range(size*moment):
        X[random_num(0, size)] += 1/(size*moment)
        Y[sys_rand_num(0, size)] += 1/(size*moment)

    '''
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
    '''
    for i in range(size):
        dist[i] = np.exp(Y[i])

    return dist





#moved to here from brownian script
def rando_sign():
    return ((-1)**(int((random()*10)%2)))

def rando():
    return random()*rando_sign()

def dist(a, b):
        return np.sqrt(a**2 + b**2)

#The following function was written for reflection at boundaries
def solver(x1, y1, x2, y2):
    a1 = x2 - x1
    b1 = y2 - y1
    c1 = a1 / b1
    a2 = c1**2 + 1
    b2 = 2*c1*(x1 - c1*y1)
    c2 = x1**2 + (c1*y1)**2 - 2*c1*y1*x1 - L**2

    len_path = dist(x2 - x1, y2 - y1)
    
    for i in range(2):
        for j in range(2):
            yc = (-b2 + (-1**(2-i))*np.sqrt(b2**2 - 4*a2*c2)) / (2 * a2)
            xc = (-1**(2-j))*np.sqrt(L**2 - yc**2)
            if dist(xc, yc) == L:
                if len_path == dist(x1 - xc, y1 - yc) + dist(x2 - xc, y2 - yc):
                    break

    v_dot_n = (xc*(xc - x1) + yc*(yc-y1))/(L * dist(xc - x1, yc - y1))

    rx = x1 - 2*v_dot_n*xc
    ry = y1 - 2*v_dot_n*yc

    return (rx, ry)
