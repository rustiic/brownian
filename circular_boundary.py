from numba import jit
import numpy as np
import scipy
from random import SystemRandom, random, randint
import matplotlib.pyplot as plt
from noise import *



#N  is number of Steps
#P  is the number of particles
#dt is time step
#R  is the radius of the particle
#T  is the temp of the bath
#Bt is the temp of the box, which I never got around to using.
#eta is the viscosity
#V  is the particle speed

#W  is the angular velocity
#W = float(input("Enter anglular Velocity : "))
def homogenous(N, P, dt, R, T, Bt, eta, V, W):
    

    kB = 1.38e-23		    # Boltzmann constant [J/K]
    gamma = 6*np.pi*R*eta	# friction coefficient [Ns/m]
    DT = kB*T/gamma		    # translational diffusion coefficient [m^2/s]
    DR = 6*DT/(8*R**2)      # rotational diffusion coefficient [rad^2/s]

    theta1 = 0         # initial conditions (angle)
    
    evo11   = np.sqrt(2*DT*dt)
    evo12   = np.sqrt(2*DR*dt)
    torque1 = dt*W
    evoc1   = dt*V

    scale = evo11*50              #Scale for the axes and plot
    L = scale                     #this is the dimension of the box

    track = 100                  #Logest lengths of the tracks of the particle


    delta_x1, delta_y1 = 0.99*scale, -0.25*scale
    delta_x1, delta_y1 = 0, -0.99*scale
    
    X1 = []
    Y1 = []

    x1, y1, xi, yi = 0, 0, 0, 0

    # Output to file, which I was initially using to plot with GNUPlot
    #f = open('data.out', 'w')
    #print("X position \t Y position \n", file = f)
    
    # Interactive on, to make plot interative
    #plt.ion()


    #Obstacles, here tiny cones which will trap the swimmers

    dist = Sample(N, 100)



    for i in range(N):
        
        # delta_x1 += evo11*rando()
        # delta_y1 += evo11*rando()
        delta_x1 += dist[sys_rand_num(0, N)]
        delta_y1 += dist[sys_rand_num(0, N)]
        theta1 += evo12*random() #*rando_sign()
        theta1 += torque1
        delta_x1 += evoc1*np.cos(theta1)*rando_sign()
        delta_y1 += evoc1*np.sin(theta1)*rando_sign()

        #This block is for reflection from circular walls
        #DOES NOT WORK YET
        '''
        if dist(delta_x1, delta_y1) > L:
            (delta_x1, delta_y1) = solver(xi, yi, delta_x1, delta_y1)
        '''

        #This block is for a square box
        #'''
        if abs(delta_x1) > L:
            delta_x1 -= 2*np.sign(delta_x1)*(abs(delta_x1) - L)
            #T = Bt + (T - Bt)*np. 
            #evo21 = evo21*np.sqrt(Bt/T)
        
        if abs(delta_y1) > L:
            delta_y1 -= 2*np.sign(delta_y1)*(abs(delta_y1) - L)
            #evo21 = evo21*np.sqrt(Bt/T)
        #'''

        X1 =  np.append(X1, delta_x1)
        Y1 =  np.append(Y1, delta_y1)

        xi = delta_x1
        yi = delta_y1


        #print("%s \t %s \n" %(delta_x1, delta_y1), file=f)

        # This following block of code plots and erases tracks interatively
        #'''
        if i > (track):
            plt.clf()
            plt.plot(X1[i-track:i], Y1[i-track:i], 'b-')
            plt.axis(( -scale, scale, -scale, scale))
            plt.pause(0.0001)
        else:
            plt.plot(X1, Y1, 'b-')
            plt.axis(( -scale, scale, -scale, scale))
            plt.pause(0.00001)
        
        #'''
    
    file = open('out.txt', 'w')
    # file.write("X \t \t Y \n")
    for i in range(N):
        file.write("{:.4E} \t {:.4E} \n".format(X1[i], Y1[i]))
    file.close()

    #this block is for plotting the graph at the end of runtime
    # plt.xlabel('X')
    # plt.ylabel('Y')
    plt.plot(X1, 'r-')
    # plt.plot(Y1, 'b-')
    # plt.axis(( -scale, scale, -scale, scale))
    plt.savefig('active_swimmer_plot.png')

    return


homogenous(N=10000, P=200 , dt=0.001, R=0.000001, T=300, Bt=400, eta=0.001, V=0.00001, W=0)
