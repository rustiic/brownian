import numpy as np
import scipy
import random
import matplotlib.pyplot as plt



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
    
    # A remnant from the first version where I'd used np arrays
    #X = np.zeros((N,2), float)  

    

    theta1 = 0         # initial conditions (angle)
    theta2 = 0
    
    evo11   = np.sqrt(2*DT*dt)
    evo12   = np.sqrt(2*DR*dt)
    torque1 = dt*W
    evoc1   = dt*V

    evo21    = np.sqrt(2*DT*dt)
    evo22    = np.sqrt(2*DR*dt)
    torque2  = dt*W
    evoc2    = dt*V

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

    # So I use random numbers a lot:
    #   The first returns either a 1 or -1
    #   The second returns a random no multiplied to -1 or 1
    def rando_sign():
        return ((-1)**(int((random.random()*10)%2)))

    def rando():
        return random.random()*rando_sign()
    
    #plt.ion()


    #Obstacles, here tiny cones which will trap the swimmers


    def dist(a, b):
        return np.sqrt(a**2 + b**2)


    def solver(x1, y1, x2, y2):
        a1 = x2 - x1
        b1 = y2 - y1
        c1 = a1 / b1
        a2 = c1**2 + 1
        b2 = 2*c1*(x1 - c1*y1)
        c2 = x1**2 + (c1*y1)**2 - 2*c1*y1*x1 - R**2

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


    for i in range(N):

        #'''
        #Particle 1
        delta_x1 += evo11*rando()
        delta_y1 += evo11*rando()
        theta1 += evo12*random.random() #*rando_sign()
        theta1 += torque1
        delta_x1 += evoc1*np.cos(theta1)*rando_sign()
        delta_y1 += evoc1*np.sin(theta1)*rando_sign()

        # Now for reflection from the walls
        if np.sqrt(delta_x1**2 + delta_y1**2) > L:
            (delta_x1, delta_y1) = solver(xi, yi, delta_x1, delta_y1)
        
        # The following code block is for another particle moving in the medium
        #   simultaneously, complete with collisions. 
        #   Since I later decided to make EVERY MODIFICATION TO THE SAME GODDAMN
        #   FILE, I had to comment whole blocks out.
        
        #This block is for a square box
        '''
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

    #this block is for plotting the graph at the end of runtime
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.plot(X1, Y1, 'r-')
    plt.axis(( -scale, scale, -scale, scale))
    plt.savefig('active_swimmer_plot.png')

    return


homogenous(N=10000, P=200 , dt=0.001, R=0.000001, T=300, Bt=400, eta=0.001, V=0.00001, W=0)
