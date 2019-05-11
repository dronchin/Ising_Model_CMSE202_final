from math import exp
from random import randrange,choice,random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from copy import deepcopy
import imageio
import os



#n is length of x-axis, T_max is highest temp and T_min is lowest temp
#Chunks is the number of pieces you want to slice array into. 
#Temperature will hold at T_min for first chunk, rise, hold at T_max for 1 chunk, fall, then hold at T_min again for last chunk.
#Thus chunks must be greater than 3, and if you want to avoid 1 or 2 rounding errors make n evenly divisible by chunks.
#Returns all of this as an array of length n


def make_gradient(n,T_min,T_max,chunks):
    a = T_min
    Temp_grad = np.zeros(n)
    divider = int(n/chunks)
    move = int(.5*(n-3*divider))


    Temp_grad[0:divider] = T_min
    Temp_grad[n-divider:n] = T_min

    if chunks%2 == 1:
        #odd
        Temp_grad[int(chunks/2)*divider:int(chunks/2+1)*divider] = T_max
    if chunks%2 == 0:
        #even
        Temp_grad[int(chunks/2-1)*divider:int(chunks/2+1)*divider] = T_max
    
    for i in range(move):
        step = (T_max-T_min)/move
        
        a += step
        
        Temp_grad[int(divider+i)] = a
        Temp_grad[-int(divider+i)-1] = a
        
        
    return(Temp_grad)


def init_ising_lattice(n,m):
    lattice = np.zeros((n,m),dtype=int)
    options = [-1,1]
    for i in range(n):
        for j in range(m):
            lattice[i,j] = choice(options)
    return lattice

def energydiff(S0,Sn,J,H):
    return 2*S0*(H+J*Sn)



def ising(isShow,n=200, m=50,nsteps=500000,H=0,J=1,T_min = 1,T_max = 4,chunks = 5, name = 'dir/temp', save = False):
    lattice = init_ising_lattice(n,m)
    Gradient = make_gradient(m,T_min,T_max,chunks)
    
    energy = 0
    energies = []
    spins = []
    spin = np.sum(lattice)
    pics = []
    filenames = []

    for step in range(nsteps):
        i = randrange(n)
        j = randrange(m)

        Sn = lattice[(i-1)%n,j]+lattice[(i+1)%n,j]+\
             lattice[i,(j-1)%m]+lattice[i,(j+1)%m]

        dE = energydiff(lattice[i,j],Sn,J,H)

        if dE < 0 or random() < exp(-dE/Gradient[j]):
            lattice[i,j] = -lattice[i,j]
            energy += dE
            energies.append(energy)
            spin += 2*lattice[i,j]

        spins.append(spin)

        if isShow:
            if step%5000 == 0:
                im = plt.imshow(lattice, animated=True)
                pics.append([im])
                if save == True:
                    filename = name + str(step) + '.png'
                    plt.savefig(filename)
                    filenames.append(filename)
            if step == nsteps-1:
                return pics, Gradient, filenames

n = 50
m=200
nsteps = 50000
H = 0
J = 1.0
T_min = 2.0
T_max = 3.5
chunks = 10
#filename = 'run1'
pics, Gradient, filenames = ising(True, n,m,nsteps,H,J,T_min,T_max,chunks,'Pics/run',True)

# spins = np.asarray(spins)/n**2

# pics = [plt.imshow([[1,1],[1,0]], animated=True), plt.imshow([[0,0],[0,1]], animated=True)]

fig = plt.figure()

# ani = animation.ArtistAnimation(fig, pics, interval=100, blit=True,
#                                 repeat_delay=0)
ani = animation.ArtistAnimation(fig, pics, interval=200, repeat_delay=1000,
                                blit=True)
plt.plot(Gradient)

# READ IMAGES AND CREATE GIF
images = []
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('animations/tempgrad1.gif', images)

# REMOVE FILES OF JUST PICS (POINTLESS AND TO AVOID CLUTTER)
for filename in filenames:
    os.remove(filename)
print("File Removed!")

# plt.plot(spins)
plt.show()
