
"""
Group: Casey Chartier, Nicholas Dronchi, Joseph Slivka, Zhiyang Yu

Taken from author: prtkrm
    https://github.com/prtkm/ising-monte-carlo
"""

"""

Main Function: def ising(isShow, n=200,nsteps=500000,H=0,J=1,T=1,name='temp',save)
"""

import sys
from math import exp
from random import randrange,choice,random
from numpy import zeros
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import imageio
import os


def init_ising_lattice(n):
    lattice = zeros((n,n),dtype=int)
    options = [-1,1]
    for i in range(n):
        for j in range(n):
            lattice[i,j] = choice(options)
    return lattice

def energydiff(S0,Sn,J,H):
    return 2*S0*(H+J*Sn)


def ising(isShow, n=200,nsteps=500000,H=0,J=1,T=1, name='dir/temp',save=False):
    """
    Parameters:
        isShow: show graphs every 10000 iterations

        n: for an n x n grid of spin states

        nsteps: number of steps to run until completion

        H, J, T : ising model parameters for external magnetic field,
                    Interaction Energy, and Temperature respectively

        name : user-defined filename to where to save (include directory)
                The pics are automatically deleted after creating the gif

        save : default False, save the file or not to animations/
                user-defined parater may need to be changed in run

    Returns:
        pics : pictures to be animated in ArtistAnimation

        spins : array of sum of spins over time

        filenames : list of files to be made into a gif

    """

    lattice = init_ising_lattice(n)
    energy = 0
    energies = []
    spins = []
    spin = np.sum(lattice)
    pics = []
    filenames = []

    fig = plt.figure()

    for step in range(nsteps):
        i = randrange(n)
        j = randrange(n)

        Sn = lattice[(i-1)%n,j]+lattice[(i+1)%n,j]+\
             lattice[i,(j-1)%n]+lattice[i,(j+1)%n]

        dE = energydiff(lattice[i,j],Sn,J,H)

        if dE < 0 or random() < exp(-dE/T):
            lattice[i,j] = -lattice[i,j]
            energy += dE
            energies.append(energy)
            spin += 2*lattice[i,j]

        spins.append(spin)

        if isShow:

            if step%10000 == 0:

                im = plt.imshow(lattice, animated=True)
                pics.append(im)
                if save == True:
                    filename = name + str(step) + '.png'
                    plt.savefig(filename)
                    filenames.append(filename)
            if step == nsteps-1:
                return pics, spins, filenames

"""
MAIN PROGRAM

Returns:
    animation of the evolution of the ising model, saved to user defined
    gif_save_file

"""

#SET PARAMS
n = 100
nsteps = 100000
H = 0
J = 1.0
T = 1.0
gif_save_file = 'animations/temp_001.gif'

#RUN ISING MODEL
pics, spins, filenames = ising(True, n,nsteps,H,J,T,'Pics/temp',True)

spins = np.asarray(spins)/n**2

#CREATE FIGURE
fig = plt.figure()

#ANIMATE FIGURE
ani = animation.ArtistAnimation(fig, pics, interval=50, repeat_delay = 0,
                            blit = True)

# READ IMAGES AND CREATE GIF
images = []
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave(gif_save_file, images)
print('GIF created as ' + gif_save_file)


# REMOVE FILES OF JUST PICS (POINTLESS AND TO AVOID CLUTTER)
for filename in filenames:
    os.remove(filename)
print("Files Removed!")
