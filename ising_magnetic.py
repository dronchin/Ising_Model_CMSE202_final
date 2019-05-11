from math import exp
from random import choice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from copy import deepcopy

def init_ising_lattice(n):
    '''Initializes a nxn grid with random spins and calculates the initial
    energy and magnetization
    Returns: lattice, initial energy, initial magnetization'''
    lattice = np.zeros((n,n))
    options = [-1.0,1.0]
    for i in range(n):
        for j in range(n):
            lattice[i,j] = choice(options)
    energy = findenergy(lattice)
    mag = findmag(lattice)
    return lattice, energy, mag

def findenergy(lattice):
    ''' finds the energy of the lattice'''
    E = 0
    for i in range(len(lattice)):
        for j in range(len(lattice[0])):
            Sn = lattice[(i+1)%n,j]+lattice[i,(j+1)%n]
            E -= lattice[i][j] *Sn
    return E

def findmag(s):
    ''' finds the magnetization of the lattice'''
    m = 0
    for i in range(len(s)):
        for j in range(len(s[0])):
            m += s[i][j]
    return m

def energydiff(S0,Sn):
    return 2*S0*Sn

def ising(n=20,mcc=1000,T=1):
    '''Create the ising model, given the lattice size and number of monte carlo
    cycles. The function return a list of updated energy and their average
    value.
    n,n2 are ranges of lattice, mcc is the number of rounds we are
    going to update, T is for Temperature.'''
    lattice, energy, mag = init_ising_lattice(n)
    E = 0
    M = 0
    M_abs = 0
    M2 = 0
    E2 = 0

    #precalulates the energy probabilities
    w = np.zeros(17)
    for de in range(-8,9,4): #include +8
        w[de+8] = np.exp(-de/T)

    #run through the experiments
    for mcstep in range(mcc):
        #n^2 steps is 1 experiment
        for step in range(n**2):
            i = np.random.randint(0,n)
            j =  np.random.randint(0,n)

            #periodic boundary condition. Edge of grid interacts with the '
            #opposite edge
            Sn = lattice[(i-1)%n,j]+lattice[(i+1)%n,j]+\
                 lattice[i,(j-1)%n]+lattice[i,(j+1)%n]

            dE = energydiff(lattice[i,j],Sn)
            if np.random.random() < w[int(dE)+8]: #Metropolis algorithm
                lattice[i,j] = -lattice[i,j]
                energy += dE
                mag += 2*lattice[i,j]
        # Averages the values over all experiments
        E += energy
        M += mag
        M_abs += np.abs(mag)
        M2 += np.abs(mag)**2
        E2 += np.abs(energy)**2

    #scales the values
    E /= float(mcc)
    M /= float(mcc)
    M_abs /= float(mcc)
    E2 /= float(mcc)
    M2 /= float(mcc)
    #calulates heat capacity and Susceptibility
    Cv = (E2 - E**2)/(float(n)**2*T**2)
    Sus = (M2 - M_abs**2)/(float(n)**2*T)
    E /= float(n)**2
    M /= float(n)**2
    M_abs /= float(n**2)
    return E, M, M_abs, Cv, Sus

n = 10
# mcc = 1
step = 0.05
T_vals = np.arange(2,2.6+step,step)

T_c = 2.269
mcc = 10000
energies = []
magnetics = []
Cvs = []
Suss = []
for T in T_vals: #runs over a range of temperatures
    print(n,T)
    start = time.time()
    energy, magnetic,magnetic_abs, Cv, Sus = ising(n,mcc,T)
    stop = time.time()
    print(stop-start)
    energies.append(energy)
    magnetics.append(magnetic_abs)
    Cvs.append(Cv)
    Suss.append(Sus)

#saves the arrays to be run by ShowGraph.py
np.save("./SavedVar/E/Energy_{}".format(n),energies)
np.save("./SavedVar/M/Magnetic_{}".format(n),magnetics)
np.save("./SavedVar/Cv/Cv_{}".format(n),Cvs)
np.save("./SavedVar/Sus/Sus_{}".format(n),Suss)

plt.plot(T_vals, energies)
plt.title("Energy vs Temperature")
plt.xlabel("Temperature")
plt.ylabel("Energy")
plt.savefig("./Pics/Energy.png")
plt.show()

plt.plot(T_vals, magnetics)
plt.title("Magnetization vs Temerature")
plt.xlabel("Temperature")
plt.ylabel("Magnetization")
plt.savefig("./Pics/Magnetization.png")
plt.show()

plt.plot(T_vals, Cvs)
plt.title("Heat Capacity vs Temperature")
plt.xlabel("Temperature")
plt.ylabel("Heat Capacity")
plt.savefig("./Pics/HeatCapacity.png")
plt.show()

plt.plot(T_vals, Suss)
plt.title("Susceptibility vs Temperature")
plt.xlabel("Temperature")
plt.ylabel("Susceptibility")
plt.savefig("./Pics/Susceptibility.png")
plt.show()
# plt.imshow(lattice)
# plt.show()







#
