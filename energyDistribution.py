from math import exp
from random import randrange,choice,random
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

def ising(n=20,n2=400,mcc=1000,T=1):
    '''Slightly modified version of the ising model found in ising_magnetic only
    calculates only the energy'''
    lattice, energy, mag = init_ising_lattice(n)
    E = 0
    E_list = []
    E2 = 0
    avgE = []

    for mcstep in range(mcc):
        for step in range(n2):
            i = randrange(n)
            j = randrange(n)

            Sn = lattice[(i-1)%n,j]+lattice[(i+1)%n,j]+\
                 lattice[i,(j-1)%n]+lattice[i,(j+1)%n]

            dE = energydiff(lattice[i,j],Sn)

            if dE < 0 or random() < exp(-dE/T):
                lattice[i,j] = -lattice[i,j]
                energy += dE
        E += energy
        E2 += energy**2
        E_list.append(energy)
        avgE.append(E2/(mcstep+1))
    return E_list, avgE

n = 20
mcc = 1000
T_vals = [0.01,2,2.269,3,4000] #np.linspace(1,4,num=5)
n_vals = [10,20,50,100]
T_c = 2.269

all_E = []
for T in T_vals:
    print(n,T)
    E_list, avgE = ising(n,n**2,mcc,T)
    all_E.append(E_list)
all_avgE = []
T = 1
# for n in n_vals:
#     print(n,T)
#     E_list, avgE = ising(n,100,mcc,T)
#     all_avgE.append(avgE)

for i in range(len(all_E)):
    plt.hist(all_E[-i-1], bins=np.linspace(-800,100,50), label="{}".format(T_vals[-i-1]))
plt.legend()
plt.title("Energy distributions at different temperatures")
plt.xlabel("Energy")
plt.ylabel("Fequency")
plt.savefig("./Pics/EnergyDistribution.png")
plt.show()

for i,E in enumerate(all_avgE):
    plt.plot(E, label="{}".format(n_vals[i]))
plt.legend()
plt.title("Energy Thermalization at different temperatures")
plt.xlabel("MC cycles")
plt.ylabel("Energy")
# plt.savefig("./Pics/EnergyThermalization.png")
plt.show()
#
