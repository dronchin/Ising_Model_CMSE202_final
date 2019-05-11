# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 20:43:52 2018

@author: Nicolas Dronchi
"""
from math import exp
from random import randrange,choice,random
import numpy as np
import matplotlib.pyplot as plt
import time

def energydiff(S0,Sn):
    # S0: state of grid at that spot
    # Sn: neighbor configuration
    return 2*S0*Sn

def init_ising_lattice(n):
    lattice = np.ones((n,n))
    options = [-1.0,1.0]
    for i in range(n):
        for j in range(n):
            lattice[i,j] = choice(options)
    energy = findenergy2x2(lattice)
    mag = findmag2x2(lattice)
    return lattice, energy, mag

def findenergy2x2(lattice):
    E = 0
    for i in range(len(lattice)):
        for j in range(len(lattice[0])):
            Sn = lattice[(i+1)%n,j]+lattice[i,(j+1)%n]
            E -= lattice[i][j] *Sn
    return E

def findmag2x2(s):
    m = 0
    for i in range(len(s)):
        for j in range(len(s[0])):
            m += s[i][j]
    return m

def findCvSus(T):
    Z = 12.0 + 2.0*np.exp(-8.0/T) + 2.0*np.exp(8.0/T)
    E = -8.0*np.sinh(8/T)/(np.cosh(8/T)+3)
    Cv = (1/T**2)*(256*np.cosh(8/T)/Z - E**2)
    M = (2*np.exp(8/T)+4)/(np.cosh(8/T)+3)
    Sus = (1/T)*((8*np.exp(8/T)+8)/(np.cosh(8/T)+3)- M**2)
    return E, M, Cv, Sus

def ising(n=2,mcc=1,T=1):
    lattice, energy, mag = init_ising_lattice(n)
    E = 0
    M = 0
    M_abs = 0
    M2 = 0
    E2 = 0
    for i in range(mcc):
        for step in range(n**2):
            i = randrange(n)
            j = randrange(n)

            Sn = lattice[(i-1)%n,j]+lattice[(i+1)%n,j]+\
                 lattice[i,(j-1)%n]+lattice[i,(j+1)%n]

            dE = energydiff(lattice[i,j],Sn)

            if dE < 0 or random() < exp(-dE/T):
                # Boltzmann distribution to flip or not flip lattice spin
                # print("flipped")
                lattice[i,j] = -lattice[i,j]
                energy += dE
                mag += 2*lattice[i,j]
        E += energy
        M += mag
        M_abs += np.abs(mag)
        M2 += np.abs(mag)**2
        E2 += energy**2

    # keeping track of values
    E /= float(mcc)
    M /= float(mcc)
    M_abs /= float(mcc)
    E2 /= float(mcc)
    M2 /= float(mcc)

    Cv = (E2 - np.abs(E)**2)/(T**2)
    Sus = (M2 - np.abs(M_abs)**2)/(T)

    return E, M_abs, Cv, Sus

n = 2
mcc = 100000
T_vals = np.linspace(1,4,num=10)
energies = []
magnetics = []
Cvs = []
Suss = []
for T in T_vals:
    energy, magnetic, Cv, Sus = ising(n,mcc,T)
    print("Numerical", energy, magnetic, Cv, Sus)
    energy_exp, magnetic_exp, Cv_exp, Sus_exp  = findCvSus(T)
    print("Analytical", energy_exp, magnetic_exp, Cv_exp, Sus_exp)
    dE = np.abs(energy - energy_exp)/energy_exp
    dM = np.abs(magnetic - magnetic_exp)/magnetic_exp
    dCv = np.abs(Cv - Cv_exp)/Cv_exp
    dSus = np.abs(Sus - Sus_exp)/Sus_exp

    energies.append(dE)
    magnetics.append(dM)
    Cvs.append(dCv)
    Suss.append(dSus)
    # energies.append(energy)
    # magnetics.append(magnetic)
    # Cvs.append(Cv)
    # Suss.append(Sus)

plt.plot(T_vals, energies)
plt.plot(T_vals, magnetics)
plt.plot(T_vals, Cvs)
plt.plot(T_vals, Suss)
plt.show()

T_c = 2.269

#start = time.time()
#mags = []
#suscepts = []
#for T in T_vals:
#    print(T)
#    plt.imshow(lattice)
#    plt.show()
    # mag = np.abs(mag/nsteps)
    # print(mag)
#    mags.append(avgMag)
#    suscepts.append(avgSuscept)

#T_vals = np.asarray(T_vals)
#mags = np.asarray(mags)
#suscepts = np.asarray(suscepts)
#print(suscepts)
#print(np.abs(mags)**2)
#Susceptibility = (suscepts - np.abs(mags)**2)/(T_vals)
##print("Susceptibility", Susceptibility)
#stop = time.time()
#total = stop - start
#print("total time", total)
#
#mags = np.abs(mags/mags.max())
#plt.plot(T_vals, mags)
#plt.show()
#
#suscepts = suscepts/suscepts.max()
#plt.plot(T_vals, suscepts)
#plt.show()
