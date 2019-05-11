# final_project_cmse202
CMSE 202 Final Project

# Project Name

### Particle Interaction --- Ising Model Based On Monte Carlo Metrapolis Method 


# Group Members

### Joseph Slivka, Nicolas Dronchi, Casey Chartier, Zhiyang Yu

# Questions We Want To Answer Through This Project 

We want to model a microscopic particles interaction. The interaction is not simple collisions but the phase transformation of a material. 

## Getting Started

### Prerequisites

```
Python 3.0
!pip install imageio --user
```

### Packages Installed

```
from math import exp
from random import choice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from copy import deepcopy

```

## What is Ising Model; What is Monte Carlo Metrapolis 

Ising Model: A lattice in 2d, 3d or 1d statement. For each grid point there is one particle and there are two spinning statement(up or down). 
For the interactions between particles we only consider the nearest four grid points: up, low, left and right ones. (2D)

Monte Carlo is a stochastic simulation which can get an approximate value for a question through computer modeling. 

Metrapolis rules a probablity for particles to have a balance spinning statement which is given by e^(-delta(T)/(kT). 
k is the Boltzzman Constant 



## Physical Constants

Energy (E), Temperature (T), Boltzmann Constant (k) , coupling constant (J), magnetization (M)


## Functions 



### init_ising_lattice(n) 

Create an initial lattice with 0 energy, magnetization (intensity) and spin. Range of spins is from -1 to 1. Return the 3 physical 
quantities. 

```
def init_ising_lattice(n):
    lattice = np.zeros((n,n))
    options = [-1.0,1.0]
    for i in range(n):
        for j in range(n):
            lattice[i,j] = choice(options)
    energy = findenergy(lattice)
    mag = findmag(lattice)
    return lattice, energy, mag
```

### findenergy(lattice)

Get the energy of a choosed particle, which is defined by -J of one particle times sum of spin of four particles around it: 
-J(i,j)* [ spin(i+1,j) + spin(i,j+1) + spin(i-1,j) + spin(i,j-1) ]   Hint: J is the coupling constant which is depends on materials

```
def findenergy(lattice):
    E = 0
    for i in range(len(lattice)):
        for j in range(len(lattice[0])):
            Sn = lattice[(i+1)%n,j]+lattice[i,(j+1)%n]
            E -= lattice[i][j] *Sn
    return E
```

### findmag(s) 

Get the magnetization (intensity) of a choosed particle, which is defined by sum of the spin(i) devided by square of N  
Hint: N is the length of the lattice.

```
def findmag(s):
    m = 0
    for i in range(len(s)):
        for j in range(len(s[0])):
            m += s[i][j]
    return m
```

### energydiff(S0,Sn)

Get the energy difference after a switch a spin, which is defined by 2 times S0 times Sn
Hint: S0 is the state of grid , Sn is the neighbor configuration 

```
def energydiff(S0,Sn):
    # S0: state of grid at that spot
    # Sn: neighbor configuration
    return 2*S0*Sn
```


### ising(n,n2,mcc,T)

Create the ising model, given the lattice range, initial energy and magnetization (intensity). Then choose a random particle(i,j).
if Energy is positive(delta E < 0), switch the spin (spin' = -spin); if Energy is negative(delta E >= 0), the particle has 
a probability P = 1/(A*e^(E/(kT))) to switch. The function return a list of updated energy and their average value. 
Hint: n,n2 are ranges of lattice, mcc is the number of rounds we are going to update, T is for Temperature. 

```
def ising(n=20,n2=400,mcc=1000,T=1):
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
```

## Plotting, Visulization and Animation
For plots, we compared the temperature with three physical quantities energy, magnetization (intensity) and heat capacity. 
If we keep updating the range of lattice (n) in these plots, we can get the approximate temperature for three physical quantities 
reach zero.

We also create bar charts to see the energy distribution at different temperatures. 

For Animation, we give the model a range of temperature and keep updating and tracing the state of ising model. 



## Further Expectation to this Project

Consider a 3D ising model. There will be more interactions between particles in this model. 



## Acknowledgments
*Magnus Holter-Sørensen Dahle(Aug. 2014) The critical temperature of the two-dimensional Ising Model, applying the gradient method. 
              
              http://folk.ntnu.no/magnud/Projects/NumFys_exam/Home_Exam_Computational_Physics_TFY4235.pdf
              
*prtkm from Github(Dec. 2014) Monte-Carlo Simulations of the 2-D Ising Model. 
              
              https://github.com/prtkm/ising-monte-carlo/blob/master/ising-monte-carlo.org 
              
