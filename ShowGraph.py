import numpy as np
import matplotlib.pyplot as plt
import glob

step = 0.05
T_vals = np.arange(2,2.6+step,step)
<<<<<<< HEAD
n_vals = [10,20,40,60,80]
=======
n_vals = [10,20,40,60]
>>>>>>> 0b90a0f8f2cc878287790c8cde76f64268bc979c
T_c = 2.269

def graphit(file_list, name):
    '''When given a list of files paths for .npy files, it reads the file and
    plots the values on the same graph'''
    print(file_list)
    for i,file in enumerate(file_list):
        A = np.load(file)
        plt.plot(T_vals, A, label="{}".format(n_vals[i]))

    plt.axvline(x=T_c)
    plt.title("{} vs Temperature".format(name))
    plt.legend()
    plt.xlabel("Temperature")
    plt.ylabel(name)
    plt.savefig("./Pics/{}_multi.png".format(name))
    plt.show()

#finds all the .npy files from each folder in the ./SavedVar folder
file_Cv = glob.glob("./SavedVar/Cv/*.npy")
file_E = glob.glob("./SavedVar/E/*.npy")
file_M = glob.glob("./SavedVar/M/*.npy")
file_Sus = glob.glob("./SavedVar/Sus/*.npy")

graphit(file_Cv, "HeatCapacity")

graphit(file_E, "Energy")

graphit(file_M, "Magnetization")

graphit(file_Sus, "Susceptibility")

#
