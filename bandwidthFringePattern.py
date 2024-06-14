## 
## bandwidthFringePattern.py
##
## Code created by: 
## William Hedegaard Langvad  (s214512)
## Nikolaj Schejbel Nielsen   (s214494)
##
## The purpose of this script it to calculate the fringe pattern and bandwidth effects 
## 


#%% Import of libraries used
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

plt.rcParams['figure.dpi'] = 300

cmap = cm.get_cmap('Greys')

#------------------------------------------------------------------------------
#%% Parameters for baseline and bandwidth setup
theta  = np.linspace(-np.pi/2, np.pi/2, num = 10000)
D      = 0.09                                            # Mindste afstand mellem to målere
dlmbda = 0.005
lmbda  = np.linspace(0.030 - dlmbda, 0.030 + dlmbda, num = 7)                                            # Bølgelængde 
N      = 23                                              # Antal målinger der laves (En måling for hver N * D/lmbda )



#------------------------------------------------------------------------------
#%% Calculating the fringe function
Fringe = np.zeros(np.size(theta)) 


fig, ax = plt.subplots()
# The for loop plots the fringe pattern for various wavelengths within the bandwidth
for i in range(np.size(lmbda)):                          
    n = np.size(lmbda)    
    Fringe_tmp = np.exp(2j*np.pi*D*np.sin(theta)/lmbda[i])
    
    if i == int(n/2):
        plt.plot(theta, np.real(Fringe_tmp), c = 'green', ls = '--', label = r'Central bandwidth')        
        print(i)
    else:
        plt.plot(theta, np.real(Fringe_tmp), c = cmap((i+2)/(n * 2)), ls = '--')
    

# The for loop calculates the bandwidth affected fringe pattern
dlmbda_fac           = (1/np.max(lmbda) - 1/np.min(lmbda))
bandwidth_factor = np.sinc(D*np.sin(theta)*dlmbda_fac) * np.exp(2j*np.pi*D*np.sin(theta)/0.03)

Fringe = np.zeros(np.size(theta)) 
lmbda = np.linspace(0.030 - dlmbda, 0.030 + dlmbda, num = 10**4)                                           
for i in range(np.size(lmbda)):                        
    n = np.size(lmbda)    
    Fringe_tmp = np.exp(2j*np.pi*D*np.sin(theta)/lmbda[i])
    Fringe = Fringe + Fringe_tmp
    
Fringe = Fringe/np.size(lmbda)  

# The rest of the code is figure parameters                          
plt.plot(theta, np.real(Fringe), c = 'red', ls = '-', label = r'Total')

plt.ylabel(r'Normalized amplitude $\:$')
plt.xlabel(r'Angular position [rad]')
plt.legend(loc='upper left')
plt.savefig('figs/temp.svg')