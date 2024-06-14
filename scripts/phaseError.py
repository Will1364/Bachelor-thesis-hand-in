## 
## phaseError.py
##
## Code created by: 
## William Hedegaard Langvad  (s214512)
## Nikolaj Schejbel Nielsen   (s214494)
##
## The purpose of this script is to determine the phase 
## errors when in the near field
## 

#%% Import of important libraries 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import glob
import pandas as pd
from scipy.interpolate import CubicSpline, interp1d
from scipy import integrate

# Sets directory as python file location
# All filepaths should be relative to this 
import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

cmap_gray = cm.get_cmap('Greys')
cmap_blue = cm.get_cmap('winter')

#------------------------------------------------------------------------------
#%% Array parameters 
lmbda = 0.03         # Wavelength [m]
d       = 9          # Distance to target [m]
D       = 0.06       # Antenna spacing [m]
N       = 16         # Number of samples 

# Antenna array spacing 
D_array = np.linspace(0, D * (N-1), num = N) + 0.5 * D  

# Source coordinate range
x1      = np.linspace(-1.25, 1.25, 100)         

# Antenna spacings to color in figure
idx     = np.array([0, 6, 9, 11, 13, 15])

#------------------------------------------------------------------------------
#%% Calculates and plots phase errors
fig, ax = plt.subplots()
for i in range(np.size(D_array)):
    # Angle towards source (far field)
    theta = np.arctan(x1/d)

    # Distance to antenna 1 and 2 (near field)
    R1 = np.sqrt((D_array[i]/2+x1)**2 + (0-d)**2)
    R2 = np.sqrt((-D_array[i]/2+x1)**2 + (0-d)**2)

    # Phase difference between antenna
    dR_near = np.abs(R1 - R2)/lmbda
    dR_far  = np.abs(D_array[i]/lmbda * np.sin(theta))

    # Difference in phase difference between near field and far field
    ddiff  = np.abs(dR_far - dR_near)  * 2 * np.pi * 180/np.pi
    
    if np.any(idx == i):
        plt.plot(np.arctan(x1/d)*180/np.pi, ddiff, c = cmap_blue(1.5-(i+4)/(N)),label = "D = {:.2f}".format(2*D_array[i]), zorder = 3)
    else:
        plt.plot(np.arctan(x1/d)*180/np.pi, ddiff, c = cmap_gray((i+1)/(N*2)))
    
plt.grid(True)
plt.ylabel('Phase difference [degrees]')
plt.xlabel('Angular position [degrees]')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
