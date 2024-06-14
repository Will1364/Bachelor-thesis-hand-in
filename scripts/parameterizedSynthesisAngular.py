## 
## parameterizedSynthesisAngular.py
##
## Code created by: 
## William Hedegaard Langvad  (s214512)
## Nikolaj Schejbel Nielsen   (s214494)
##
## The purpose of this program is to simulate visibilites and 
## reconstructed images for the near field setup.  
##

#%% Import of important libraries
import matplotlib.pyplot as plt
import numpy as np


# Sets directory as python file location
# All filepaths should be relative to this 
import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

plt.rcParams['figure.dpi'] = 300

# Import custom function
from synthesisFunctions import visibilitySim, VisibilitySimAngular, IntensityReconstructionAngular, hogbomCLEAN

#%% Additional functions
def placePointSource(_ang, _d, _Sf, _intensity):
    # Places source at angle _ang in source field _Sf, with _intensity
    _ang   = _ang * np.pi/180 
    _theta = np.arctan(Sf[0, :]/d) 
    _idx   = np.sum(_theta < _ang)
    
    _Sf[2,  _idx] = _intensity
    return
    
#------------------------------------------------------------------------------
#%% Setup of dimensions used in the simulations

d       = 9.0        # Distance to target
D       = 0.06       # Mindste afstand mellem to målere [m]
lmbda   = 0.03       # Bølgelængde [m]
N       =   16       # Number of measurements
w       = 2          # width of the target, [m]
n       = 10000      # amount of points in field

# Setting up antenna coordinates and saving into antenna_coord array
antenna_coord = np.zeros([2, N])  
D_array = np.linspace(0, D * (N-1), num = N) + 0.5 * D  
for i in range(N):
    x = D_array[i]
    y = d
    antenna_coord[:, i] = [x, y]
    
# Source field given in x, y, source intensity
Sf        = np.zeros([3, n])  
Sf[0, :]  = np.linspace(-w, w, n) * 0.5
Sf[1, :]  = 0
Sf[2, :] = 0

# Point source output
#Sf[2, 5000] = 1


# Extended sim used
Sf[2,    1:2000] = 1
Sf[2, 3500:5000] = 1
Sf[2, 6500:7000] = 1
Sf[2, 9750:9999] = 1



# Multiple sources output
#placePointSource(-5, 9, Sf, 1)
#placePointSource(-2, 9, Sf, 0.5)
#placePointSource( 3, 9, Sf, 1)

# placePointSource(0, 9, Sf, 1)

# Source filed used for emulating actual target
if False:
    # Model of Experiment source field
    Sf        = np.zeros([3, n])  # Source field given in x, y, source intensity
    Sf[0, :]  = np.linspace(-w, w, n) * 0.5
    Sf[1, :]  = 0
    Sf[2, :] = 10

    Sf[2,    0:375]  = 10 * (1 - 0.1) + (0.1 * (273 + 15))
    Sf[2,  375:1125] = 273 + 15
    Sf[2, 1125:1300] = 10 * (1 - 0.1) + (0.1 * (273 + 15))  # Width of roughly 15 cm
    Sf[2, 1300:2050] = 273 + 15
    Sf[2, 2050:2400] = 10 * (1 - 0.1) + (0.1 * (273 + 15))
    Sf[2, 2400:3150] = 273 + 15
    Sf[2, 3150:3850] = 10 * (1 - 0.1) + (0.1 * (273 + 15))
    Sf[2, 3850:4600] = 273 + 15
    Sf[2, 4600:6000] = 10 * (1 - 0.1) + (0.1 * (273 + 15))
    Sf[2, 6000:6750] = 273 + 15
    Sf[2, 6750:8850] = 10 * (1 - 0.1) + (0.1 * (273 + 15))
    Sf[2, 8850:9600] = 273 + 15
    Sf[2, 9600:]     = 10 * (1 - 0.1) + (0.1 * (273 + 15))


# Noise field used to emulate noise, same setup as Sf
dNf  =  50
wNf  = 400
Nf   = np.zeros([3, n])
Nf[0, :]  =  np.linspace(-wNf, wNf, n)
Nf[1, :]  = -dNf
Nf[2, :]  =  np.random.rand(n) * 0

# Field of view
fov   = np.array([-lmbda/(D), lmbda/(D)]) / 4

# Angular coordinates
theta = np.linspace(-lmbda/(D), lmbda/(D), num = n)/4   * 3
theta_t = np.linspace(-np.pi/2, np.pi/2, num = n)



#------------------------------------------------------------------------------
#%% Simulating visibilites and intensity
# Zero bandwidth
Vu  = visibilitySim(lmbda, antenna_coord, Sf, Nf, _bw = 0)
I   = IntensityReconstructionAngular(lmbda, theta, antenna_coord, Vu)

# 2 GHz bandwidth
Vu_02G = visibilitySim(lmbda, antenna_coord, Sf, Nf, _bw = 02e+09)
I_02G  = IntensityReconstructionAngular(lmbda, theta, antenna_coord, Vu_02G)

#------------------------------------------------------------------------------
#%% Figures used in the report
if True:
    # Visibility function
    fig, ax = plt.subplots()
    plt.xlabel('Antenna distance [wavelengths]')  

    ax.plot(2 * antenna_coord[0,:]/lmbda, np.real(Vu)/(np.max(np.abs(Vu))), marker = 'x', c = '#069AF3', label = 'Real(Vu)')
    secax  = ax.twinx()
    secax.plot(2 * antenna_coord[0,:]/lmbda, np.angle(Vu), marker = 'x', c = 'gray', label = 'Phase', ls = '--')
    
    ax.grid(True, axis = 'y')    
    ax.set_ylabel('Real part of visibility [Normalized]', c =  '#069AF3')
    secax.set_ylabel('Phase of visibility [rad]', c =  'gray')
    
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = secax.get_legend_handles_labels()

    plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)

    plt.show()
    plt.close()


if True:
    # Plot of the image
    fig, ax = plt.subplots()
    
    plt.plot(theta * 180/np.pi, (np.real(I)) /(1 + 0* np.max(np.abs(np.real(I)))), c = '#069AF3', label = 'Image, BW = 0 Hz')
    
    plt.plot(np.arctan(Sf[0,:]/d) * 180/np.pi, Sf[2,:] / np.max(Sf[2,:]), c = 'gray', ls = '--', label = 'Source distribution')
    plt.scatter(fov * 180/np.pi, np.array([0, 0]), color='black', marker = 'x', zorder = 3)
    
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Normalised intensity')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-fov[1] * 1.1 * 180/np.pi, fov[1] * 1.1 * 180/np.pi)
    plt.ylim(bottom = -0.3)
  
    plt.show()
    plt.close()
    
  
#------------------------------------------------------------------------------
#%% CLEAN image and figure
clean_im, x_coord = hogbomCLEAN(lmbda, theta, fov, antenna_coord, I, Sf)

 
if True:
    # Plot of the CLEAN image
    fig, ax = plt.subplots()
    
    plt.plot(x_coord * 180/np.pi, (np.real(clean_im)) / np.max(np.abs(np.real(clean_im))), c = '#069AF3', label = 'CLEAN image')
    # plt.plot(x_coord * 180/np.pi, (np.real(clean_im_02G)) / np.max(np.abs(np.real(clean_im))), c = '#E50000', label = 'Image, BW = 2 GHz', ls = '-')
    
    plt.plot(np.arctan(Sf[0,:]/d) * 180/np.pi, Sf[2,:] / np.max(Sf[2,:]), c = 'gray', ls = '--', label = 'Source distribution')
    plt.scatter(fov * 180/np.pi, np.array([0, 0]), color='black', marker = 'x', zorder = 3)
    
    # plt.title('Multiple sources image')
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Normalised CLEAN intensity')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    lim = fov[1] * 1.1 * 180/np.pi
    plt.xlim(-lim, lim)
    plt.ylim(bottom = -0.3)
    
    plt.show()
    plt.close()

#------------------------------------------------------------------------------
#%% Figure with errors on distance between antennas
if False:
    
    itt = 20                    # Number of itterations
    
    fig, ax = plt.subplots()
    
    # Plot intensity wihtout errors on antenna spacing
    plt.plot(theta * 180/np.pi, (np.real(I)) / np.max(np.abs(np.real(I))), c = '#069AF3', label = 'Image', zorder = 3)
    
    print("Creating noise image")
    for i in range(itt):
        print("\rItteration {:>4}/{:>4}".format(i+1, itt), end = '')
        
        # Add noise to coordinate
        ns_antenna_coord = antenna_coord
        ns_antenna_coord[0, :] = ns_antenna_coord[0, :] + np.random.normal(0, 0.001, size = N)
        
        # Calculate new Vu and intensity
        ns_Vu  = visibilitySim(lmbda, ns_antenna_coord, Sf, Nf, _bw = 0)
        ns_I   = IntensityReconstructionAngular(lmbda, theta, antenna_coord, ns_Vu)
        
        # Plot the new intensities
        if i == 0:
            plt.plot(theta * 180/np.pi, (np.real(ns_I)) / np.max(np.abs(np.real(I))), c = '#A9A9A9', label = 'Noise images')
        else:
            plt.plot(theta * 180/np.pi, (np.real(ns_I)) / np.max(np.abs(np.real(I))), c = '#A9A9A9')
    
    plt.plot(np.arctan(Sf[0,:]/d) * 180/np.pi, Sf[2,:] / np.max(Sf[2,:]), c = 'gray', ls = '--', label = 'Source distribution', zorder = 2)
    plt.scatter(fov * 180/np.pi, np.array([0, 0]), color='black', marker = 'x', zorder = 3)
    
    # plt.title('Multiple sources image')
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Normalised intensity')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-fov[1] * 1.1 * 180/np.pi, fov[1] * 1.1 * 180/np.pi)
    plt.ylim(bottom = -0.3)

    plt.show()
    plt.close()


    
#------------------------------------------------------------------------------
#%% Figure with constant sky intensity 
# Sky with constant intensity 
sky   = np.zeros(np.size(theta_t)) + np.random.normal(100, 0, size = np.size(theta_t))
    
if True:
    fig, ax = plt.subplots()
    
    # Calculate visibility and intensity
    ns_Vu  = VisibilitySimAngular(lmbda, 0, antenna_coord, sky, theta_t)
    ns_I   = IntensityReconstructionAngular(lmbda, theta_t, antenna_coord, ns_Vu)
    
    
    plt.plot(theta_t * 180/np.pi, (np.real(ns_I)) / np.max(np.abs(np.real(ns_I))), c = '#069AF3', label = 'Image')
    plt.scatter(fov *  180/np.pi, np.array([0, 0]), color='black', marker = 'x', zorder = 3)
    
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Normalised intensity')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-fov[1] * 1.1 * 180/np.pi, fov[1] * 1.1 * 180/np.pi)
    plt.ylim(bottom = -0.3)

    plt.show()
    plt.close()