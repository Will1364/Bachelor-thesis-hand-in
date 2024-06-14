## 
## antennaDirectivity.py
##
## Code created by: 
## William Hedegaard Langvad  (s214512)
## Nikolaj Schejbel Nielsen   (s214494)
##
## The purpose of this script it to calculate the antenna patterns 
## 

#%% Import of libraries used
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from scipy.interpolate import interp1d

# Sets directory as python file location
# All filepaths should be relative to this 
import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


#------------------------------------------------------------------------------
#%% Coordinate and field setup
res   = 1000
theta = np.linspace(0.0000000000001, np.pi, num = res)         
phi   = np.linspace(0.0000000000001, 2 * np.pi, num = res)     


Theta, Phi = np.meshgrid(theta, phi)
dtheta = np.abs(theta[0] - theta[1])
dphi   = np.abs(phi[0] - phi[1])

# Electromagnetic parameters 
E0    = 1
eta_0 = 4 * np.pi * 29.97924 

# Antenna size
a = 0.02286
b = 0.01016

# Wavelengths
lmbda = 0.03
k     = 2*np.pi / lmbda

# Coordinates for calculating coefficients
D = 0.06
n = 10000
m = np.linspace(-lmbda/(D), lmbda/(D), num = n * 3) * 3 / 4


#------------------------------------------------------------------------------
#%% Antenna far-field pattern function
X  = k * a / 2 * np.sin(Theta) * np.cos(Phi)
Y  = k * b / 2 * np.sin(Theta) * np.sin(Phi)
Cf = a * b * k * E0 /(2 * np.pi)

F_theta = Cf/2 * np.sin(Phi) * (1 + np.cos(Theta)) * np.sin(X)/X * np.sin(Y)/Y
F_phi   = Cf/2 * np.cos(Phi) * (1 + np.cos(Theta)) * np.sin(X)/X * np.sin(Y)/Y

F = np.abs(np.sqrt(F_theta**2 + F_phi**2))

#------------------------------------------------------------------------------
#%% Directivity, found with numerical integration
tmp =  np.sum(F * np.sin(Theta), axis = 1) * dtheta
tmp =  np.sum(tmp) * dphi
D = 4 * np.pi * F / tmp

# Max directivity
D_0 = 4 * np.pi * (a * b / lmbda**2)


#------------------------------------------------------------------------------
#%% Array correction parameters 
# Vertical correction is theta integration over target / integration over whole theta
target_size = 2 * np.arctan(1/(2*np.sqrt(2))/9)
mask        = np.abs(theta) < target_size/2

D_ttarget   = np.sum(D[:, mask], axis = 1) * dtheta

mask        = theta < np.pi/2
D_ttotal    = np.sum(D, axis = 1) * dtheta

n = int(res/4)

v_coeff = D_ttarget/D_ttotal

# H coefficient is just directivity at theta = 0
h_coeff = D[n,:]


#------------------------------------------------------------------------------
#%% Interpolate of coefficients within reconstructed range
idx     = np.sum(np.abs(phi) < np.max(m)*1.1)
v_coeff = v_coeff[:idx]
v_coeff = np.append(np.flip(v_coeff),v_coeff)
v_coord = phi[:idx]
v_coord = np.append(-np.flip(v_coord),v_coord)

f = interp1d(v_coord, v_coeff)
vD = f(m)


idx     = np.sum(np.abs(theta) < np.max(m)*1.1)
h_coeff = h_coeff[:idx]
h_coeff = np.append( np.flip(h_coeff),h_coeff)
h_coord = theta[:idx]
h_coord = np.append(-np.flip(h_coord),h_coord)

f  = interp1d(h_coord, h_coeff)
hD = f(m)


#------------------------------------------------------------------------------
#%% Figures antenna correction coefficient

# Calculating the antenna coefficient
G_0   = 10**(6/10)              # Gain coefficient, ~6dB
hD    = hD * G_0 / np.max(hD)   # Scale directivity to gain
coeff = 1/(hD * vD)

if True:
    fig, ax = plt.subplots()
    
    plt.plot(m * 180/np.pi,1/(hD * vD), c = '#069AF3', label = 'Antenna correction coefficient')
    # plt.scatter(m * 180/np.pi, np.array([0, 0]), color='black', marker = 'x', zorder = 3)
    
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Correction factor')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    # plt.xlim(-fov[1] * 1.1 * 180/np.pi, fov[1] * 1.1 * 180/np.pi)
    # plt.ylim(bottom = -0.3 * np.max(np.real(I * coeff)))
  
    plt.show()
    plt.close()




#------------------------------------------------------------------------------
#%% Figures with 2D cuts

if True:
    # Vertical cut
    n = int(res/2)
    
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})                                                                    
    ax.plot(theta, (np.real(D[n,:])/np.max(D[n,:]))**2, c = 'blue')   
    ax.plot(-theta, (np.real(D[n,:])/np.max(D[n,:]))**2,  c = 'blue')                                              # for hvert nyt plot skriver vi plt.figure(n) hvor n er det nummer som det nye plot vil være i rækken
    ax.set_rmax(1)
    ax.set_rticks([])  # Less radial ticks
    ax.grid(True)
    ax.set_theta_zero_location('N')
    ax.set_thetalim([-np.pi, np.pi])
    ax.set_rgrids([0.2, 0.4, 0.6, 0.8, 1])
    # plt.legend(loc = 'lower left')
    plt.savefig('figs/temp.svg')
    
    
    
if True:
    # Horizontal cut
    n = int(res/4)
    
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})                                                                    
    ax.plot(theta, (np.real(D[n,:])/np.max(D[n,:]))**2, c = 'blue')   
    ax.plot(-theta, (np.real(D[n,:])/np.max(D[n,:]))**2,  c = 'blue')                                              # for hvert nyt plot skriver vi plt.figure(n) hvor n er det nummer som det nye plot vil være i rækken
    ax.set_rmax(1)
    ax.set_rticks([])  # Less radial ticks
    ax.grid(True)
    ax.set_theta_zero_location('N')
    ax.set_thetalim([-np.pi, np.pi])
    ax.set_rgrids([0.2, 0.4, 0.6, 0.8, 1])
    # plt.legend(loc = 'lower left')
    plt.savefig('figs/temp.svg')
    

