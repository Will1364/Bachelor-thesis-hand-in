## 
## simulatedSynthesis.py
##
## Code created by: 
## William Hedegaard Langvad  (s214512)
## Nikolaj Schejbel Nielsen   (s214494)
##
## The purpose of this script is to simulate 
## an interferometer under far-field conditions
##


#%% Import of libraries used 
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.dpi'] = 300

#------------------------------------------------------------------------------
#%% Defining custom funcitons

def placeSource(_sky, _theta, _intensity):
    # Function to place sources at a specified anglular position (in degrees)
    # Input:
    # _sky        =  Array containing sky information
    # _theta      =  Angular postion of source [degrees]
    # _intensity  =  Intensity of placed source
    _theta = _theta * np.pi/180
    _i = np.sum(theta < _theta)
    _sky[_i] = _intensity
    return


def VisibilitySim(_lmbda, _bw, _D, _N, _sky, _theta):
    # Simulation of visibility sampling. Assumes antenna spacing goes as u = n * D
    # Input: 
    # _lmdba = wavelength [m]
    # _bw    = bandwidth [Hz]
    # _D     = Antenna spacing [m]
    # _N     = Number of samples
    # _sky   = Simulated sky
    # _theta = Angular coordinate 

    
    _Vu      = np.zeros(_N, dtype = complex) 
    _D_noise = np.zeros(_N)
    
    for _i in range(_N):
        # Fringe function with bandwidth errors
        _Fringe      = np.exp(-2j*np.pi*(_i + 1)*(_D) * np.sin(_theta)/_lmbda)   *   np.sinc(((_i + 1)*(_D)) * np.sin(_theta) * _bw/3e+08 )        
        # Visibility for givne antenna spacing
        _Vu[_i] = np.sum(_Fringe * _sky)
        
    return _Vu 

def IntensityReconstructionAngular(_lmbda, _phi, _D, _Vu):
    # Reconstruction of intensity from complex visibility function, _Vu
    # Input:
    # _lmbda            = Wavelength [m]
    # _phi              = Angular coordinate [rad]
    # _D                = Antenna spacing [m]
    # _Vu               = Complex visibilities
    # Return: 
    # _I                = Recontstructed intensities over _phi

    _n     = np.size(_phi)                                  # Number of points to reconstruct over
    _u     = (np.arange(0, np.size(_Vu)) + 1) * _D/_lmbda   # Array with antenna spacings
    _I     = np.zeros(_n, dtype=complex)                    # Emtpy intensity array
    _I_tmp = _I
    
    # Reconstructs intesities based on angle, using DFT
    for _i in range(_n):
        _fexp = np.exp(2j * np.pi * _phi[_i] * _u)
        _I_tmp[_i] = 1/(N) *( np.sum(_Vu * _fexp )) 
    _I = _I_tmp + _I

    return _I


#------------------------------------------------------------------------------
#%% Defining the sky
# Vinkel mellem -Pi/2 og Pi/2 der definerer vores himmel opdelt i 10000 punkter
theta = np.linspace(-np.pi/2, np.pi/2, num = 10000)
dtheta = np.abs(theta[0] - theta[1])


# Simulering af signal sources på himlen 
Sky = np.zeros(np.size(theta))

#placeSource(Sky, 0, 1)              # Coordinates used for single point source figure

placeSource(Sky, -7.5, 1)          # Coordinates used for multisource plot
placeSource(Sky, -2.5, 0.5)
placeSource(Sky, 5, 1)

# 55 points ~ one degree
# middle at 5000
#Sky[4400-2*55:4400+2*55] = 1        # Coordinates used for extended source figure
#Sky[4800-2*55:4800+2*55] = 1
#Sky[5400-2*55:5400+4*55] = 1


#------------------------------------------------------------------------------
#%% Array setup

# Parameters for setup 
D     = 0.06        # Minimum distance between antennas [m]
lmbda = 0.03        # Wavelength [m]
N     = 16          # Number of samples 



# Simulation of the measured correlations / Visibility function (V(u))
Vu     = VisibilitySim(lmbda, 0, D, N, Sky, theta)
Vu_27M = VisibilitySim(lmbda, 27e+06, D, N, Sky, theta)
Vu_02G = VisibilitySim(lmbda, 02e+09, D, N, Sky, theta)
 



#------------------------------------------------------------------------------
#%% Intensity reconstructions
n = theta

I     = IntensityReconstructionAngular(lmbda, theta, D, Vu)
I_27M = IntensityReconstructionAngular(lmbda, theta, D, Vu_27M)
I_02G = IntensityReconstructionAngular(lmbda, theta, D, Vu_02G)

fov   = np.array([-lmbda/(D), lmbda/(D)]) * 90/np.pi   # Field of view of the array [degrees]

#------------------------------------------------------------------------------
#%% Figures

if True:
    # Simulated visibility
    fig, ax = plt.subplots()
    
    plt.xlabel('Antenna distance [wavelengths]')  

    ax.plot((np.arange(N) + 1)*2, np.real(Vu)/(np.max(np.abs(Vu))), marker = 'x', c = '#069AF3', label = 'Real(Vu)')
    secax  = ax.twinx()
    secax.plot((np.arange(N) + 1)*2, np.angle(Vu), marker = 'x', c = 'gray', label = 'Phase', ls = '--')
       
    ax.grid(True, axis = 'y')
    
    ax.set_ylabel('Real part of visibility [Normalized]', c =  '#069AF3')
    secax.set_ylabel('Phase of visibility [rad]', c =  'gray')
    
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = secax.get_legend_handles_labels()
    plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)



if True:
    # Plot of the reconstructed image
    fig, ax = plt.subplots()
    
    plt.plot(n * 180/np.pi, (np.real(I)) / np.max(np.abs(np.real(I))), c = '#069AF3', label = 'Image, BW = 0 Hz')
    # plt.plot(n * 180/np.pi, (np.real(I_27M)) / np.max(np.abs(np.real(I))), c = '#E50000', label = 'Image, BW = 27 MHz', ls = '--')
    plt.plot(n * 180/np.pi, (np.real(I_02G)) / np.max(np.abs(np.real(I))), c = '#E50000', label = 'Image, BW = 2 GHz', ls = '-')
    
    plt.plot(n * 180/np.pi, Sky/np.max(Sky), c = 'grey', ls = '--', label = "Source distribution")
    
    plt.scatter(fov, np.array([0, 0]), color='grey', marker = 'x')
    # plt.title('Multiple sources image')
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Normalised intensity')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-15, 15)
    



#------------------------------------------------------------------------------
#%% Point source response

u = (np.arange(0, N) + 1) * D/lmbda
B = np.zeros(np.size(n))
for i in range(N):
    # Simplified reconstruction
    B_tmp = np.cos(2 * np.pi * u[i] * n)
    B = B + B_tmp
    
if True:
    # Point source response figure
    fig, ax = plt.subplots()
    plt.plot(n * 180/np.pi, B/np.max(B),c = '#069AF3', label = "Point source response")
    plt.xlim(-15, 15)
    plt.grid(True)
    # plt.title('Point source response')
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Normalised intensity')
    
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    

#------------------------------------------------------------------------------
#%% Bandwidth scaling at max baseline 
# Scaling function:
Rbw_27M = np.sinc((2*N-2)*D * np.sin(n/2) * 27e+06/3e+08 )
Rbw_10G = np.sinc(2*N*D * np.sin(n) * 02e+09/3e+08 )

if True:
    fig, ax = plt.subplots()
    
    plt.plot(n * 180/np.pi, Rbw_27M, c = '#069AF3', label = 'Scaling, BW = 27 MHz')
    plt.plot(n * 180/np.pi, Rbw_10G, c = '#E50000', label = 'Scaling, BW = 2 GHz', ls = '-')


    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Scaling of fringe function')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-15*2, 15*2)
    
   
#------------------------------------------------------------------------------
#%%  # Bandwidth scaling at multiple baselines 
if True: 
    import matplotlib.cm as cm
    cmap_blue = cm.get_cmap('winter')    
    cmap_gray = cm.get_cmap('Greys')
    
    # Indexes to colorize and label
    idx     = np.array([0, 6, 9, 11, 13, 15])
    # Antenna distancess
    D_array = np.linspace(0, D * (N-1), num = N) + 0.5 * D  
    
    fig, ax = plt.subplots()
    for i in range(N):
        # Scaling function:
        Rbw_27M = np.sinc((2*N-2)*D_array[i] * np.sin(n/2) * 27e+06/3e+08)
        
        if np.any(idx == i):
            plt.plot(n * 180/np.pi, Rbw_27M, c = cmap_blue(1.5-(i+4)/(N)),label = "D = {:.2f}".format(2*D_array[i]), zorder = 3)
        else:
            plt.plot(n * 180/np.pi, Rbw_27M, c = cmap_gray((i+1)/(N*2)))
            
    plt.grid(True)
    plt.ylabel('Bandwidth attenuation')
    plt.xlabel('Angular position [degrees]')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-8, 8)
    #plt.ylim([0.94, 1.01])
