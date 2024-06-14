## 
## dataRead.py
##
## Code created by: 
## William Hedegaard Langvad  (s214512)
## Nikolaj Schejbel Nielsen   (s214494)
##
## The purpose of this script is to analyse the radiometer data
## using functions from synthesisFunctions.py, and present the data
## in various figures. 
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

# Import custom function from synthesisFunctions.py
from synthesisFunctions import loadVu, IntensityReconstructionAngular, hogbomCLEAN, antennaCorrection, hm_to_ssm

# Set dpi for sharper figures
plt.rcParams['figure.dpi'] = 300

#------------------------------------------------------------------------------
#%% Read data files and calculate Vu

# File paths for the data, relative to file placement 
dir_path_T   = '../data/240508-Measurements/with_target/meas_*_out/'  
dir_path_B   = '../data/240508-Measurements/without_target/meas_*_out/' 
 
dir_Temp     = '../data/240508-Measurements/'


# Timestamps and measured temperatures of CalH. Mainly used for 'corrCoeffT' setting used in earlier drafts. 

# Time and temperatures for 240507 measurements
#time  = np.array(['14:29', '14:38', '14:55', '15:01', '15:12', '15:18', '15:23', '15:28', '15:33', '15:41', '15:46', '15:51', '15:55', '16:00', '16:11', '16:18', '16:24', '16:28', '16:33', '16:37', '16:41', '16:46', '16:51', '16:56', '17:00', '17:05', '17:10', '17:14', '17:19', '17:25'])
#time  = hm_to_ssm(time)
#temp2 = np.array([38.5, 38.6, 38.7, 38.5, 38.6, 38.7, 38.6, 38.7, 38.7, 38.7, 38.7, 38.6, 38.8, 38.6, 38.7, 38.6, 38.6, 38.6, 38.6, 38.7, 38.7, 38.5, 38.5, 38.7, 38.5, 38.6, 38.6, 38.5, 38.5, 38.7]) + 273.15
#temp8 = np.array([37.3, 37.4, 37.4, 37.4, 37.4, 37.6, 37.4, 37.5, 37.4, 37.4, 37.4, 37.4, 37.4, 37.5, 37.4, 37.4, 37.4, 37.5, 37.6, 37.4, 37.4, 37.4, 37.4, 37.3, 37.4, 37.3, 37.4, 37.3, 37.5, 37.4]) + 273.15


# Time and temperatures for 240508 measurements
time  = np.array(['13:03', '13:08', '13:12', '13:16', '13:19', '13:29', '13:33', '13:43', '13:49', '13:54', '13:58', '14:02', '14:09', '14:15', '14:18', '14:24', '14:40', '14:45', '14:48', '14:53', '15:03', '15:12', '15:15', '15:19', '15:23', '15:27', '15:32', '15:37', '15:41', '15:45', '15:52', '15:58'])  
time  = hm_to_ssm(time)
temp2 = np.array([ 37.1,  37.5,  37.6,  37.9,  38.1,  38.3,  38.4,  38.5,  38.4,  38.5,  38.4,  38.6,  38.6,  38.6,  38.5,  38.6,  38.5,  38.5,  38.6,  38.7,  38.7,  38.7,  38.7,  38.7,  38.6,  38.5,  38.6,  38.6,  38.6,  38.6,  38.5,  38.7]) + 273.15
temp8 = np.array([ 35.9,  36.5,  36.6,  36.8,  36.8,  37.4,  37.3,  37.3,  37.3,  37.3,  37.3,  37.4,  37.3,  37.4,  37.4,  37.3,  37.4,  37.4,  37.5,  37.4,  37.3,  37.3,  37.4,  37.3,  37.4,  37.5,  37.5,  37.5,  37.4,  37.4,  37.4,  37.4]) + 273.15

# Time and temperatures for 240515 measurements 
# The timing in the data is wrong, so a 28294 second delay needs to be added to match the timing in the radiometer output files.
#time  = np.array(['23:36:25', '23:40:30', '23:45:02', '23:49:57', '23:53:37', '23:58:07', '00:03:10', '00:06:52', '00:11:51', '00:16:31', '00:20:50', '00:25:03', '00:29:40', '00:35:00', '00:39:36', '00:54:23', '00:58:05', '01:02:05', '01:07:26', '01:11:56', '01:16:36', '01:20:43', '01:34:34', '01:38:41', '01:42:41', '01:46:47', '01:52:14', '01:56:46', '02:01:08', '02:05:21', '02:09:14', '02:14:25']) 
#time  = hm_to_ssm(time)   
#temp2 = np.zeros(np.size(time))
#temp8 = np.array([ 36.2, 36.3, 36.3, 36.4, 36.5, 36.5, 36.6, 36.6, 36.5, 36.6, 36.5, 36.6, 36.5, 36.7, 36.6, 36.5, 36.5, 36.5, 36.5, 36.5, 36.6, 36.5, 36.6, 36.7, 36.6, 36.6, 36.6, 36.7, 36.7, 36.7, 36.7, 36.9]) + 273.15

temp_calH = np.array([time, temp2, temp8])



# Visibilites, timestamps, mean CalH values and antenna temperatures are calculated for target and background measurements
Vu_T, measTime_T, mean_calH_T, tempUsed_T = loadVu(dir_path_T, dir_Temp, temp_calH, setting = 'corrCoeffTInterp', plotMeasuremetnsCombined = False)
Vu_B, measTime_B, mean_calH_B, tempUsed_B = loadVu(dir_path_B, dir_Temp, temp_calH, setting = 'corrCoeffTInterp')


# Correcting the visibilites for cable loss
l    = 10**(-0.82*1.5/10)   # Cable loss
Vu_T = Vu_T/l
Vu_B = Vu_B/l

# Calculating the difference in visibilites 
Vu = Vu_T - Vu_B

#------------------------------------------------------------------------------
#%% Setup of dimensions used in the experiment

d       = 9.0        # Distance to target [m]
D       = 0.06       # Minimum distance between two antennas [m]
lmbda   = 0.03       # Bølgelængde [m]
N       =   16       # Number of measurements
w  =   2             # Width of the target field, [m]
n  = 10000           # Resolution of reconstructions

# Calculating antenna coordinates and saving into antenna_coord array
antenna_coord = np.zeros([2, N])  
D_array = np.linspace(0, D * (N-1), num = N) + 0.5 * D  
for i in range(N):
    x = D_array[i]
    # y = np.sqrt(d**2 - x**2)
    y = d
    antenna_coord[:, i] = [x, y]
    
# Parameters of source field, used for layover in figures
# Source field given in x, y, source intensity
Sf        = np.zeros([3, n])  
Sf[0, :]  = np.linspace(-w, w, n) * 0.5
Sf[1, :]  = 0
Sf[2, :] = 10
if True:
    Sf[2,    0:375]  = 10 * (1 - 0.1) + (0.1 * (273 + 15))
    Sf[2,  375:1125] = 273 + 15
    Sf[2, 1125:1300] = 10 * (1 - 0.1) + (0.1 * (273 + 15))
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

#------------------------------------------------------------------------------
#%% Reconstructing intensityand making last corrections 

# Calculating Field of view, and defing coordinates to reconstruct over
fov = np.array([-lmbda/(D), lmbda/(D)]) / 4                     # fov, wavelength over smallest distance between two antennas
phi = np.linspace(-lmbda/(D), lmbda/(D), num = n * 3) * 3 / 4   # Angular coordinate


# Reconstructed intensity and CLEANed image
I = IntensityReconstructionAngular(lmbda, phi, antenna_coord, Vu)
clean_im, x_coord = hogbomCLEAN(lmbda, phi, fov, antenna_coord, I, Sf)


# Determining target offset through convolution
m  = np.linspace(-w, w, n * 3) * 0.5 * 3                         # Target coordinates
conv = np.correlate(np.real(I)/np.max(I), Sf[2,:]/np.max(Sf[2,:]), mode = 'same')
idx  = np.argmax(conv)
offset = m[idx]

# Antenna correction coeffictient
vD, hD = antennaCorrection(lmbda, phi)

G_0   = 10**(6/10)              # Gain coefficient, ~6dB
Ae    = 0.0002322576
hD    = hD * G_0 / np.max(hD)   # Scale directivity to gain

coeff = 1/(hD * vD)

#------------------------------------------------------------------------------
#%% Figures used in the report
if True:
    # Plot of the visibility function
    fig, ax = plt.subplots()
    
    ax.plot(2 * antenna_coord[0,:]/lmbda, np.real(Vu), marker = 'x', c = '#069AF3', label = 'Real(Vu)')
    secax  = ax.twinx()
    secax.plot(2 * antenna_coord[0,:]/lmbda, np.angle(Vu), marker = 'x', c = 'gray', label = 'Phase', ls = '--')
    
    ax.grid(True, axis = 'y')
    ax.set_xlabel('Antenna distance [wavelengths]')
    ax.set_ylabel('Real part of visibility', c =  '#069AF3')
    secax.set_ylabel('Phase of visibility [rad]', c =  'gray')
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = secax.get_legend_handles_labels()
    # plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc=0)
    plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)

    plt.show()
    plt.close()


if True:
    # Background and target visibility    
    # Target
    fig, ax = plt.subplots()
    
    ax.plot(2 * antenna_coord[0,:]/lmbda, np.real(Vu_T), marker = 'x', c = '#069AF3', label = 'Real(Vu)')
    secax  = ax.twinx()
    secax.plot(2 * antenna_coord[0,:]/lmbda, np.angle(Vu_T), marker = 'x', c = 'gray', label = 'Phase', ls = '--')
    
    ax.grid(True, axis = 'y')
    ax.set_xlabel('Antenna distance [wavelengths]')
    ax.set_ylabel('Real part of visibility', c =  '#069AF3')
    secax.set_ylabel('Phase of visibility [rad]', c =  'gray')
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = secax.get_legend_handles_labels()
    # plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc=0)
    plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)

    plt.show()
    plt.close()
    
    
    # Background
    fig, ax = plt.subplots()
    
    ax.plot(2 * antenna_coord[0,:]/lmbda, np.real(Vu_B), marker = 'x', c = 'red', label = 'Real(Vu)')
    secax  = ax.twinx()
    secax.plot(2 * antenna_coord[0,:]/lmbda, np.angle(Vu_B), marker = 'x', c = 'gray', label = 'Phase', ls = '--')
    
    ax.grid(True, axis = 'y')
    ax.set_xlabel('Antenna distance [wavelengths]')
    ax.set_ylabel('Real part of visibility', c =  'red')
    secax.set_ylabel('Phase of visibility [rad]', c =  'gray')
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = secax.get_legend_handles_labels()
    # plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc=0)
    plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)

    plt.show()
    plt.close()


if True:
    # Plot of the reconstructed image
    fig, ax = plt.subplots()
    
    plt.plot(phi * 180/np.pi, np.real(I) * coeff, c = '#069AF3', label = 'Image')
    plt.plot((np.arctan((Sf[0,:]+offset)/d)) * 180/np.pi , Sf[2,:] / np.max(Sf[2,:]) * np.max(np.real(I * coeff)), c = 'gray', ls = '--', label = 'Source distribution')
    plt.scatter(fov * 180/np.pi, np.array([0, 0]), color='black', marker = 'x', zorder = 3)
    
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Intensity')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-fov[1] * 1.1 * 180/np.pi, fov[1] * 1.1 * 180/np.pi)
    plt.ylim(bottom = -0.3 * np.max(np.real(I * coeff)))
  
    plt.show()
    plt.close()
    

if True:
    # Plot of the reconstructed image without source distribution
    fig, ax = plt.subplots()
    
    plt.plot(phi * 180/np.pi, np.real(I) * coeff, c = '#069AF3', label = 'Image')
    plt.scatter(fov * 180/np.pi, np.array([0, 0]), color='black', marker = 'x', zorder = 3)
    
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Intensity')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-fov[1] * 1.1 * 180/np.pi, fov[1] * 1.1 * 180/np.pi)
    plt.ylim(bottom = -0.3 * np.max(np.real(I * coeff)))
  
    plt.show()
    plt.close()
    
if True:
    # Background and target reconstruction 
    I_t = IntensityReconstructionAngular(lmbda, phi, antenna_coord, Vu_T)
    I_b = IntensityReconstructionAngular(lmbda, phi, antenna_coord, Vu_B)
    
    fig, ax = plt.subplots()

    plt.plot(phi * 180/np.pi, np.real(I_t) * coeff, c = '#069AF3', label = 'With target', zorder = 3)    
    plt.plot(phi * 180/np.pi, np.real(I_b) * coeff, c = 'red', label = 'Without target')
    # plt.plot((np.arctan((Sf[0,:]+offset)/d)) * 180/np.pi , Sf[2,:] / np.max(Sf[2,:]) * np.max(np.real(I * coeff)), c = 'gray', ls = '--', label = 'Source distribution')
    plt.scatter(fov * 180/np.pi, np.array([0, 0]), color='black', marker = 'x', zorder = 3)
    
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Intensity')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-fov[1] * 1.1 * 180/np.pi, fov[1] * 1.1 * 180/np.pi)
    plt.ylim(bottom = -0.3 * np.max(np.real(I_t * coeff)))
  
    plt.show()
    plt.close()
    

if True:
    # Plot of the CLEANed image
    fig, ax = plt.subplots()
    
    plt.plot(x_coord * 180/np.pi, (np.real(clean_im)) / (1+0 * np.max(np.abs(np.real(clean_im)))), c = '#069AF3', label = 'CLEAN Image')
    plt.plot((np.arctan((Sf[0,:]+offset)/d)) * 180/np.pi , Sf[2,:] / np.max(Sf[2,:]), c = 'gray', ls = '--', label = 'Source distribution')
    plt.scatter(fov * 180/np.pi, np.array([0, 0]), color='black', marker = 'x', zorder = 3)
    
    plt.xlabel('Angular position [degrees]')
    plt.ylabel('Normalised intensity')
    plt.grid(True)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    plt.xlim(-fov[1] * 1.1 * 180/np.pi, fov[1] * 1.1 * 180/np.pi)
    plt.ylim(bottom = -0.3)
  
    plt.show()
    plt.close()


if True:
    # Plot of CalH counts / Actual temperature of CalH
    fig, ax = plt.subplots()

    ax.scatter(measTime_T - measTime_T[0], mean_calH_T, marker = 'x', c = '#069AF3', label = 'W. Target')
    ax.scatter(measTime_B - measTime_T[0], mean_calH_B, marker = 'x', c = '#E50000', label = 'Background')
    
    secax  = ax.twinx()
    secax.scatter(temp_calH[0,:] - temp_calH[0,0], temp_calH[2,:], marker = '+', c = 'grey', label = 'Measured temperature')

    ax.set_ylabel('Counts')
    secax.set_ylabel('Measuered temperature [Kelvin]', c =  'gray')
    ax.set_xlabel('Seconds since start')    
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = secax.get_legend_handles_labels()
    plt.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    ax.grid(True, axis = 'y')
    
    plt.show()
    plt.close



if True: 
    # Antenna temperatures used for visivility correction
    fig, ax = plt.subplots()
    
    ax.scatter(np.arange(N) + N + 1, tempUsed_B, marker = 'x', c = '#E50000', label = 'Background')
    ax.scatter(np.arange(N) + 1, tempUsed_T, marker = 'x', c = '#069AF3', label = 'W. Target')
    
        
    plt.grid('on')
    ax.set_ylabel('Antenna temperature')
    ax.set_xlabel('Measurement #')  
    
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)

#------------------------------------------------------------------------------
#%% Additional figures, showing the reconstructed intensity as an ´actual´ image


if True:
    # Flat reconstruction of the intensities, 
    y_data = np.tile(np.real(I[np.abs(phi) < fov[1]]),(2,1))
    x_data = phi[np.abs(phi) < fov[1]]
    x_data = x_data / np.max(x_data)
    
    plt.figure()
    plt.title("pixel_plot") 
    pixel_plot = plt.imshow(y_data, cmap='gray_r', interpolation='nearest', origin='lower',extent = [np.min(x_data),np.max(x_data),0,1],aspect=1, vmin=0) 
    plt.colorbar(pixel_plot) 
    pixel_plot.axes.get_yaxis().set_visible(False)

if True:
    # Flat reconstruction of the CLEANed intensities,
    # Image used for front page
    y_data = np.tile(np.real(clean_im[np.abs(x_coord) < fov[1]]),(2,1))
    x_data = x_coord[np.abs(x_coord) < fov[1]] * 1/np.max(x_coord)
    
    plt.figure()
    # plt.title("pixel_plot") 
    pixel_plot = plt.imshow(y_data, cmap='inferno', interpolation='nearest', origin='lower',extent = [np.min(x_data),np.max(x_data),0,1],aspect=1.5, vmin=0.2, vmax = 1) 
    #plt.colorbar(pixel_plot) 
    #pixel_plot.axes.get_yaxis().set_visible(False)
    plt.axis('off')
    

    
    