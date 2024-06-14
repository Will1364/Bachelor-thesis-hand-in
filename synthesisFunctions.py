## 
## synthesisFunctions.py
##
## Code created by: 
## William Hedegaard Langvad  (s214512)
## Nikolaj Schejbel Nielsen   (s214494)
##
## This file contains all the important functions used across scripts in
## simulations and data analysis. These have been collected in one file to 
## make the remaining scripts easier to work with
##

#%% Import of important libraries used
import matplotlib.pyplot as plt
import numpy as np
import glob
import pandas as pd
import os
from scipy.interpolate import CubicSpline, interp1d


#------------------------------------------------------------------------------
#%% Functions used in the script
def visibilitySim(_lmbda, _antenna_coord, _Sf, _Nf, _bw = 0):
    # Simulates the visibility function (V(u)) based on the parameterized antenna setup and source field. 
    # Assumes that the antenna coordinates should be mirrored around y-axis.  
    # Input: 
    # _lmbda            = Wavelength [m]
    # _bw               = Instrument bandwidth [Hz]
    # _antenna_coord    = Coordinates of antennas, assumes that the antenna coords are mirrored around middle of the array
    # _Sf               = Array with the sourcefield (x, y, intensity)
    # _Nf               = Array with the noise field (x, y, intensity), must be same size as _Sf
    # Return: 
    # _Vu               = Array with complex visibilities, same length as antenna_coord
    
    _N  = np.size(_antenna_coord, 1)        # Number of measurements
    _n  = np.size(_Sf, 1)                   # Size of source (and noise) field
    _Vu = np.zeros(_N, dtype = complex)     # Creates array to store visibility function
    
    
    # Determines the visibility for each of the antenna pairs
    for _i in range(_N):
        
        _fringe   = np.zeros(_n, dtype = complex)       # Creates array to store fringe function for the source field
        _fringeNf = np.zeros(_n, dtype = complex)       # Creates array to store fringe function for the noise field

        # Calculates the distance from the two antenna to source, finds the difference and determines the complex correlated signal (Fringe function)
        for _j in range(_n):
            _d1     =  np.sqrt((_Sf[0, _j] + _antenna_coord[0, _i])**2 + (_Sf[1, _j] - _antenna_coord[1, _i])**2)
            _d2     =  np.sqrt((_Sf[0, _j] - _antenna_coord[0, _i])**2 + (_Sf[1, _j] - _antenna_coord[1, _i])**2)
            _d_diff =  (_d2 - _d1)
        
            _fringe[_j] = np.exp(2j * np.pi * _d_diff / _lmbda) * np.sinc( _d_diff * _bw/3e+08 )
            
        # Calculates the distance from the two antenna to the noise field, finds the difference and determines the complex correlated signal (Fringe function)  
        for _j in range(_n):
            _d1     =  np.sqrt((_Nf[0, _j] + _antenna_coord[0, _i])**2 + (_Nf[1, _j] - _antenna_coord[1, _i])**2)
            _d2     =  np.sqrt((_Nf[0, _j] - _antenna_coord[0, _i])**2 + (_Nf[1, _j] - _antenna_coord[1, _i])**2)
            _d_diff =  (_d2 - _d1)
        
            _fringeNf[_j] = np.exp(2j * np.pi * _d_diff / _lmbda) 
            
        # Calculates Vu based on field intensities and the fringe functions. 
        _Vu[_i] = np.sum(_Sf[2, :] * _fringe + _Nf[2, :] * _fringeNf)
        
    # Returns the simulated fringe function
    return _Vu

def VisibilitySimAngular(_lmbda, _bw, _antenna_coord, _sky, _theta, _noise = 0):
    # Simulates the visibility function for sources in the far field (sky)
    # Input: 
    # _lmbda            = Wavelength [m]
    # _bw               = Instrument bandwidth [Hz]
    # _antenna_coord    = Coordinates of antennas, assumes that the antenna coords are mirrored around middle of the array
    # _theta            = Angular coordinates used [rad]
    # _sky              = Array with information about sources intensity at _theta
    # _noise            = Adds noise in antenna placement, with errors of +-_noise
    # Return: 
    # _Vu               = Array with complex visibilities, same length as antenna_coord
    
    # Simulering af det målte signal / Visibility funktionen (V(u))
    _N       = np.size(_antenna_coord[0,:])
    _Vu      = np.zeros(_N, dtype = complex) 
    _D_noise = np.zeros(_N)
    
    for _i in range(_N):
        # Antenna placement errors
        _n           = (np.random.rand(1) * 2 - 1) * _noise         # Amplitude of the noise, distributed evenly around zero    
        _D_noise[_i] = _n[0]                                        # Saves the noise into an array so that it can be looked at again
        
        # Determines the fringe function for each spacing, and applies bandwidth effect.
        _Fringe = np.exp(-2j*np.pi*(2 * _antenna_coord[0, _i] ) * np.sin(_theta)/_lmbda)   *   np.sinc((2 * _antenna_coord[0, _i] ) * np.sin(_theta) * _bw/3e+08 )        
        
        # Calculates Visibilitues
        _Vu[_i] = np.sum(_Fringe * _sky)
        
    return _Vu


def IntensityReconstruction(_lmbda, _m, _antenna_coord, _Vu, plotIndividual = False):
    # Reconstruction of intensity from complex visibility function, _Vu, for a parameterized setup
    # Corrects for the effects of the near-field. 
    # Input:
    # _lmbda            = wavelength (meter)
    # _m                = Coord to reconstruct over
    # _antenna_coord    = Coordinates of antennas, assumes that the antenna coords are mirrored around middle of the array
    # _Vu               = Complex visibility function
    # plotIndividual    = Creates a plot with the individual signals for Vu
    # Output: 
    # _I                = Recontructed intensities
    
    
    _N  = np.size(_antenna_coord, 1)        # Number of measurements 
    # Coordinate correction for the antenna array, used to 'focus' the measurements on the source field
    _xi = 2 * _antenna_coord[0, :] / (np.sqrt(_antenna_coord[0, :]**2 + (_antenna_coord[1, :])**2)) 
    _u  = 1/_lmbda 
    
    # Empty array to save intensity into
    _I  = np.zeros(np.size(_m), dtype=complex) 
    _temp_I = np.zeros(np.size(_m), dtype=complex) 
    
    if plotIndividual:
        fig, ax = plt.subplots()


    # Reconstructs intensities over the target (in coordinates on target)
    for _j in range(np.size(_xi)):
        for _i in range(np.size(_m)):
            _fexp = np.exp(2j * np.pi * (_m[_i] * _xi[_j] * _u))
            # 
            _temp_I[_i] = 1/_N * (np.sum(_Vu[_j] * _fexp ))
            
        if plotIndividual:
            plt.plot(_m, np.real(_temp_I))  
        _I = _I + _temp_I
        
    if plotIndividual: 
        plt.show()
        plt.close()
        
    return _I

def IntensityReconstructionAngular(_lmbda, _phi, _antenna_coord, _Vu):
    # Reconstruction of intensity from complex visibility function, _Vu
    # Input:
    # _lmbda            = Wavelength [m]
    # _phi              = Angular coordinate
    # _antenna_coord    = Coordinates of antennas, assumes that the antenna coords are mirrored around middle of the array. 
    # _Vu               = Complex visibilites, array. 
    # Output: 
    # _I                = Recontructed intensities
    
    _n  = np.size(_phi)                     # Number of coordinate points
    _u  = 2 * _antenna_coord[0,:] / _lmbda  # Antenna spacing in wavelengths
    _I  = np.zeros(_n, dtype=complex)       # Empty array to save intensity into

    # Reconstructs intesities based on angle
    for _i in range(_n):
        _fexp = np.exp(2j * np.pi * _phi[_i] * _u)
        _I[_i] = 1/(np.size(_u)) *( np.sum(_Vu * _fexp ))
 
    return _I

def dirtyBeam(_lmbda, _m, _antenna_coord, _Sf_input):
    # Simulates a dirty beam / point source response for the antenna setup. The point source is centered in the middle of the source field. 
    # Used in hogbomCLEAN function. 
    # Input:
    # _lmbda            = wavelength (meter)
    # _m                = Coord to reconstruct over
    # _antenna_coord    = Coordinates of antennas, assumes that the antenna coords are mirrored around middle of the array
    # _Sf_input         = Array with the sourcefield (x, y, intensity) - used to reconstruct same size source field but with single point source
    # Output: 
    # _dBeam            = Dirty beam / point source response, normalized. 
    
    # Moves _Sf_input into a temporary variable so it doesn't get overwritten
    _Sf = np.array(_Sf_input)
    
    _N  = np.size(_antenna_coord, 1)        # Number of measurements
    _n  = np.size(_Sf, 1)                   # Source field array size
    _Vu = np.zeros(_N, dtype = complex)     # Visibility function
    _Sf[2, :] = 0                           # Ressetting source field and inserting point source 
    _Sf[2, int(np.size(_Sf, 1)/2)] = 1      # 
    
    # Determining Fringe and Visibility function
    for _i in range(_N):
        _fringe   = np.zeros(_n, dtype = complex)
        
        # Calculates the distance from the two antenna to source, finds the difference and determines the complex correlated signal (Fringe function)
        for _j in range(_n):
            _d1     =  np.sqrt((_Sf[0, _j] + _antenna_coord[0, _i])**2 + (_Sf[1, _j] - _antenna_coord[1, _i])**2)
            _d2     =  np.sqrt((_Sf[0, _j] - _antenna_coord[0, _i])**2 + (_Sf[1, _j] - _antenna_coord[1, _i])**2)
            _d_diff =  (_d2 - _d1)
        
            _fringe[_j] = np.exp(2j * np.pi * _d_diff / _lmbda) 

        # Calculated visibility
        _Vu[_i] = np.sum(_Sf[2, :] * _fringe)
    
    # Legacy near field reconstruction
    # Coordinate correction for the antenna array, used to 'focus' the measurements on the source field
    #_xi = 2 * _antenna_coord[0, :] / (np.sqrt(_antenna_coord[0, :]**2 + (_antenna_coord[1, :])**2)) 
    #_u  = 1/_lmbda 
    #_dBeam  = np.zeros(np.size(_m), dtype=complex)  # Dirty beam 
    # Reconstructing the dirty beam 
    #for _j in range(np.size(_xi)):
    #    for _i in range(np.size(_m)):
    #        _fexp = np.exp(2j * np.pi * (_m[_i] * _xi[_j] * _u))
    #        # Rekonstrueret intensitet på nattehimlen
    #        _dBeam[_i] = _dBeam[_i] + 1/_N * (np.sum(_Vu[_j] * _fexp ))
    #
    # Normalization of amplitude
    #_dBeam = np.real(_dBeam) * 1/np.max(np.real(_dBeam))
    
    # Reconstruction for angular coordinates
    _n  = np.size(_m)
    _u  = 2 * _antenna_coord[0,:] / _lmbda
    _dBeam  = np.zeros(_n, dtype=complex) 
    
    # Reconstructs intesities based on angle
    for _i in range(_n):
        _fexp = np.exp(2j * np.pi * _m[_i] * _u)
        _dBeam[_i] = 1/(np.size(_u)) *( np.sum(_Vu * _fexp ))
    
    _dBeam = np.real(_dBeam) * 1/np.max(np.real(_dBeam))        # Scales dirty beam to unity
    # plt.plot(_dBeam)
    
    return _dBeam

def hogbomCLEAN(_lmbda, _m, _fov, _antenna_coord, _I, _Sf):
    # Function to run the Högbom CLEAN algorithm. The function is limited close to the FOV as it doesn't handle negative intensities. 
    # Input:
    # _lmbda            = wavelength (meter)
    # _m                = Array with coord to reconstruct over (same as used in reconstruction of _I)
    # _antenna_coord    = Coordinates of antennas, assumes that the antenna coords are mirrored around middle of the array
    # _I                = Reconstructed intensities (dirty image)
    # _Sf               = Array with the sourcefield (x, y, intensity) - used to reconstruct same size source field but with single point source
    # Output: 
    # _clean_im         = CLEANed image 
    # _x_coord          = Corresponding x_coordinates for the image
    
    print('CLEANing image')
    
    _dirty_I = np.real(_I[np.abs(_m) < _fov[1] * 1.2])  # Dirty image, limited in size based on the fov of the array
    _max_I   = np.max(_dirty_I)                         # Maximum value of the image (a variable value is used now for better accuracy)
    _x_coord = _m[np.abs(_m) < _fov[1] * 1.2]           # x_coordinates for the limited dirty image
    
    _dirty_Beam = dirtyBeam(_lmbda, _m, _antenna_coord, _Sf)    # dirty beam determined with dirtyBeam function
    
    _im_size    = np.size(_dirty_I)     # size of dirty image array, used for positioning dirty beam and clean beam
    _beam_size  = np.size(_dirty_Beam)  # Dirty beam array size, used for positioning dirty beam and clean beam

    _clean_im     = np.zeros(_im_size)  # Array for storing clean image
    _clean_im_dlt = np.zeros(_im_size)  # Array for storing delta functions for clean image. (Should be convolded with the clean beam, although it is not done here)
    
    _FWHM = 2*np.abs(_m[int((np.size(_dirty_Beam) - np.sum(_dirty_Beam > 0.5))/2-1)])   # Full width half maximum of the dirty beam, used to determine size of CLEAN beam
    _a    = 4*np.log(2)/(_FWHM**2)                                                      # Temporary variable for calculating clean beam
    _clean_bm     = np.exp(-(_m**2)*_a)                                                 # CLEAN beam determined from FWHM
    
    _rms   = np.sqrt(np.mean(np.square(_dirty_I)))  # RMS of the dirty image 
    _gamma = 0.10                                   # Loop gain, used to scale the dirty beam and clean beam in algorithm
    
    _jj = 0     # Itteration counter 
    
    
    # The Högbom clean algorithm 
    # Finds maximum of dirty image, subtracts scaled dirty beam and inserts clean beam in the cleaned image at the same place. 
    # Repeats until the rms is no longer lowered or OR if the maximum value of the dirty image is less than the rms (can happen if there are negative values in the image)
    while True:
        _i = np.argmax(_dirty_I)    # Location of maxima in dirty image
        
        _cut = np.array([_beam_size/2 - _i, _beam_size/2 + (_im_size - _i)], dtype = int)   # Coordinates used to move dirty and clean beam to the correct location
        _scale = (np.max(_dirty_I) * _gamma)                                                # Scaling factor used, based on maximum of dirty image
        
        _dirty_I = _dirty_I - _scale  * _dirty_Beam[_cut[0]:_cut[1]]                        # Subtracts scaled dirty beam from dirty image       
        
        _rms_new = np.sqrt(np.mean(np.square(_dirty_I)))                                    # Calculates new RMS
        
        # Breaking conditions. The algorithm stops if it no longer lowers the rms of the image OR if the maximum value is less than the rms (happens if there are negative values in the image)
        if _rms_new >= _rms:
            _dirty_I = _dirty_I + _scale  * _dirty_Beam[_cut[0]:_cut[1]]        # Adds back scaled dirty beam to dirty image       
            print('\nBreak condition:  new rms of dirty image >= old rms')
            
            break
        if np.max(_dirty_I) < _rms:
            _dirty_I = _dirty_I + _scale  * _dirty_Beam[_cut[0]:_cut[1]]        # Adds back scaled dirty beam to dirty image       
            print('\nBreak condition:  dirty image amplitude < rms')            # Print break condition
            break
        
        _rms = _rms_new                                                         # Saves new rms
        
        _clean_im_dlt[_i] = _clean_im_dlt[_i] + _scale                          # Adds scaled delta function to clean image dlt (Can be used for convolution, although not done here)
        _clean_im        = _clean_im + _scale * _clean_bm[_cut[0]:_cut[1]]      # Adds clean beam to clean image
        
        _jj += 1 # Counter foes up by one
        
        print('\r# of Itterations: %i' %_jj, end = '')                          # Print break condition

    _clean_im = 1/np.max(_clean_im) * _clean_im + 1/_max_I * _dirty_I           # Scales CLEAN image
    
    
    # Optional figure of dirty and clean beam
    if False:
        fig, ax = plt.subplots()
        plt.plot(_m, _dirty_Beam, label = 'Dirty beam', c = '#069AF3')
        plt.plot(_m, _clean_bm, label = 'Clean beam', c = 'red')
        # plt.scatter(_Sf[0,:], np.zeros(np.size(_Sf[0,:])), c = _Sf[2,:], s=0.1)

        plt.title('Dirty og clean beam')
        plt.xlabel('Koordinat på target [m]')
        plt.ylabel('Normaliseret intensitet')
        plt.xlim(_fov[0]*0.5, _fov[1]*0.5)
        plt.legend(loc='upper left')
        plt.grid(True)
    
        plt.savefig('figs/temp_a.svg')
        
        
    print('done\n')
    return _clean_im, _x_coord

def linLSQ(_d, _t): 
    # LSQ linear function,  d = t * m[0] + m[1]
    # Input:
    # _d = Measurement/intensity
    # _t = Variable/time
    # Output:
    # _m = Linear function parameters 
        
    _G   = np.column_stack((_t, np.ones(np.size(_t)))) 
    _m   = np.matmul( np.linalg.inv(np.matmul(_G.T, _G)) , np.matmul(_G.T, _d.T))
    
    return _m

def hm_to_ssm(_time):
    # Hour minute in string format ('hh:mm' or 'hh:mm:ss') to seconds since midnight
    # For compatibility with radiometer file time formating
    # _time = array of time as strings 
    # _ssm  = array of time as seconds since midnight
    
    _n   = np.size(_time)
    _ssm = np.zeros(_n)
    
    # Checks for 'hh:mm' format
    if len(_time[0]) == 5:
        for i in range(_n):
            _ssm[i] = int(_time[i][0:2]) * (60**2) + int(_time[i][3:5]) * 60
    # Checks for 'hh:mm:ss' format
    elif len(_time[0]) == 8:
        for i in range(_n):
            _ssm[i] = int(_time[i][0:2]) * (60**2) + int(_time[i][3:5]) * 60 + int(_time[i][6:8])
    # Error message if format is wrong
    else:
        print("Error in time formatting, please inserting as 'hh:mm' or 'hh:mm:ss'")
        
    # Adds 24 hours if the minimum is not the first value (to account for going past midnight)
    _idx = np.argmin(_ssm)
    if _idx != 0:
        _ssm[_idx:] = _ssm[_idx:] + 86400
        
    return _ssm

def interpTemp(path_T2, path_T8, interpType):
    # This function interpolates temperature measurements for a given set of times
    # Input: 
    # File path of T2 and  T8 files, containing 
    #   Tidspunktet hvor forsøget startede [sekunder siden midnat]
    #   Tidsserie for antennemålinger [sekunder siden midnat]
    # Interpolation type ["linear" or "spline"]
    #   Linear is a straight line from point to point
    #   Cubic draws curved lines 
    # Output: Temperatures interpolated to the measurement times
    
    data_T2 = pd.read_csv(path_T2, delimiter=';')
    
    
    T2 = data_T2.T2.values  
    T2 = np.append(T2[0], np.append(T2, T2[np.size(T2)-1])) + 273.15
    time2 = data_T2.Time.values
    time2 = hm_to_ssm(time2)
    time2 = np.append(0, np.append(time2, time2[np.size(time2)-1]*2))


    data_T8 = pd.read_csv(path_T8, delimiter=';')
    
    T8 = data_T8.T8.values    
    T8 = np.append(T8[0], np.append(T8, T8[np.size(T8)-1])) + 273.15
    time = data_T8.Time.values
    time = hm_to_ssm(time)
    time = np.append(0, np.append(time, time[np.size(time)-1]*2))
          
    if interpType == "linear":
        interp_T2 = interp1d(time2, T2)
        interp_T8 = interp1d(time, T8)

    if interpType == "spline":
        interp_T8 = CubicSpline(time, T8)
        interp_T2 = CubicSpline(time2, T2)
        
    
    return interp_T2, interp_T8

def loadVu(dir_path, temp_path = '', tempMeas = None, setting = 'corrCoeff', plotMeasurements = False, plotMeasuremetnsCombined = False, specificIdx = None, prevCalH = None, prevCalV = None):
    # Loads data from radiometer measurement and calculates visibility function
    # dir_path is path containing data, temp_path path of temperature data
    # Settings: 
    #   'corrCoeff'         calculates Vu from the correlation coefficient 
    #   'corrCoeffT'        adds antenna temperature to the visibities
    #   'corrCoeffTInterp'  uses interpolated temperatures for CalH
    # plotMeasurements: 
    #   if True, plots the individual measurements used for computing Vu
    # plotMeasuremetnsCombined
    #   Plots all the measurement data in a combined figure
    # specificIdx:
    #   Insert specific order/measurements to load. Defaults to follow sorted order loaded by glob.glob
    # prevCalH/prevCalV:
    #   Makes it possible to input values for calibration using LSQ if a measurement is missing CalH/Nd data 
        
    print("Processing Vu for filepath '{}'".format(dir_path))
    
    # file extensions used by 'Unpack2.exe' 
    _ext_antA   = '.r01'  
    _ext_calH   = '.c01'
    _ext_antANd = '.n01'
    _ext_calHNd = '.n03'
    
    # Finding all measurements using glob, and then sorting correctly 
    _filenames_antA = glob.glob(dir_path + '*' + _ext_antA)
    _filenames_antA.sort()
    _filenames_calH = glob.glob(dir_path + '*' + _ext_calH)
    _filenames_calH.sort()
    _filenames_antANd = glob.glob(dir_path + '*' + _ext_antANd)
    _filenames_antANd.sort()
    _filenames_calHNd = glob.glob(dir_path + '*' + _ext_calHNd)
    _filenames_calHNd.sort()
    
    # Finds measurements where there is data from antA
    if np.size(specificIdx) == 1:
        _index = np.array([], dtype=int)
        for i in range(np.size(_filenames_antA)):
            if os.stat(_filenames_antA[i]).st_size != 0:
                _index = np.append(_index, [i])
    else: 
        _index = specificIdx
        
    # Print message depending on provided settings
    print('Setting "{}" chosen for Vu'.format(setting))
    
    if np.size(tempMeas) == 1 and (setting == 'corrCoeffT'):
        print('No temperature measurements provided. Using 300 K for all measurements')
        tempMeas = np.array([[0], [300], [300]])
    
    # Empty function to store visibility function and other output variables
    _Vu        = np.zeros(np.size(_index), dtype = complex)
    _mean_calH = np.zeros(np.size(_index), dtype = float)
    _measTime  = np.zeros(np.size(_index), dtype = float)
    _tempUsed  = np.zeros(np.size(_index), dtype = float)
    
    # Arrays for plotMeasuremetnsCombined 
    if plotMeasuremetnsCombined:
        data_antAt   = np.zeros([1, 2])
        data_calHt   = np.zeros([1, 2])
        data_calHNdt = np.zeros([1, 2])
    
    # Linear interpolation of temperature
    _file_Temp2 = temp_path + 'T2.csv'
    _file_Temp8 = temp_path + 'T8.csv'
    _lin_T2, _lin_T8 = interpTemp(_file_Temp2, _file_Temp8, 'linear')
    
    # Determining the visibility function based on the provided settings 
    for i in _index:
        # Progress message
        print('\rProcessing meas {:>3}/{:>3}'.format(i+1, np.size(_filenames_antA)), end='')
        
        # Clear variables for each loop
        _data_antA   = None
        _data_calH   = None
        _data_antANd = None
        _data_calHNd = None
        
        # Loads data into arrays
        _data_antA   = np.loadtxt(_filenames_antA[i])
        if os.stat(_filenames_calH[i]).st_size != 0:
            _data_calH   = np.loadtxt(_filenames_calH[i])
        #if os.stat(_filenames_antANd[i]).st_size != 0:
        #    _data_antANd = np.loadtxt(_filenames_antANd[i])
        if os.stat(_filenames_calHNd[i]).st_size != 0:
            _data_calHNd = np.loadtxt(_filenames_calHNd[i])
        
        # Save data into one long array if 'plotMeasuremetnsCombined' is chosen
        if plotMeasuremetnsCombined:
            if i == _index[0]:
                data_antAt =  np.column_stack((_data_antA[:,0],_data_antA[:,7]))
                data_calHt =  np.column_stack((_data_calH[:,0],_data_calH[:,7]))
                data_calHNdt = np.column_stack((_data_calHNd[:,0],_data_calHNd[:,7]))
            else: 
                data_antAt    = np.concatenate((data_antAt, np.column_stack((_data_antA[:,0],_data_antA[:,7]))))
                data_calHt    = np.concatenate((data_calHt, np.column_stack((_data_calH[:,0],_data_calH[:,7]))))
                data_calHNdt  = np.concatenate((data_calHNdt, np.column_stack((_data_calHNd[:,0],_data_calHNd[:,7]))))
          
        # Time of measurement  
        _measTime[i] = _data_antA[0,0]  
        
        # Calculates Vu depending on chosen setting
        if setting == 'corrCoeff':
            # Scaled based on correlation coefficient
            _Vu[i]   = np.mean((_data_antA[:,1] - 1j * _data_antA[:,2])/np.sqrt(_data_antA[:,7] * _data_antA[:,8]))
       
        
       
        elif setting == 'corrCoeffT':
            # Scaled based on corrCoeff and antenna temperatures 
            # LSQ for finding T_ant
            # First, correct calH temperatures are determined, based on closest time provided in 'tempMeas': 
            _idx = np.argmin(np.abs(tempMeas[0,:] - _data_antA[0,0]))
            _tempCalH    = tempMeas[2, _idx]
            _tempCalHNd  = _tempCalH + 150
            _tempUsed[i] = _tempCalH

            # Tests if there is data and calculates correction parameters 
            if (np.size(_data_calH) > 1) and (np.size(_data_calHNd) > 1):
                d = np.concatenate((_data_calH[:,7], _data_calHNd[:,7]))
                t = np.concatenate((np.ones(np.size(_data_calH,0)) * _tempCalH, np.ones(np.size(_data_calHNd,0)) * _tempCalHNd))
                _m_temp_H = linLSQ(d, t)
                
                d = np.concatenate((_data_calH[:,8], _data_calHNd[:,8]))
                t = np.concatenate((np.ones(np.size(_data_calH,0)) * _tempCalH, np.ones(np.size(_data_calHNd,0)) * _tempCalHNd))
                _m_temp_V = linLSQ(d, t)
            # Uses provided LSQ values if no data is found: 
            elif prevCalH != None and prevCalV != None:
                _m_temp_V = prevCalH
                _m_temp_H = prevCalV
            # Prints error message if data is missing, stops function
            else:
                print('Error: Missing calibration data, use corrCoef setting or insert values for prevCalH and prevCalV ')
                break
            
            # Corrects Vu with correlation coefficient 
            _Vu[i]   = _Vu[i]   = np.mean((_data_antA[:,1] - 1j * _data_antA[:,2]) / np.sqrt(_data_antA[:,7]*_data_antA[:,8]))
            
            # Calculates antenna temperature, saves it, and applies to Vu 
            _T_antAH = (np.mean(_data_antA[:,7]) - _m_temp_H[1])/_m_temp_H[0]
            _T_antAV = (np.mean(_data_antA[:,8]) - _m_temp_V[1])/_m_temp_V[0]
            _T_tmp   = np.sqrt(_T_antAH*_T_antAV)
            _tempUsed[i] = _T_tmp
            
            _Vu[i]   = _Vu[i] * np.sqrt(_T_antAH*_T_antAV)

            _tempUsed[i] = _T_tmp

        elif setting == 'corrCoeffTInterp':
            # CalH temperatures are found based on linear interpolation
            _tempCalH_H   = _lin_T8(_data_calH[0,0] + 28294)
            _tempCalHNd_H = _tempCalH_H + 150
            
            _tempCalH_V   = _lin_T2(_data_calH[0,0] + 28294)
            _tempCalHNd_V = _tempCalH_V + 150
            
            # Finds LSQ parameters for T_ant
            if (np.size(_data_calH) > 1) and (np.size(_data_calHNd) > 1):
                d = np.concatenate((_data_calH[:,7], _data_calHNd[:,7]))
                t = np.concatenate((np.ones(np.size(_data_calH,0)) * _tempCalH_H, np.ones(np.size(_data_calHNd,0)) * _tempCalHNd_H))
                _m_temp_H = linLSQ(d, t)
                
                d = np.concatenate((_data_calH[:,8], _data_calHNd[:,8]))
                t = np.concatenate((np.ones(np.size(_data_calH,0)) * _tempCalH_V, np.ones(np.size(_data_calHNd,0)) * _tempCalHNd_V))
                _m_temp_V = linLSQ(d, t)
            else:
                # Error message if data is missing
                print('Error: Missing calibration data, use corrCoef setting or corrCoefT setting and insert values for prevCalH and prevCalV ')
                break
            
            
            # Calculates antenna temperature, saves it, and applies to Vu 
            _Vu[i]   = _Vu[i]   = np.mean((_data_antA[:,1] - 1j * _data_antA[:,2])/np.sqrt(np.mean(_data_antA[:,7])*np.mean(_data_antA[:,8])))
            _T_antAH = (np.mean(_data_antA[:,7]) - _m_temp_H[1])/_m_temp_H[0]
            _T_antAV = (np.mean(_data_antA[:,8]) - _m_temp_V[1])/_m_temp_V[0]
            _T_tmp   = np.sqrt(_T_antAH*_T_antAV)
            _tempUsed[i] = _T_tmp
            
            _Vu[i]   = _Vu[i] * _T_tmp
            
        else: 
            # Warning message if setting provided is wrong, defaults to 'corrCoef'
            if i == _index[0]:
                print('\rError in applied setting for Vu, defaulting to "corrCoef" setting')
            _Vu[i]   = np.mean((_data_antA[:,1] - 1j * _data_antA[:,2])/np.sqrt(np.mean(_data_antA[:,7])*np.mean(_data_antA[:,8])))
        
        
        # Calculates and saves mean measured calH value
        _mean_calH[i] = np.nan
        if np.size(_data_calH) > 1:
            _mean_calH[i] = np.mean(_data_calH[:,7])
        
        
        if plotMeasurements:
            # Plots each individual measurement for easy check of the data
            fig, ax = plt.subplots()
    
            if np.size(_data_antA) > 1:
                plt.scatter(_data_antA[:,0], _data_antA[:,7], color = 'blue', label = 'antA')
            if np.size(_data_calH) > 1:
                plt.scatter(_data_calH[:,0], _data_calH[:,7], color = 'red', label = 'calH')
            if np.size(_data_antANd) > 1:
                plt.scatter(_data_antANd[:,0], _data_antANd[:,7], color = 'green', label = 'antA + nd')
            if np.size(_data_calHNd) > 1: 
                plt.scatter(_data_calHNd[:,0], _data_calHNd[:,7], color = 'purple', label = 'calH + nd')
                
            plt.title('measurement {}'.format(i + 1))
            plt.legend()
            plt.show()
            plt.close()
     
    # Plots all the data used if 'plotMeasuremetnsCombined' is chosen
    if plotMeasuremetnsCombined:
        fig, ax = plt.subplots()
        
        
        plt.scatter(data_antAt[:,0], data_antAt[:,1], color = '#069AF3', marker = ',', label = 'Antenna A')
        plt.scatter(data_calHt[:,0], data_calHt[:,1], color = 'red', marker = ',', label = 'CalH')
        plt.scatter(data_calHNdt[:,0], data_calHNdt[:,1], color = '#89fe05', marker = ',', label = 'CalH + Nd')
        
        plt.xlabel('Time [seconds since midnight]')
        plt.ylabel('Counts')
        
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
        
        plt.show()
        plt.close()
        
    print('\ndone\n')
    return _Vu, _measTime, _mean_calH, _tempUsed



def antennaCorrection(_lmbda, _m):
    # Calculates antenna correction parameters in horizontal and vertical direction
    # For open ended rectangular waveguide
    # Input: 
    # _lmbda  = wavelength [m]
    # _m      = coordinate to calculate over [rad]
    # Return: 
    # _vD     = Vertical correction parameter
    # _hD     = Horizontal correction parameter
    
    # Coordinate system setup
    _res   = 1000
    _theta = np.linspace(0.0000000000001, np.pi, num = _res)         
    _phi   = np.linspace(0.0000000000001, 2 * np.pi, num = _res)     

    _Theta, _Phi = np.meshgrid(_theta, _phi)
    _dtheta = np.abs(_theta[0] - _theta[1])
    _dphi   = np.abs(_phi[0] - _phi[1])
    
    # Field setup 
    _E0    = 1
    
    _a = 0.02286   # Antenna size, meters
    _b = 0.01016
    
    _k     = 2*np.pi / _lmbda   # Wavenumber
    
    
    # Antenna field parameters as described in Balanis Antenna Theory
    _X  = _k * _a / 2 * np.sin(_Theta) * np.cos(_Phi)
    _Y  = _k * _b / 2 * np.sin(_Theta) * np.sin(_Phi)
    _Cf = _a * _b * _k * _E0 /(2 * np.pi)
    
    _F_theta = _Cf/2 * np.sin(_Phi) * (1 + np.cos(_Theta)) * np.sin(_X)/_X * np.sin(_Y)/_Y
    _F_phi   = _Cf/2 * np.cos(_Phi) * (1 + np.cos(_Theta)) * np.sin(_X)/_X * np.sin(_Y)/_Y
    
    # Far-field pattern function
    _F = np.abs(np.sqrt(_F_theta**2 + _F_phi**2))
    
    # Directivity 
    _tmp =  np.sum(_F * np.sin(_Theta), axis = 1) * _dtheta
    _tmp =  np.sum(_tmp) * _dphi
    _D = 4 * np.pi * _F / _tmp
    
    
    # Calculate size of target in antenna field (correction for vertical line integration)
    _target_size = 2 * np.arctan(1/(2*np.sqrt(2))/9)
    _mask        = np.abs(_theta) < _target_size/2
    
    _D_ttarget   = np.sum(_D[:, _mask], axis = 1) * _dtheta
    
    _mask        = _theta < np.pi/2
    _D_ttotal    = np.sum(_D, axis = 1) * _dtheta
    
    _n = int(_res/4)
    
    _v_coeff = _D_ttarget/_D_ttotal
    _h_coeff = _D[_n,:]
    
    
    # Interpolate within reconstructed range, for horizontal and vertical direction
    _idx     = np.sum(np.abs(_phi) < np.max(_m)*1.1)
    _v_coeff = _v_coeff[:_idx]
    _v_coeff = np.append(np.flip(_v_coeff),_v_coeff)
    _v_coord = _phi[:_idx]
    _v_coord = np.append(-np.flip(_v_coord),_v_coord)
    
    _f = interp1d(_v_coord, _v_coeff)
    _vD = _f(_m)
    
    
    _idx     = np.sum(np.abs(_theta) < np.max(_m)*1.1)
    _h_coeff = _h_coeff[:_idx]
    _h_coeff = np.append( np.flip(_h_coeff),_h_coeff)
    _h_coord = _theta[:_idx]
    _h_coord = np.append(-np.flip(_h_coord),_h_coord)
    
    _f  = interp1d(_h_coord, _h_coeff)
    _hD = _f(_m)
    
    return _vD, _hD