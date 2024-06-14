This Repository contains all code and data collected in the bacholor thesis 'Improving Radiometric Measurements using
Interferometry'. The scripts have been written by Nikolaj Schejbel Nielsen and William Hedegaard Langvad for the data analysis related to the project. The scripts are made for pyhton 3.10. 
 
To make the scripts work as is, a folder containing a 'data' and 'scripts' folder should be made. All scripts should be placed in the scripts folder, and unzipped data directly in the data folder, but a simple change in the specified directory in 'dataRead.py' should be enough to make it work with any folder, as long as the scripts remain in the same folder. 


 
**Quick rundown of the scripts**:

'antennaDirectivity.py' calculates and plots the antenna pattern for an open ended rectangular waveguide antenna. 

'bandwidthFringePattern.py' plots the bandwith limited fringe pattern. 

'dataRead.py' is used for data analysis, and was used for all the figures in the results section. 

'parameterizedSynthesisAngular.py' is used for simulating the setup for the experiment. 

'phaseError.py' calculates the phase errors for the measurement setup. 

'simulatedSynthesis.py' is used for a simple interferometer simulation. 
'synthesisFunction.py' holds a lot of functionality for most of the other scripts, as multiple functions are defined here. 
