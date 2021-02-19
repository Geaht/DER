# DER algorithm
Readme for the Duplicate Event Removal algorithm to implement the MATLAB code from:  
https://github.com/Geaht/DER

Version 1.0, released Februar 2021  

This project is licensed under the terms of the Mozilla Public License Version 2.0.  
Copyright (c) by Gert Dehnen, Marcel S. Kehl, Florian Mormann and the University of Bonn Medical Center  
This software was tested with MATLAB (R2018a) on Linux, Windows and MacOS X

------------------------------------------------------------------------------------------
The DER algorithm is written in MATLAB2018a. To run the code the following MATLAB packages are required:  

* Statistics and Machine Learning Toolbox
* Wavelet Toolbox  


## How to use the code 

Download the source code from this repository and add it along with the subdir to your MATLAB path. 
To run the DER algorithm for automated artefact detection you need to execute the DER.m in MATLAB.  
Please define the spike-sorting algorithm used ('Wave_clus' or 'Combinato') and the path to the data as:  

```
DER(data_path,'Wave_clus')  
```

The supported data structure are outputs of different cluster algorithms. Currently compatible are  
Combinato Spike Sorting and Wave_clus. In order to apply the DER algorithm to a different 
spike sorting algorithms you can adjust der_get_spikeInfos.m and der_save_spikeinfos.m to your data structure. 
Within this pipeline the data structure needed for the detection is created by the der_get_spikeInfos.m.  


## Structure of the code

* Part I - Detection of artifacts within diffrerent wire bundles  
The detection of artifacts across different bundles is done by '''der_detectArtifacts.m'''.  

* Part II - Detection of spikes within channels of same wire bundle  
The detection of artifacts within bundles is done by '''der_detectDuplicateSpikes'''.  

* Part III - Detection of suspicious cross-correlations  
Calculation of all cross-correlations is performed by '''der_cal_spike_cross_corr_mat.m '''. 
Suspicious central bins of the corss-correlograms are identified by '''der_detect_cross_corr_spikes.m'''.

Detected spike events are labeled individually for each part of  the detection pipeline using prime factor labels.
The detection label is given by the procuct of the following prime factors:

* 2 - detected across bundles (Part I)
* 3 - detected within the same channel (Part I)  
* 5 - detected within the same bundle (but different channel) (Part I)  
* 7 - detected based on the cross-correlograms (Part III)


## References

Cite as: TBA


