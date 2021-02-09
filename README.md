# DER 
Readme for the Duplicate Event Removal algorithm to implement the MATLAB code from:
https://github.com/Geaht/DER

Version 1.0, released Februar 2021
This project is licensed under the terms of the Mozilla Public License Version 2.0.

Copyright (c) by Gert Dehnen, Marcel S. Kehl and the University of Bonn Medical Center

This software is tested with MATLAB on Linux, Windows and MacOS ....

------------------------------------------------------------------------------------------
The DER algorithm is written in MATLAB2018a. Download the source code from this repository
is all you need to use the DER algorithm. To run the code two packages are needed:

Statistics and Machine Learning Toolbox
Wavelet Toolbox


To run the DER algorithm for the artefact detection you need to execute the DER.m in MATLAB.
Within this pipeline the data structure needed for the detection is created by the der_get_spikeInfos.m.
The supported data structure are Neuralynx binary files (NCS). Is you want to run the DER algorithm
on outputs from different cluster algorithms than the currently compatible onces (Combinato, WaveClus,
OSort), you only need to adjust the der_get_spikeInfos.m and der_save_spikeinfos.m to your data structure.

Part I - Detection of artifacts within diffrerent bundles
The detection of artifacts within different bundles is done by der_detectArtifacts.m.

Part II - Detection within the same channel of bundle
The detection of artifacts within different bundles is done by der_detectDuplicateSpikes.

Part III - Detection of suspicious cross-correlations
Analysing the central bin of all combinations of cross-correlations is done by
der_cal_spike_cross_corr_mat.m (preparing the data) and der_detect_cross_corr_spikes.m.

Every spike events that is detected by the DER algorithm is labeled individually for each part of
the detection pipeline. This was done multiplying prime numbers to the initial detection label of 1:

2: detected in Part I
3: detected within the same channel (Part I)
5: detected within the same bundle (but different channel) (Part I) 
7: detected in Part III


