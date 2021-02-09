# DER
Duplicate Event Removal algorithm - Artifact detection in human single unit recordings

This project is licensed under the terms of the Mozilla Public License Version 2.0.


The DER algorithm is written in MATLAB2018a. Download the source code from this repository
is all you need to use the DER algorithm. To run the code you need two packages in MATLAB:

Statistics and Machine Learning Toolbox
Wavelet Toolbox


The supported data structure are Neuralynx binary files (NCS). Is you want to run the DER algorithm
on outputs from different cluster algorithms than the currently compatible onces (Combinato, WaveClus,
OSort), you only need to adjust the der_get_spikeInfos.m and der_save_spikeinfos.m to your data structure.
