function [] = DER(dataPath, clusterAlgorithm, save_structures)
%DER algorithm (Duplicate Event Removal)
%	DER algorithm (Duplicate Event Removal) identifies spike events
%   recorded several times (for details see Dehnen, Kehl et al.: Duplicate
%   detection of spike events: A relevand problem in human single-unit
%   recordings, Brain Science 2021)
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.


if exist('dataPath','var') && ~isempty(dataPath)
    cd(dataPath)
end
if ~exist('save_structures','var')
    save_structures = true;
end


% collects list of region, bundleID, channelID, threshold, clusterID, unitClass, 
% spike-times and spike-shapes for each spike in a session; select the
% cluster algorithm to Combianto, Wave_clus
[spikeInfos] = der_get_spikeInfos(clusterAlgorithm);


%% Part I: spike events within different bundles
% deletes duplicate spike-shapes over different bundles
[spikeInfos] = der_detectArtifacts(spikeInfos);

%% Part II: spike events within the same bundle
% search for dublicate spikes in all channels of a bundle; duplicate spikes
% in the same channel (because of positive and negative clustering) are
% also labeled
[spikeInfos] = der_detectDuplicateSpikes(spikeInfos);

%% Part III: cross-correlation
% Analysing the central bin of cross-correlations of clusters in different 
% channels (wires). If the central bin exceeds a threshold of z = 5, these
% spike events are labeled as artifacts

% first calculate a matrix containing all cross-correlations
[cross_corr_mat_info] = der_cal_spike_cross_corr_mat(spikeInfos);

[spikeInfos] = der_detect_cross_corr_spikes(spikeInfos, cross_corr_mat_info);


%% Saving data
if save_structures
    save spikeInfos spikeInfos
end
% deleting spike events labeled as artificial
spikeInfos=spikeInfos(spikeInfos.detectionLabel == 1,:);

% saving CSC*_spikes.mat and times_CSC*.mat from spikeInfos after deleting
% spikes detected multiple times 
der_save_spikeTimes(spikeInfos, clusterAlgorithm);

end
