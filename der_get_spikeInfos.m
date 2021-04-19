function [spikeInfos] = der_get_spikeInfos(clusterAlgorithm, nr_chBundle)
%get_spikeInfos
%   get_spikeInfos collects list of region, bundleID, channelID, threshold,
%   clusterID, unitClass, spike-times and spike-shapes for each spike in a 
%   session. Adjusted to data structure from Neuralynx Inc.
%
%
%   Input:
%   clusterAlgorithm: define the cluster-algorithm with that the data are 
%       preprocessed; Combinato, Wave_clus
%   nr_chBundle: number of channels per bundle; default is 8
%
%   Output: spikeInfos (table) containing the following information for each
%   spike:
%   region: hemisphere and region of each channel
%   bundleID: number of bundle
%   channelID: number of channel
%   threshold: threshold for spike-detection in µV
%   clusterID: number of cluster of current channel
%   unitClass: classification of unit in SU (single-), MU (multi-unit) or A
%       (artifact)
%   index_TS: time of amplitude of a spike in milliseconds
%   SpikeShapes: shape of spike (in out setup using combinato: 64 samples)
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.
dbstop if error
if ~exist('clusterAlgorithm','var')
    error('not enough input argumemts: "clusterAlgorithm" is missing')
end

% set default number of channels per bundle
if ~exist('nr_chanPerBundle','var')
    nr_chBundle = 8;
end

% get list of bundle per channel
NCSFiles = struct2table(dir('*.ncs'));
no_channels = size(NCSFiles,1);

chnname = cell(no_channels,1);
channels = nan(no_channels,1);

for chan = 1:no_channels
    % get header infos of current channel
    currFile = NCSFiles.name{chan};
    channels(chan) = str2double(currFile(isstrprop(currFile,'digit')));

    fileID = fopen(currFile);
    header = textscan(fileID,'%s',55);
    header = header{1};
    fclose(fileID);
    chnname{chan} = header{find(strcmp('-AcqEntName',header))+1};
end
[channels, idx] = sort(channels);
chnname = chnname(idx);


% get bundle identity
chnno = cell2mat(cellfun(@(x) str2double(x(end)), chnname,'uniformoutput',false));
index_lastChPBd = find(chnno == nr_chBundle);

index_bundle = ones(no_channels,1);
currBundle = 1;
index_firstCh = 1;
for bndl = 1:length(index_lastChPBd)

    index_lastCh = index_lastChPBd(bndl);
    index_bundle(index_firstCh:index_lastCh) = currBundle;
    currBundle = currBundle + 1;
    index_firstCh = index_lastCh + 1;

end   


% get all spikeshapes, peaktimes and indices of bundle
region = [];
bundleID = [];
channelID = [];
threshold = [];
clusterID = [];
unitClass = [];
timeStamps = [];
SpikeShapes = [];


switch clusterAlgorithm
    case 'Combinato'
        load cluster_info cluster_info

        % loop over channels
        for chan = 1:no_channels
            timesfile = sprintf('times_CSC%d.mat',channels(chan));
            currBundle = index_bundle(chan);

            if exist(timesfile,'file')
                load(timesfile, 'cluster_class', 'spikes')
                timeStamps = [timeStamps; cluster_class(:,2)]; %#ok<*AGROW,*IDISVAR,*NODEF> % time-stamp of spike event amplitude (sample 20)
                SpikeShapes = [SpikeShapes; spikes]; % shape of individual spike events
                channelID = [channelID; ones(size(cluster_class,1),1) + (channels(chan)-1)]; % channel number of each individual spike event
                bundleID = [bundleID; zeros(size(cluster_class,1),1) + currBundle]; % bundle number of each individual spike event
                clusterID = [clusterID; cluster_class(:,1)]; % cluster number of each individual spike event

                currChnname = cell(size(spikes,1),1);
                currChnname(:) = {chnname{chan}(1:end-1)};
                region = [region; currChnname];

                % get unit-class for each spike
                cd(fullfile(sprintf('CSC%d',channels(chan))))
                thresholds = h5read(sprintf('data_CSC%d.h5', channels(chan)),'/thr');
                currThreshold = nanmedian(thresholds(3,:));
                cd ..
                currUnitID = cluster_info{1,chan};

                threshold = [threshold; zeros(size(cluster_class,1),1) + currThreshold];
                currUnitClasses = cell(size(cluster_class,1),1);
                for uc = 1:length(currUnitID)
                    index_currUC = cluster_class(:,1) == uc;
                    currUC = currUnitID(uc);

                    if currUC == -1
                        currUC = {'A'};
                    elseif currUC == 1
                        currUC = {'MU'};
                    elseif currUC == 2
                        currUC = {'SU'};
                    end

                    currUnitClasses(index_currUC) = currUC;

                end
                unitClass = [unitClass; currUnitClasses];

                % set unitClass for artifacts
                unclassified = cellfun(@isempty,unitClass,'Uniformoutput', false);
                index_unclassi = cell2mat(unclassified);
                unitClass(index_unclassi) = {'A'};

            end
        end
        
        
    case 'Wave_clus'
        % loop over channels
        for chan = 1:no_channels
            timesfile = sprintf('times_CSC%d.mat',channels(chan));
            currBundle = index_bundle(chan);

            if exist(timesfile,'file')
                load(timesfile, 'cluster_class', 'spikes')
                timeStamps = [timeStamps; cluster_class(:,2)]; %#ok<*AGROW,*IDISVAR,*NODEF> % time-stamp of spike event amplitude (sample 20)
                SpikeShapes = [SpikeShapes; spikes]; % shape of individual spike events
                channelID = [channelID; ones(size(cluster_class,1),1) + (channels(chan)-1)]; % channel number of each individual spike event
                bundleID = [bundleID; zeros(size(cluster_class,1),1) + currBundle]; % bundle number of each individual spike event
                clusterID = [clusterID; cluster_class(:,1)]; % cluster number of each individual spike event

                currChnname = cell(size(spikes,1),1);
                currChnname(:) = {chnname{chan}(1:end-1)};
                region = [region; currChnname];

                % get unit-class for each spike
                spikefile = sprintf('CSC%d_spikes.mat',channels(chan));
                load(spikefile,'threshold_all')
                currThreshold = nanmedian(threshold_all);
                currUnitID = unique(cluster_class(:,1));
        %             currUnitID = cluster_info(chan,:);

                threshold = [threshold; zeros(size(cluster_class,1),1) + currThreshold];
                currUnitClasses = cell(size(cluster_class,1),1);
                for uc = 1:length(currUnitID)
                    index_currUC = cluster_class(:,1) == uc;
  
                    currUC = {'MU'};
                    disp('DER Version 1.0: using WaveClus will only separate into multi-units and artifacts')

                    currUnitClasses(index_currUC) = currUC;

                end
                unitClass = [unitClass; currUnitClasses];

                % set unitClass for artifacts
                unclassified = cellfun(@isempty,unitClass,'Uniformoutput', false);
                index_unclassi = cell2mat(unclassified);
                unitClass(index_unclassi) = {'A'};

            end
        end        
end


detectionLabel = ones(size(timeStamps,1),1);
spikeInfos = table(region, bundleID, channelID, threshold, clusterID, unitClass, timeStamps, SpikeShapes, detectionLabel);
[~, idxsort] = sort(spikeInfos.timeStamps);
spikeInfos = spikeInfos(idxsort,:);

end

