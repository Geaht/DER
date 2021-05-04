function [spikeInfos, EuclDis] = der_detectArtifacts(spikeInfos, threshold, deltaT, minNoEventsToCompaire)
%der_detectArtifacts
%	der_detectArtifacts deletes duplicate spike-shapes over different bundles
%
%	Input: 
%   spikeInfos (table) containing the following information for each spike:
%   bundleID: number of bundle
%   channelID: number of channel
%   clusterID: number of cluster of current channel
%   timeStamps: time of amplitude of a spike
%   SpikeShapes: shape of spike (in out setup using combinato: 64 samples)
%
%   threshold: threshold for median eucledian distance to detect duplicate spikes
%   deltaT: max time diff between two spikes to calculate their eucledian
%       distance of wavelet-decomposition
%
%
%   Output:
%   spikeInfos without detected duplicate spikes
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.

EuclDis = [];

% set parameter
if ~exist('threshold', 'var') || isempty(threshold)
    threshold = 14.4;
end

if ~exist('deltaT','var') || isempty(deltaT)
    deltaT = .05; % threshold of time-diff between spikes within different bundles in ms
%     gitterSamples = 1; % time-gitter as number of samples
end

% % set some variables
% compaire spike shapes from first to a chosen sample after amplitude;
% default is 40, see details in Dehnen & Kehl et al., 2020,....
no_samples = 40;

if ~exist('minNoEventsToCompaire','var') || isempty(minNoEventsToCompaire)
    minNoEventsToCompaire = 3;
end
minNoPairs = size(nchoosek([1:minNoEventsToCompaire],2),1);

% sort all spiketimes
if issorted(spikeInfos.timeStamps)
    TS = spikeInfos.timeStamps;
else
    [TS, idxsort] = sort(spikeInfos.timeStamps);
    spikeInfos = spikeInfos(idxsort,:);
    warning('Input not sorted according to time stamps!');
    fprintf('Sorting input spikes... \n')
end

% find possible duplicate spikes
DiffIndexTS = diff(TS);
duplicateSpikes = find(DiffIndexTS <= deltaT);
diffDuSp = diff(duplicateSpikes);


% loop over time-stamp-cluster
index_duplicateSpikes = [];
currStatus = 0;
for currSpike = 1:length(diffDuSp)
    index_DuSp = [];

    progress = round((currSpike / length(diffDuSp))*100);
    if progress > currStatus
        sprintf('%d%%',progress)
        currStatus = currStatus + 1;
    end
    
    if diffDuSp(currSpike) == 1
        ds = currSpike;
        while ds < size(diffDuSp,1) && diffDuSp(ds) == 1
            index_DuSp = [index_DuSp duplicateSpikes(ds)]; %#ok<*AGROW>
            ds = ds + 1;
            if ds > size(diffDuSp,1)
                break
            end
        end
        index_DuSp = [index_DuSp duplicateSpikes(ds) duplicateSpikes(ds)+1];

        % correct for time diff to first spike of current index_DuSp
        currTimeDiff = DiffIndexTS(index_DuSp(1:end-1));
        if sum(currTimeDiff) > deltaT
            currTD = 1;
            while sum(currTimeDiff(1:(currTD+1))) <= deltaT
                if currTD < length(currTimeDiff)
                    currTD = currTD + 1;
                else
                    break
                end
            end
            diffLength = length(currTimeDiff) - (currTD);
            index_DuSp = index_DuSp(1:end-diffLength);
        end


        index_chToCompare = nchoosek(index_DuSp,2);
    else
        index_DuSp = [duplicateSpikes(currSpike) duplicateSpikes(currSpike)+1];
        index_chToCompare = index_DuSp;
    end

    
    % just compaire within different bundles
    diffBundles = diff([spikeInfos.bundleID(index_chToCompare(:,1)) spikeInfos.bundleID(index_chToCompare(:,2))],1,2);
    index_chToCompare = index_chToCompare(diffBundles ~=0,:);

    if ~isempty(index_chToCompare) && size(index_chToCompare,1) >= minNoPairs
        no_spikesToCompaire = size(index_chToCompare,1);
        possDuplSpi1 = spikeInfos.SpikeShapes(index_chToCompare(:,1),1:no_samples);
        possDuplSpi2 = spikeInfos.SpikeShapes(index_chToCompare(:,2),1:no_samples);
        
        if ~isempty(possDuplSpi1)
            % shift the spike just in case of bipolar setting            
            index_sign = find(sign(possDuplSpi1(20)) ~=sign(possDuplSpi2(20)));
            possDuplSpi1(index_sign) = possDuplSpi1(index_sign)*-1;
        end


        euclDis = nan(no_spikesToCompaire,1);
        for spi = 1:no_spikesToCompaire

            spikeShapes = [possDuplSpi1(spi,:);possDuplSpi2(spi,:)];
            [coeff] = der_wavedec(spikeShapes);
            euclDis(spi) = sqrt(sum(diff(coeff).^2));

        end

        if median(euclDis) <= threshold
            index_duplicateSpikes = [index_duplicateSpikes index_DuSp];
        end
        EuclDis = [EuclDis; euclDis];
    end
end
index_duplicateSpikes = reshape(index_duplicateSpikes,[],1);
index_duplicateSpikes = unique(index_duplicateSpikes);
spikeInfos.detectionLabel(index_duplicateSpikes) = spikeInfos.detectionLabel(index_duplicateSpikes)*2;
end
