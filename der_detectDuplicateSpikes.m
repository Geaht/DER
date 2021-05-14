function [spikeInfos, duplicateSpikes_neg, duplicateSpikes_sameSign, EuclDis] = der_detectDuplicateSpikes(spikeInfos, nr_chBundle, index_ampSpi, thresholdShape, minDiffSpike)
%der_detectDuplicateSpikes
%   der_detectDuplicateSpikes search for dublicate spikes in all channels 
%   of a bundle; duplicate spikes in the same channel (because of positive
%   and negative spike event extraction) are also deleted
%
%
%	Input:
%   spikeInfos (table) containing the following information for each spike:
%   bundleID: number of bundle
%   channelID: number of channel
%   threshold: threshold for spike-detection in µV
%   clusterID: number of cluster of current channel; 0 defines artifacts and
%   	unassigned spikes
%   unitClass: classification of unit in SU (single-), MU (multi-unit) or A
%   	(artifact)
%   timeStamps: time of amplitude of a spike in milliseconds
%   SpikeShapes: shape of spike (in out setup using combinato: 64 samples)
%
%   nr_chBundle: number of channels per bundle; default is 8
%   thresholdRho: threshold for median eucledian distance to detect duplicate spikes
%
%   Output:
%   corrected spikeInfos
%   duplicateSpikes_neg
%   duplicateSpikes_sameSign
%   EuclDis
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.

dbstop if error

if ~exist('nr_chBundle','var') || isempty(nr_chBundle)
    nr_chBundle = 8;
end

if ~exist('index_ampSpi','var') || isempty(index_ampSpi)
    index_ampSpi = 20;
end

if ~exist('thresholdShape', 'var') || isempty(thresholdShape)
    thresholdShape = 8.5;
end

% compaire spike shapes von first to a chosen sample after amplitude;
% default is 40, see details in Dehnen & Kehl et al., 2020,....
no_samples = 40;

% sort all spiketimes
if ~issorted(spikeInfos.timeStamps)
    [~, idxsort] = sort(spikeInfos.timeStamps);
    spikeInfos = spikeInfos(idxsort,:);
    warning('Input not sorted according to time stamps!');
    fprintf('Sorting input spikes... \n')
end

% set some variables
nr_channels = max(spikeInfos.channelID);

duplicateSpikes_neg = 0;

% diff of timestamps in same channel must be below 
minDiffSpikeSC = 0.65;

% min difference of times of amplitudes in different channels in �s
if ~exist('minDiffSpike','var') || isempty(minDiffSpike)
    minDiffSpike = .05;
end

index_doubleSpikesAll = [];
EuclDis = [];
% loop over bundle
for channel = 1:nr_chBundle:nr_channels
    disp(num2str(channel))
    % find all existing channels in current bundle

    index = 1;
    index_ChBundle = [];

    for ch = channel:channel+(nr_chBundle-1)
        
        index_currSp = find(spikeInfos.channelID == ch);
    
        if sum(index_currSp) > 0
            spikes = spikeInfos.SpikeShapes(index_currSp,:);
            timeStamps = spikeInfos.timeStamps(index_currSp);
            clusterID = spikeInfos.clusterID(index_currSp);
            
            index_ChBundle(index) = ch; %#ok<*AGROW>
            index = index + 1;
    
            %%%%%%%%%%%%%%% compare spikes in same channel %%%%%%%%%%%%%%%%
            % delete duplikate spikes compairing positive and negative spikes
            % in the same channel
        
            % diff of timestamps must be below minDiffSpikeSC
            index_duplicateSpikes = find(abs(diff(timeStamps)) <= minDiffSpikeSC);
            index_duplicateSpikes = [index_duplicateSpikes, index_duplicateSpikes+1];

            
            if ~isempty(index_duplicateSpikes)
                % compared spikes must have opposite signs to be the same spike
                % (this is a necessary criterium because each double detected
                % spike is one time extracted as a negative and one time as a
                % positive spike)
                signesMaxima = [sign(spikes(index_duplicateSpikes(:,1),index_ampSpi)), sign(spikes(index_duplicateSpikes(:,2),index_ampSpi))];
                index_oppositeSign = abs(diff(signesMaxima,[],2)) == 2;
                               
                index_duplicateSpikes = index_duplicateSpikes(index_oppositeSign,:);

                % count duplicate spikes
                duplicateSpikes_neg = duplicateSpikes_neg + size(index_duplicateSpikes,1);

                % decide which spike should survive
                cluster_nr = clusterID(index_duplicateSpikes(:,1));
                cluster_nr2 = clusterID(index_duplicateSpikes(:,2));

%                 index_artifactsToDelete = sort(reshape(index_duplicateSpikes(cluster_nr == 0 | cluster_nr2 == 0,:),[],1));
                index_artifactsToDelete = sort([index_duplicateSpikes(cluster_nr == 0,1); index_duplicateSpikes(cluster_nr2 == 0,2)]);
                index_duplicateSpikes(cluster_nr == 0 | cluster_nr2 == 0,:) = [];


                % kill all artifacts, if spike which is compared to is not in
                % artifact cluster
                timeStamps(index_artifactsToDelete) = nan;
                cluster_nr = [clusterID(index_duplicateSpikes(:,1)) clusterID(index_duplicateSpikes(:,2))];

                % delete spike of MU, if same spike is detected in SU and MU;
                % else delete spike of cluster with smaller signal to noise
                unitClass = cell(2,1);
                for spi = 1:size(index_duplicateSpikes,1)

                    unitClass(1) = spikeInfos.unitClass(index_currSp(index_duplicateSpikes(spi,1)));
                    unitClass(2) = spikeInfos.unitClass(index_currSp(index_duplicateSpikes(spi,2)));

                    % if unit class is equal keep spike of cluster with better
                    % signal / noise
                    if isequal(unitClass(1),unitClass(2))

                        index_ch = clusterID == cluster_nr(spi,1);
                        amplitude = nanmedian(spikes(index_ch,index_ampSpi));
                        threshold = spikeInfos.threshold(index_currSp(1));
                        signalToNoiseClA = abs(amplitude) / threshold;

                        index_ch = clusterID == cluster_nr(spi,2);
                        amplitude = nanmedian(spikes(index_ch,index_ampSpi));
                        threshold = spikeInfos.threshold(index_currSp(1));
                        signalToNoiseClB = abs(amplitude) / threshold;

                        if signalToNoiseClA <= signalToNoiseClB
                            timeStamps(index_duplicateSpikes(spi,1)) = nan;
                        else
                            timeStamps(index_duplicateSpikes(spi,2)) = nan;
                        end

                    elseif strcmp(unitClass{1},'SU')
                        timeStamps(index_duplicateSpikes(spi,2)) = nan;
                    else
                        timeStamps(index_duplicateSpikes(spi,1)) = nan;
                    end              
                end

                % delete duplicate spikes
                index_DS = isnan(timeStamps);
                index_duSp = index_currSp(index_DS);
                spikeInfos.detectionLabel(index_duSp) = spikeInfos.detectionLabel(index_duSp)*3;
            end
        end
    end
    
    
    %%%%%%%%%%%%%% compare spikes in different channels %%%%%%%%%%%%%%%%%%%
    % all possible combinations of existing channels
    if length(index_ChBundle) > 1
        
        index_currSpikes = ismember(spikeInfos.channelID,index_ChBundle);
        currSpikeInfos = spikeInfos(index_currSpikes,:);
        DiffIndexTS = diff(currSpikeInfos.timeStamps);
        duplicateSpikes = find(DiffIndexTS <= minDiffSpike);
        diffDuSp = diff(duplicateSpikes);
        
        index_duplicateSpikes = [];
        for currSpike = 1:length(diffDuSp)
            index_DuSp = [];
            
            if diffDuSp(currSpike) == 1
                ds = currSpike;
                while diffDuSp(ds) == 1 && ds < size(diffDuSp,1)
                    index_DuSp = [index_DuSp duplicateSpikes(ds)]; %#ok<*AGROW>
                    ds = ds + 1;
                    if ds > size(diffDuSp,1)
                        break
                    end
                end
                index_DuSp = [index_DuSp duplicateSpikes(ds) duplicateSpikes(ds)+1];
                
                % corr for time diff to first spike of current index_DuSp
                currTimeDiff = DiffIndexTS(index_DuSp(1:end-1));
                if sum(currTimeDiff) > minDiffSpike
                    currTD = 1;
                    while sum(currTimeDiff(1:(currTD+1))) <= minDiffSpike
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

            % just compaire within different channels
            diffChannels = diff([currSpikeInfos.channelID(index_chToCompare(:,1)) currSpikeInfos.channelID(index_chToCompare(:,2))],1,2);
            index_chToCompare = index_chToCompare(diffChannels ~=0,:);

            if ~isempty(index_chToCompare)
%                 no_spikesToCompaire = size(index_chToCompare,1);
                possDuplSpi1 = currSpikeInfos.SpikeShapes(index_chToCompare(:,1),1:no_samples);
                possDuplSpi2 = currSpikeInfos.SpikeShapes(index_chToCompare(:,2),1:no_samples);


                spikeShapes = [possDuplSpi1; possDuplSpi2];
                [inspk, inputs] = der_wavedec(spikeShapes);
                no_spikes = size(inspk,1);
                coefficients = nan(no_spikes/2, inputs, 2);
                coefficients(:,:,1) = inspk(1:no_spikes/2,:);
                coefficients(:,:,2) = inspk(no_spikes/2+1:end,:);

                % calculating euclidean distance
                euclDis = sqrt(sum((diff(coefficients,[],3)).^2,2));

                EuclDis = [EuclDis; euclDis];

                if median(euclDis) <= thresholdShape
                    index_duplicateSpikes = [index_duplicateSpikes; index_chToCompare];
                end
            end
        end
        
        for spi = 1:size(index_duplicateSpikes,1)

            unitClass(1) = currSpikeInfos.unitClass(index_duplicateSpikes(spi,1));
            unitClass(2) = currSpikeInfos.unitClass(index_duplicateSpikes(spi,2));

            % if unit class is equal keep spike of channel with larger
            % signal / noise
            % if one unit is labeled as an artifact both spikes
            % will be marked to be deleted
            if isequal(unitClass(1),'A') || isequal(unitClass(2),'A')
                
            elseif isequal(unitClass(1),unitClass(2))

                index_ch = currSpikeInfos.channelID(index_duplicateSpikes(spi,1));
                index_ch = currSpikeInfos.channelID == index_ch;
                index_cl = currSpikeInfos.clusterID(index_duplicateSpikes(spi,1));
                index_cl = currSpikeInfos.clusterID(index_ch) == index_cl;
                amplitude = nanmedian(currSpikeInfos.SpikeShapes(index_cl,index_ampSpi));
                threshold = currSpikeInfos.threshold(index_duplicateSpikes(spi,1));
                signalToNoiseChA = abs(amplitude) / threshold;
                
                index_ch = currSpikeInfos.channelID(index_duplicateSpikes(spi,2));
                index_ch = currSpikeInfos.channelID == index_ch;
                index_cl = currSpikeInfos.clusterID(index_duplicateSpikes(spi,2));
                index_cl = currSpikeInfos.clusterID(index_ch) == index_cl;
                amplitude = nanmedian(currSpikeInfos.SpikeShapes(index_cl,index_ampSpi));
                threshold = currSpikeInfos.threshold(index_duplicateSpikes(spi,2));
                signalToNoiseChB = abs(amplitude) / threshold;
                
                if signalToNoiseChA <= signalToNoiseChB
                    index_duplicateSpikes(spi,2) = nan;
                else
                    index_duplicateSpikes(spi,1) = nan;
                end

            elseif strcmp(unitClass{1},'SU')
                index_duplicateSpikes(spi,2) = nan;
            else
                index_duplicateSpikes(spi,1) = nan;
            end
        end            

        % delete duplicate spikes
        index_duplicateSpikes = reshape(index_duplicateSpikes,[],1);
        index_currSpikes = find(index_currSpikes~=0);
        index_doubleSpikesAll = [index_doubleSpikesAll; index_currSpikes(index_duplicateSpikes(~isnan(index_duplicateSpikes)))];

    end
end

index_doubleSpikesAll = unique(index_doubleSpikesAll);
duplicateSpikes_sameSign = length(index_doubleSpikesAll);
spikeInfos.detectionLabel(index_doubleSpikesAll) = spikeInfos.detectionLabel(index_doubleSpikesAll)*5;
end