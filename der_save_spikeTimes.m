function [] = der_save_spikeTimes(spikeInfos, clusterAlgorithm)
%save_spikeTimes
%   saving information of extracted and sorted spike shapes from 
%   spikeInfos after deleting spikes detected multiple times
%
%
%   Input: spikeInfos (table) containing the following information for each
%   spike:
%      region: hemisphere and region of each channel
%      bundleID: number of bundle
%      channelID: number of channel
%      threshold: threshold for spike-detection in µV
%      clusterID: number of cluster of current channel
%      unitClass: classification of unit in SU (single-), MU (multi-unit) or A
%       (artifact)
%      index_TS: time of amplitude of a spike in milliseconds
%      SpikeShapes: shape of spike (in out setup using combinato: 64 samples)
%   clusterAlgorithm: define the cluster-algorithm with that the data are 
%       preprocessed; Combinato, Wave_clus
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.


if ~exist('clusterAlgorithm','var')
    error('not enough input argumemts: "clusterAlgorithm" is missing')
end

% get list of bundle per channel
NCSFiles = struct2table(dir('*.ncs'));
no_channels = size(NCSFiles,1);

chnname = cell(no_channels,1);
channels = nan(no_channels,1);

for idx = 1:no_channels
    % get header infos of current channel
    currFile = NCSFiles.name{idx};
    channels(idx) = str2double(currFile(isstrprop(currFile,'digit')));

    fileID = fopen(currFile);
    header = textscan(fileID,'%s',55);
    header = header{1};
    fclose(fileID);       
    chnname{idx} = header{find(strcmp('-AcqEntName',header))+1};
end
[channels, idx] = sort(channels);
chnname = chnname(idx);

switch clusterAlgorithm
    case 'Combinato'
        cluster_info = cell(3,no_channels);
        cluster_info(2,:) = chnname;
        
        % save new spike-times
        
        % get template of times* and spikes* names
        timesFiles = dir('times_*.mat');
        [idx_start, idx_stop] = regexp(timesFiles(1).name, '[0-9]{1,3}');
        temp_timesName = [timesFiles(1).name(1:idx_start-1) '%i' timesFiles(1).name(idx_stop+1:end)];

        spikesFiles = dir('*_spikes.mat');
        [idx_start, idx_stop] = regexp(spikesFiles(1).name, '[0-9]{1,3}');
        temp_spikesName = [spikesFiles(1).name(1:idx_start-1) '%i' spikesFiles(1).name(idx_stop+1:end)];
        
        % loop over channels
        for idx = 1:no_channels
            chan = channels(idx);
            timesfile = sprintf(temp_timesName,chan);
            spikesfile = sprintf(temp_spikesName,chan);
            
            if exist(timesfile,'file') && exist(spikesfile,'file')
                
                index_currSp = spikeInfos.channelID == chan;

                if sum(index_currSp) > 0

                    % collect data of current CSC*_spikes.mat and times_CSC*.mat
                    spikes = spikeInfos.SpikeShapes(index_currSp,:); %#ok<*NASGU>
                    index_ts = spikeInfos.timeStamps(index_currSp);
                    clusterID = spikeInfos.clusterID(index_currSp);
                    cluster_class = [clusterID index_ts];

                    % restore cluster_info
                    [clusterIDList, index_first] = unique(clusterID);
                    unitClasses = spikeInfos.unitClass(index_currSp);
                    unitClasses = unitClasses(index_first);

                    cluster_info(3,chan) = {ones(1,size(clusterIDList(clusterIDList~=0),1))};

                    unitClasses = unitClasses(~contains(unitClasses,'A'));

                    unitID = [];
                    for unit = 1:length(unitClasses)
                        if contains(unitClasses(unit),'MU')
                            unitID = [unitID 1];
                        elseif contains(unitClasses(unit),'SU')
                            unitID = [unitID 2];
                        end
                    end

                    cluster_info(1,chan) = {unitID};

                    % save data of current channel
                    save(spikesfile,'index_ts','spikes')
                    save(timesfile,'spikes','cluster_class')
                    
                else % if no spike events remain in a channel after DER, previous files of spike events will be deleted
                    
                    delete(timesfile, spikesfile)
                    
                end
                
            elseif exist(timesfile,'file') || exist(spikesfile,'file')
                
                if exist(timesfile,'file')
                    error(['Missmatch: ' spikesfile ' is missing'])
                else
                    error(['Missmatch: ' timesfile ' is missing'])
                end
                
            end
        end

        % save cluster_info
        label_info = '1 = MU 2 = SU-1 = Artif.Refers to "cluster_class"-values 1 and up.Ignores Unassigned (value 0)';
        save cluster_info cluster_info label_info
    
    case 'Wave_clus'

        % save new spike-times
        
        % get template of times* and spikes* names
        timesFiles = dir('times_*.mat');
        [idx_start, idx_stop] = regexp(timesFiles(1).name, '[0-9]{1,3}');
        temp_timesName = [timesFiles(1).name(1:idx_start-1) '%i' timesFiles(1).name(idx_stop+1:end)];

        spikesFiles = dir('*_spikes.mat');
        [idx_start, idx_stop] = regexp(spikesFiles(1).name, '[0-9]{1,3}');
        temp_spikesName = [spikesFiles(1).name(1:idx_start-1) '%i' spikesFiles(1).name(idx_stop+1:end)];
        
        for idx = 1:no_channels
            
            timesfile = sprintf(temp_timesName,idx);
            spikesfile = sprintf(temp_spikesName,idx);
            
            if exist(timesfile,'file') && exist(spikesfile,'file')
                
                index_currSp = spikeInfos.channelID == idx;

                if sum(index_currSp) > 0

                    % collect data of current CSC*_spikes.mat and times_CSC*.mat
                    spikes = spikeInfos.SpikeShapes(index_currSp,:); %#ok<*NASGU>
                    index_ts = spikeInfos.timeStamps(index_currSp);
                    clusterID = spikeInfos.clusterID(index_currSp);
                    cluster_class = [clusterID index_ts];
                    index_ts = index_ts';

                    % get current filenames
                    spikesfile = sprintf('CSC%d_spikes.mat',idx);
                    timesfile = sprintf('times_CSC%d.mat',idx);

                    % save data of current channel
                    save(spikesfile,'index_ts','spikes','-append')
                    save(timesfile,'spikes','cluster_class','-append')
                    
                else % if no spike events remain in a channel after DER, previous files of spike events will be deleted
                    
                    delete(timesfile, spikesfile)
                    
                end
                
            else
                
                if exist(timesfile,'file')
                    error(['Missmatch: ' spikesfile ' is missing'])
                else
                    error(['Missmatch: ' timesfile ' is missing'])
                end
                
            end
        end  
end
