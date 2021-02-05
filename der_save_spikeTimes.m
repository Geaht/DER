function [] = der_save_spikeTimes(spikeInfos)
%save_spikeTimes
%   saving CSC*_spikes.mat, times_CSC*.mat and cluster_info.mat from 
%   spikeInfos after deleting spikes detected multiple times

%   no_channels = length(unique(spikeInfos.channelID));
%   channels = get_channels_gd;
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.

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
    currChanName = regexp(header{find(strcmp('Name', header))+1},'\','split');
    chnname{chan} = currChanName{3}(1:end-4);
end
[~, idx] = sort(channels);
chnname = chnname(idx);


cluster_info = cell(3,no_channels);

chanPerBundle = max(unique(spikeInfos.channelID)) / max(unique(spikeInfos.bundleID));
localChanNo = 0;

% save new spike-times
for chan = 1:no_channels
    index_currSp = spikeInfos.channelID == chan;

    if localChanNo < chanPerBundle
        localChanNo = localChanNo + 1;
    else
        localChanNo = 1;
    end
    
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

        channame = spikeInfos.region(index_currSp);
        currBundle = spikeInfos.bundleID(index_currSp);
        
        channame = [channame{1} num2str(localChanNo)];
        cluster_info(2,chan) = {channame};
        
        % get current filenames
        spikesfile = sprintf('CSC%d_spikes.mat',chan);
        timesfile = sprintf('times_CSC%d.mat',chan);
        
        % save data of current channel
        save(spikesfile,'index_ts','spikes')
        save(timesfile,'spikes','cluster_class')
    end
end

% save cluster_info
label_info = '1 = MU 2 = SU-1 = Artif.Refers to "cluster_class"-values 1 and up.Ignores Unassigned (value 0)';
save cluster_info cluster_info label_info

end