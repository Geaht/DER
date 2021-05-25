function [spikeInfos,portion_of_spikes_found ] = der_detect_cross_corr_spikes(spikeInfos,cross_corr_mat_info, ...
                                 z_thr,option_for_detection,index_ampSpi,session_name,do_plots,outputpath)
%%  [spikeInfos,portion_of_spikes_found ] = dds_detect_cross_corr_spikes(spikeInfos,cross_corr_mat_info, ...
%                z_thr,option_for_detection,index_ampSpi,session_name,do_plots,outputpath)
%
% Function for identifieing and labeling spikes occuring within the central bin 
% of a highly significant zero-lag cross-correlation
% 
% Inputs:
%
% spikeInfos: Table containing spike information
%             (output of the function get_spikeInfos)
%
% cross_corr_mat_info: Matrix containing the cross-correlations of all cluster combination in the
%          session (output of the function dds_cal_spike_cross_corr_mat)
%
% z_thr: Threshold value of the spike-count in the central time-bin in the 
%        z-normalized Cross-Correlation matrix. Only for counts larger than z_thr, spikes 
%        spikes will be deleted. (default value 20)
%
% option_for_detection: 1 --> detects significant zero lag Cross-Correlations 
%                             only within the same bundle on diferent channels
%                       2 --> detects significant zero lag Cross-Correlations over all channels
%
% index_ampSpi: Index of time bin to which spike waveforms are alligend.
%               (default 20 for Combinato)
%
% do_plots: If do_plots==1 Cross-Correlation plots of the units exceeding the z-value 
%           threshold will be ploted
%
% outputpath: Path where data output plots and data shall be stored.
%             (default current directory)
%
% Output:
%
% spikeInfos: Table containing spike information with spikes marked for
%             deletion (detection label 7)            
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.


%% Input data
if ~exist('z_thr','var') || isempty(z_thr)
    z_thr=5; % to be finalized       
end
% index of spike wave max
if ~exist('index_ampSpi','var') || isempty(index_ampSpi)
    index_ampSpi = 20;
end

if ~exist('outputpath','var') || isempty(outputpath)
    outputpath=pwd; 
end

if ~exist('do_plots','var') || isempty(do_plots)
    do_plots=0; 
end

if ~exist('delete_artifacts','var') || isempty(delete_artifacts)
    delete_artifacts=1;
end

if ~exist('session_name','var')  || isempty(session_name)
    session_name='temp'; 
end

if ~exist('cross_corr_mat_info','var')  || isempty(cross_corr_mat_info)
    error('No cross-correlation matrix provided');
else 
    % get data from Cross-Correlation calculation
    cross_corr_mat=cross_corr_mat_info.cross_corr_mat;
    spikes_in_central_bin=cross_corr_mat_info.spikes_in_central_bin;
    bin_width=cross_corr_mat_info.bin_width;
    hist_limit=cross_corr_mat_info.hist_limit;
end

if  ~exist('option_for_detection','var')
    option_for_detection=2;
end

% In casere there are no input data are provided, 
% try to load spikes data from current directory

if  ~exist('spikeInfos','var')
    
    try
        fprintf('Loading spikes from current folder! \n')
        spike_info_files=dir('spikeInfos*.mat');
        if numel(spike_info_files)
            load(spike_info_files(1).names);
        else
            error('More than one spikeInfo file found!');
        end
    catch
        error('No spike data profided or found');
    end
end

% Check if the spike data are already sorted acording to spike times, if
% not sort them.
if ~issorted(spikeInfos.timeStamps)
    [~,idxsort] = sort(spikeInfos.timeStamps);
    spikeInfos=spikeInfos(idxsort,:);
    warning('Input not sorted according to spike times!');
    fprintf('Sorting input spikes... \n')
end

% get list of all cluster in the session
cluster_list=unique([spikeInfos.channelID spikeInfos.clusterID],'rows');
central_bin_idx = floor(size(cross_corr_mat,3)/2)+1;

%% detect spikes with high zero-lag cross-correlation 

% get all bunlde and channels
all_bundleIDs=unique(spikeInfos.bundleID);
% svae idx of spikes to delete
spikes_to_delete=zeros(numel(spikeInfos.timeStamps),1);


%% option 1: Detect spikes with high zero-lag cross-correlations in same bundle
if option_for_detection==1
    
    tic
    % loop over all bundels
    for idx_bundle=1:numel(all_bundleIDs)

        fprintf('Analysing bunlde %i \n',all_bundleIDs(idx_bundle));  
        
        % get channels for this bundle 
        temp_channels_per_bundle=unique(spikeInfos.channelID( ...
                            spikeInfos.bundleID==all_bundleIDs(idx_bundle)));  

        % generate all possible pairs of two channels with in this bundle  
        temp_channel_pairs=nchoosek(temp_channels_per_bundle,2);
        
        % if there is only one channel, continue
        if numel(temp_channels_per_bundle)==1
            continue
        end

        % otherwise loop over all channel pairs
        for channel_pair_idx = 1:length(temp_channel_pairs(:,1))

            ch1=temp_channel_pairs(channel_pair_idx,1);
            ch2=temp_channel_pairs(channel_pair_idx,2);

            % get all clusters for booth channels
            temp_clusters_ch1=unique(spikeInfos.clusterID( ...
                                        spikeInfos.channelID==ch1));
            temp_clusters_ch2=unique(spikeInfos.clusterID( ...
                                        spikeInfos.channelID==ch2));


            % generate all possible pairs of two cluster booth channels
            temp_cluster_pairs=combvec(temp_clusters_ch1(:)', ...
                                  temp_clusters_ch2(:)')';    

            % loop over pairs of clusters that are compared
            for idx_clus_pair=1:length(temp_cluster_pairs(:,1))

                % get current cluster numbers to compare
                clus_nr1=temp_cluster_pairs(idx_clus_pair,1);
                clus_nr2=temp_cluster_pairs(idx_clus_pair,2);

                % find correct row in the cross correlation Matrix                
                cluster_number_1 = ismember(cluster_list,[ch1 clus_nr1],'rows');
                cluster_number_2 = ismember(cluster_list,[ch2 clus_nr2],'rows');

                % get bincounts from matrix
                bincounts = squeeze(cross_corr_mat(cluster_number_1,cluster_number_2,:))';

                % calculate number spikes within 1 ms (or 1 bin_width) in booth cluster
                N_central_bin = bincounts(central_bin_idx);

                % mean of all bins except the central bin
                % may consider using median_bincount=median(bincounts(1:end ~= central_bin_idx));
                mean_bincount = mean(bincounts(1:end ~= central_bin_idx));
                std_bincount  =  std(bincounts(1:end ~= central_bin_idx));

                % In case there are no spikes, go to next cluster
                if std_bincount==0 
                    z_central_bin=0;
                else
                    z_central_bin=(N_central_bin-mean_bincount)/std_bincount;
                end
                
                
                %% if z-value of central bin exceeds the threshold, 
                %  find spikes to delete
                if  z_central_bin>z_thr

                    % check if clusts corespond to single unit or multi unit
                    clus_kind_1=unique(spikeInfos.unitClass(spikeInfos.channelID==ch1 & ...
                                       spikeInfos.clusterID==clus_nr1));
                    clus_kind_2=unique(spikeInfos.unitClass(spikeInfos.channelID==ch2 & ...
                                       spikeInfos.clusterID==clus_nr2));

                    idx_clus_1=spikeInfos.channelID==ch1 & ...
                               spikeInfos.clusterID==clus_nr1;

                    idx_clus_2=spikeInfos.channelID==ch2 & ...
                               spikeInfos.clusterID==clus_nr2;                

                    % if one of the clusters is marked as an Artifact and
                    % the option delete_artifacts is set to 1, central spikes
                    % within booth clusters will be deleted
                    if (strcmp(clus_kind_1,'A')||strcmp('A',clus_kind_2)) && ...
                            delete_artifacts==1

                        % loop over spikes and label spikes for deletion
                        A1_delete_idx = spikeInfos.channelID==ch1 & ...
                                        spikeInfos.clusterID==clus_nr1  & ...
                                        spikes_in_central_bin;

                        A2_delete_idx = spikeInfos.channelID==ch2 & ...
                                        spikeInfos.clusterID==clus_nr2  & ...
                                        spikes_in_central_bin;

                        % identify spikes in that are in the central bin of 
                        % the two selected cluster
                        A1_spike_times=spikeInfos.timeStamps(A1_delete_idx)';
                        A2_spike_times=spikeInfos.timeStamps(A2_delete_idx)';
                        A1_delete_idx_numbers=find(A1_delete_idx);
                        A2_delete_idx_numbers=find(A2_delete_idx);
                        
                        A2_delete_idx_new=zeros(size(A2_delete_idx));
                        
                        % loop over spikes in the first cluster
                        for idx1=1:numel(A1_spike_times)
                            
                             % find spikes in the secound cluster that are
                             % within the bin-width
                             temp_idx_A2_delete = find(abs(A1_spike_times(idx1) - A2_spike_times)< bin_width/2);
                            
                             % if there are spikes of cluster 2
                             % in the same time bin delete them
                             if numel(temp_idx_A2_delete) > 0
                                  A2_delete_idx_new(A2_delete_idx_numbers(temp_idx_A2_delete))=1;
                                  
                             % if there is no spike of clus 2 in the same
                             % time-bin do not delete the spike
                             else
                                  A1_delete_idx(A1_delete_idx_numbers(idx1)) = 0;
                             end
                        end
                       
                       % save for late which spikes should be deleted
                       spikes_to_delete(A1_delete_idx)=spikes_to_delete(A1_delete_idx)+1;
                       spikes_to_delete(logical(A2_delete_idx_new))=spikes_to_delete(logical(A2_delete_idx_new))+1;

                       % sanity check
                       temp_delta_t= pdist2(spikeInfos.timeStamps(logical(A2_delete_idx_new)),...
                           spikeInfos.timeStamps(logical(A1_delete_idx)));
                       
                       assert(sum(sum(temp_delta_t < bin_width/2)) == N_central_bin,...
                       'Number of spikes to delete do not match the cross-correlations!')
                       
                       if N_central_bin < sum(A1_delete_idx) ||  N_central_bin < sum(A2_delete_idx_new) 
                           warning('Number of spikes to delete do not match the cross-correlations!');
                       end
                       
                       fprintf('Deleting %i artefact spikes in channel %i cluster %i! \n', ...
                            sum(A1_delete_idx),ch1,clus_nr1)
                       fprintf('Deleting %i artefact spikes in channel %i cluster %i! \n', ...
                            sum(A2_delete_idx_new),ch2,clus_nr2)
                    end

                    % If a SU and a MU are compared the dublicated spikes 
                    % in the central bin of the MU are deleted
                    
                    if (strcmp(clus_kind_1,'MU') && strcmp(clus_kind_2,'SU'))||...
                        (strcmp(clus_kind_1,'SU') && strcmp(clus_kind_2,'MU'))
                    
                        % Find spikes in the central bin
                        % These are all spikes within the time bin_width
                        if strcmp(clus_kind_1,'MU') && strcmp(clus_kind_2,'SU')
                            ch_MU=ch1; clus_MU=clus_nr1; ch_SU=ch2; clus_SU=clus_nr2; 
                        elseif strcmp(clus_kind_2,'MU') && strcmp(clus_kind_1,'SU')
                            ch_MU=ch2; clus_MU=clus_nr2; ch_SU=ch1; clus_SU=clus_nr1; 
                        end

                        % loop over spikes and delete MU spikes 
                        MU_delete_idx = spikeInfos.channelID==ch_MU    & ...
                                        spikeInfos.clusterID==clus_MU  & ...
                                        spikes_in_central_bin;
                                    
                        SU_keep_idx   = spikeInfos.channelID==ch_SU    & ...
                                        spikeInfos.clusterID==clus_SU  & ...
                                        spikes_in_central_bin;

                        MU_spike_times = spikeInfos.timeStamps(MU_delete_idx)';
                        SU_spike_times = spikeInfos.timeStamps(SU_keep_idx)';
                        MU_delete_idx_numbers = find(MU_delete_idx);
                        
                        % loop over spikes of the MU and check if they are within
                        % the same time bin as spikes from the SU 
                        for idxMU = 1:numel(MU_spike_times)
                            
                            temp_min_time_diff=min(abs(MU_spike_times(idxMU) - ...
                                                             SU_spike_times));
                             if temp_min_time_diff >= bin_width/2
                                  MU_delete_idx(MU_delete_idx_numbers(idxMU)) = 0;
                             end
                        end

                        if sum(MU_delete_idx) > N_central_bin
                           error('Number of spikes to delete not matches cross-correlations!');
                        end
                        % save for late which spikes should be deleted
                        spikes_to_delete(MU_delete_idx)=spikes_to_delete(MU_delete_idx)+1;
                        fprintf('Deleting %i Spikes in channel %i cluster %i! \n', ...
                            sum(MU_delete_idx),ch_MU,clus_MU)

                    % If booth clusters are MU or booth SU, delete spikes
                    % in the cluster with a lower signal / noise 
                    elseif (strcmp(clus_kind_1,'SU') && strcmp(clus_kind_2,'SU')) || ...
                           (strcmp(clus_kind_1,'MU') && strcmp(clus_kind_2,'MU'))
                       
                        % calculate the median amplitude of the clusters 
                        amplitude_clus_1=nanmedian(spikeInfos.SpikeShapes(idx_clus_1,index_ampSpi));
                        amplitude_clus_2=nanmedian(spikeInfos.SpikeShapes(idx_clus_2,index_ampSpi));

                        % get the noise level of the clusters
                        noise_clus_1=spikeInfos.threshold(find(idx_clus_1,1,'first'));
                        noise_clus_2=spikeInfos.threshold(find(idx_clus_2,1,'first'));

                        signal_to_noise_clus_1 = abs(amplitude_clus_1)/noise_clus_1;
                        signal_to_noise_clus_2 = abs(amplitude_clus_2)/noise_clus_2;

                        if signal_to_noise_clus_1 < signal_to_noise_clus_2
                            ch_Low=ch1; clus_Low=clus_nr1; ch_High=ch2; clus_High=clus_nr2; 
                        elseif signal_to_noise_clus_2 < signal_to_noise_clus_1
                            ch_Low=ch2; clus_Low=clus_nr2; ch_High=ch1; clus_High=clus_nr1; 
                        end                           

                        % loop over spikes and label spikes in the cluster
                        % with lower S/N
                        Low_delete_idx = spikeInfos.channelID==ch_Low & ...
                                         spikeInfos.clusterID==clus_Low& ...
                                         spikes_in_central_bin;

                        High_keep_idx = spikeInfos.channelID==ch_High & ...
                                        spikeInfos.clusterID==clus_High  & ...
                                        spikes_in_central_bin;

                        Low_spike_times  = spikeInfos.timeStamps(Low_delete_idx)';
                        High_spike_times = spikeInfos.timeStamps(High_keep_idx)';
                        Low_delete_idx_numbers = find(Low_delete_idx);

                        % loop over spikes of the lower S/N unit and check
                        % if they are spikes of the other unit with the
                        % bin_width
                        for idxLow=1:numel(Low_spike_times)
                            temp_min_time_diff=min(abs(Low_spike_times(idxLow) - ...
                                                             High_spike_times));
                             if temp_min_time_diff >= bin_width/2
                                  Low_delete_idx(Low_delete_idx_numbers(idxLow)) = 0;
                             end
                        end

                        if  sum(Low_delete_idx) > N_central_bin
                            warning('Number of spikes to delete not matches cross_corr mat');
                        end

                        % save for late which spikes should be deleted
                        spikes_to_delete(Low_delete_idx)=spikes_to_delete(Low_delete_idx)+1;        
                        fprintf('Deleting %i Spikes in channel %i cluster %i! \n', ...
                                sum(Low_delete_idx),ch_Low,clus_Low)
                    end

                    % plot and save cross-correlogramm
                    if do_plots
                        fprintf('z-Score of the central bin %f (chan %i and %i) \n', ...
                                     z_central_bin,ch1,ch2);
                        der_plot_cross_corr_example(spikeInfos,bincounts,bin_width,hist_limit,outputpath,...
                                           session_name,ch1,ch2,clus_nr1,clus_nr2,...
                                           idx_clus_1,idx_clus_2,clus_kind_1,clus_kind_2); 
                    end
                end
            end

        end
    end



%% option 2: Delete spikes with high zentral bin in the cross-correlogramm 
%  for all channel combinations

elseif option_for_detection == 2
    
    % get all possible clustercombinations
    N_cluster=size(cluster_list,1);
    cluster_combinations =  nchoosek(1:N_cluster,2); 
    N_dist_pairs=size(cluster_combinations,1);
    
    % add pairs with same cluster
    for ii=1:N_cluster
        cluster_combinations(N_dist_pairs+ii,1)=ii;
        cluster_combinations(N_dist_pairs+ii,2)=ii;
    end

    tic
    % loop over all pairs
    for idx_clus_pair=1:size(cluster_combinations,1)

        % get current cluster numbers to compare
        ch1 = cluster_list(cluster_combinations(idx_clus_pair,1),1);
        ch2 = cluster_list(cluster_combinations(idx_clus_pair,2),1);

        clus_nr1 = cluster_list(cluster_combinations(idx_clus_pair,1),2);
        clus_nr2 = cluster_list(cluster_combinations(idx_clus_pair,2),2);

        % find correct row in the cross correlation Matrix                
        cluster_number_1 = ismember(cluster_list,[ch1 clus_nr1],'rows');
        cluster_number_2 = ismember(cluster_list,[ch2 clus_nr2],'rows');
                              
        % get bincounts from matrix
        bincounts=squeeze(cross_corr_mat(cluster_number_1,cluster_number_2,:))';

        % calculate number spikes within 1 ms (or 1 bin_width) in booth cluster
        N_central_bin=bincounts(central_bin_idx);

        % mean of all bins except the central one
        % or use median_bincount=median(bincounts(1:end ~= central_bin_idx));
        mean_bincount = mean(bincounts(1:end ~= central_bin_idx));
        std_bincount  = std(bincounts(1:end ~= central_bin_idx));

        % In case there are no spikes, go to next cluster
        if std_bincount==0
            z_central_bin=0;
        else
            z_central_bin=(N_central_bin-mean_bincount)/std_bincount;
        end

        %% if z-value of central bin exceeds the threshold, 
        %  find spikes to delete
        if  z_central_bin>z_thr

            % check if clusters are single units or multi units
            clus_kind_1 = unique(spikeInfos.unitClass(spikeInfos.channelID==ch1 & ...
                               spikeInfos.clusterID==clus_nr1));
            clus_kind_2 = unique(spikeInfos.unitClass(spikeInfos.channelID==ch2 & ...
                               spikeInfos.clusterID==clus_nr2));

            idx_clus_1=spikeInfos.channelID==ch1 & ...
                       spikeInfos.clusterID==clus_nr1;

            idx_clus_2=spikeInfos.channelID==ch2 & ...
                       spikeInfos.clusterID==clus_nr2;  
                   
                   
            % get the bundle IDs of booth cluster
            bundleID_1=unique(spikeInfos.bundleID(spikeInfos.channelID==ch1 & ...
                               spikeInfos.clusterID==clus_nr1));

            bundleID_2=unique(spikeInfos.bundleID(spikeInfos.channelID==ch2 & ...
                               spikeInfos.clusterID==clus_nr2));

            % if one of the clusters is marked as an Artifact and
            % the option delete_artifacts is set to 1, central spikes
            % within booth clusters will be deleted
            if ((strcmp(clus_kind_1,'A')||strcmp('A',clus_kind_2)) && delete_artifacts==1)||...
                (bundleID_1~=bundleID_2 && delete_artifacts==1)

                A1_delete_idx= spikeInfos.channelID==ch1 & ...
                               spikeInfos.clusterID==clus_nr1  & ...
                               spikes_in_central_bin;

                A2_delete_idx= spikeInfos.channelID==ch2 & ...
                               spikeInfos.clusterID==clus_nr2  & ...
                               spikes_in_central_bin;
                          
                % identify spikes in that are in the central bin of 
                % the two selected cluster
                A1_spike_times=spikeInfos.timeStamps(A1_delete_idx)';
                A2_spike_times=spikeInfos.timeStamps(A2_delete_idx)';
                A1_delete_idx_numbers=find(A1_delete_idx);
                A2_delete_idx_numbers=find(A2_delete_idx);
                A2_delete_idx_new=zeros(size(A2_delete_idx));

                for idx1=1:numel(A1_spike_times)
                     tem_idx_A2_delete = find(abs(A1_spike_times(idx1) - A2_spike_times) < bin_width/2);
                     
                     if numel(tem_idx_A2_delete) > 0
                          A2_delete_idx_new(A2_delete_idx_numbers(tem_idx_A2_delete))=1;
                     else
                          A1_delete_idx(A1_delete_idx_numbers(idx1)) = 0;
                     end
                end

               % save for late which spikes should be deleted
               spikes_to_delete(A1_delete_idx)=spikes_to_delete(A1_delete_idx)+1;
               spikes_to_delete(logical(A2_delete_idx_new))=spikes_to_delete(logical(A2_delete_idx_new))+1;

               % sanity check
               temp_delta_t= pdist2(spikeInfos.timeStamps(logical(A2_delete_idx_new)),...
                   spikeInfos.timeStamps(logical(A1_delete_idx)));
            
               % changed from smaller equal to smaller to avoid error if dist is exactly bin_width/2
               assert(sum(sum(temp_delta_t< bin_width/2)) == N_central_bin,...
               'Number of spikes to delete do not match the cross-correlations!');
           
               if N_central_bin < sum(A1_delete_idx) ||  N_central_bin < sum(A2_delete_idx_new) 
                   warning('Number of spikes to delete not matches cross_corr mat');
               end

               fprintf('Deleting %i artefact spikes in channel %i cluster %i! \n', ...
                    sum(A1_delete_idx),ch1,clus_nr1)
               fprintf('Deleting %i artefact spikes in channel %i cluster %i! \n', ...
                    sum(A2_delete_idx_new),ch2,clus_nr2)
                
            end

            % If a SU and a MU are compared, we delete in the same bundle the 
            % dublicated spikes in the central bin of the MU
            if ((strcmp(clus_kind_1,'MU') && strcmp(clus_kind_2,'SU')) ||...
               (strcmp(clus_kind_1,'SU') && strcmp(clus_kind_2,'MU'))) && ...
                   bundleID_1==bundleID_2
                    
                % Find spikes in the central bin
                % These are all spikes within the time bin_width
                if strcmp(clus_kind_1,'MU') && strcmp(clus_kind_2,'SU')
                    ch_MU=ch1; clus_MU=clus_nr1; ch_SU=ch2; clus_SU=clus_nr2; 
                elseif strcmp(clus_kind_2,'MU') && strcmp(clus_kind_1,'SU')
                    ch_MU=ch2; clus_MU=clus_nr2; ch_SU=ch1; clus_SU=clus_nr1; 
                end

                % loop over spikes and delete MU spikes 
                MU_delete_idx= spikeInfos.channelID==ch_MU & ...
                            spikeInfos.clusterID==clus_MU  & ...
                            spikes_in_central_bin;
                SU_keep_idx = spikeInfos.channelID==ch_SU & ...
                            spikeInfos.clusterID==clus_SU  & ...
                            spikes_in_central_bin;

                MU_spike_times=spikeInfos.timeStamps(MU_delete_idx)';
                SU_spike_times=spikeInfos.timeStamps(SU_keep_idx)';
                MU_delete_idx_numbers=find(MU_delete_idx);

                for idxMU=1:numel(MU_spike_times)
                    temp_min_time_diff=min(abs(MU_spike_times(idxMU) - ...
                                                     SU_spike_times));
                     if temp_min_time_diff >= bin_width/2
                          MU_delete_idx(MU_delete_idx_numbers(idxMU)) = 0;
                     end
                end

               if sum(MU_delete_idx) > N_central_bin
                   warning('Number of spikes to delete not matches cross_corr mat');
               end
                % save for late which spikes should be deleted
                spikes_to_delete(MU_delete_idx)=spikes_to_delete(MU_delete_idx)+1;
                fprintf('Deleting %i Spikes in channel %i cluster %i! \n', ...
                    sum(MU_delete_idx),ch_MU,clus_MU)

            % If booth clusters are MU or booth SU in the same bundle, delete spikes
            % in the cluster with a lower signal / noise 
            elseif ((strcmp(clus_kind_1,'SU') && strcmp(clus_kind_2,'SU')) || ...
                   (strcmp(clus_kind_1,'MU') && strcmp(clus_kind_2,'MU'))) && ...
                   bundleID_1==bundleID_2
                   

               % calculate the median amplitude of the clusters 
                amplitude_clus_1 = nanmedian(spikeInfos.SpikeShapes(idx_clus_1,index_ampSpi));
                amplitude_clus_2 = nanmedian(spikeInfos.SpikeShapes(idx_clus_2,index_ampSpi));

                % get the noise level of the clusters
                noise_clus_1 = spikeInfos.threshold(find(idx_clus_1,1,'first'));
                noise_clus_2 = spikeInfos.threshold(find(idx_clus_2,1,'first'));

                signal_to_noise_clus_1 = abs(amplitude_clus_1)/noise_clus_1;
                signal_to_noise_clus_2 = abs(amplitude_clus_2)/noise_clus_2;

                if signal_to_noise_clus_1 < signal_to_noise_clus_2
                    ch_Low=ch1; clus_Low=clus_nr1; ch_High=ch2; clus_High=clus_nr2; 
                elseif signal_to_noise_clus_2 < signal_to_noise_clus_1
                    ch_Low=ch2; clus_Low=clus_nr2; ch_High=ch1; clus_High=clus_nr1; 
                end                           
                
                % for comparisons within each channel
                if ch1==ch2 && clus_nr1==clus_nr2
                    ch_Low=ch1; clus_Low=clus_nr1; ch_High=ch2; clus_High=clus_nr2; 
                end

                % loop over spikes and delete MU spikes 
                Low_delete_idx=spikeInfos.channelID==ch_Low & ...
                            spikeInfos.clusterID==clus_Low  & ...
                            spikes_in_central_bin;

                High_keep_idx = spikeInfos.channelID==ch_High & ...
                            spikeInfos.clusterID==clus_High   & ...
                            spikes_in_central_bin;

                Low_spike_times  = spikeInfos.timeStamps(Low_delete_idx)';
                High_spike_times = spikeInfos.timeStamps(High_keep_idx)';
                Low_delete_idx_numbers=find(Low_delete_idx);

                if ch1==ch2 && clus_nr1==clus_nr2
                   % delete the secound spike that apears within less than
                   temp_secound_spike_in_same_chan=find(bin_width/2  > diff(Low_spike_times))+1;
                   Low_delete_idx(Low_delete_idx_numbers(setdiff(1:numel(Low_delete_idx_numbers),temp_secound_spike_in_same_chan)))=0;
                   
                   if ~(N_central_bin/2==sum(Low_delete_idx))
                       if sum(Low_delete_idx) > N_central_bin/2
                         warning('Number of spikes to delete not matches cross_corr mat');
                       end
                  end
                   
                else

                    for idxMU=1:numel(Low_spike_times)

                        temp_min_time_diff=min(abs(Low_spike_times(idxMU) - ...
                                                         High_spike_times));

                         if temp_min_time_diff >= bin_width/2
                              Low_delete_idx(Low_delete_idx_numbers(idxMU)) = 0;
                         end
                    end
                   
                   
                   if  sum(Low_delete_idx) > N_central_bin
                         warning('Number of spikes to delete not matches cross_corr mat');
                   end
                end

                % save for late which spikes should be deleted
                spikes_to_delete(Low_delete_idx)=spikes_to_delete(Low_delete_idx)+1;        
                fprintf('Deleting %i Spikes in channel %i cluster %i! \n', ...
                    sum(Low_delete_idx),ch_Low,clus_Low)
            end

            % plot and save cros-corr
            if do_plots
                fprintf('z-Score of the central bin %f (chan %i and %i) \n', ...
                             z_central_bin,ch1,ch2);
               der_plot_cross_corr_example(spikeInfos,bincounts,bin_width,hist_limit,outputpath,...
                                           session_name,ch1,ch2,clus_nr1,clus_nr2,...
                                           idx_clus_1,idx_clus_2,clus_kind_1,clus_kind_2);             
            end
        end
    end
end



%% Delete double recorded spikes
portion_of_spikes_found=100*sum(spikes_to_delete>0)/numel(spikes_to_delete);
fprintf('Time for identifying  a total of %i spikes for deletion: %.1f sec \n \n',sum(spikes_to_delete>0),toc)
fprintf('These are %.1f%% of all spikes! \n',100*sum(spikes_to_delete>0)/numel(spikes_to_delete))
fprintf('These are %.1f%% of all spikes in the central time-bin! \n \n',100*sum(spikes_to_delete>0)/sum(spikes_in_central_bin))

if sum(spikes_to_delete>0) > sum(spikes_in_central_bin)
    error('Something went wrong! To many spikes marked for deletion!')
end
% mark spikes that should be deleted due to a significant cross_corr with the prim factor 7
spikeInfos.detectionLabel(spikes_to_delete>0)=spikeInfos.detectionLabel(spikes_to_delete>0)*7;

% print some information on the labeled spikes
N_art_del7 = sum(strcmp(spikeInfos.unitClass(mod(spikeInfos.detectionLabel,7)==0),'A'));
N_MU_del7  = sum(strcmp(spikeInfos.unitClass(mod(spikeInfos.detectionLabel,7)==0),'MU'));
N_SU_del7  = sum(strcmp(spikeInfos.unitClass(mod(spikeInfos.detectionLabel,7)==0),'SU'));

N_art = sum(strcmp(spikeInfos.unitClass,'A'));
N_MU  = sum(strcmp(spikeInfos.unitClass,'MU'));
N_SU  = sum(strcmp(spikeInfos.unitClass,'SU'));

fprintf(' %.2f %% of all SU spikes (%i) were labeled. \n', 100*N_SU_del7/N_SU,N_SU)
fprintf(' %.2f %% of all MU spikes (%i) were labeled. \n', 100*N_MU_del7/N_MU,N_MU)
fprintf(' %.2f %% of all artifact spikes (%i) were labeled. \n', 100*N_art_del7/N_art,N_art)


if do_plots
    %% save output data
    save(sprintf('spikeInfos_outputcross_corr_%s.mat',session_name), 'spikeInfos');
end

end % of function

