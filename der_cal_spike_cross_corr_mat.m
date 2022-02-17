function [cross_corr_mat_info] = ...
            der_cal_spike_cross_corr_mat(spikeInfos,hist_limit,bin_width,outputpath,do_plots,session_name,skip_labeled_spikes)
%% output_spikes=dds_cal_spikes_corr_mat(spikeInfos,hist_limit,bin_width,outputpath,do_plots,save_data_to_file)
% Function to calculate matrix with the cross-correlations of all clusters in a session
%
% Inputs:
%
% spikeInfos: Table containing spike information
%             (output of the function get_spikeInfos)
%
% hist_limit: Maximal lag of spikes counted in the cross_corr matrix 
%             (Default 20 ms)
%
% bin_width: Width of time bins within the cross_corr matrix
%            (Default 0.5 ms)
%
% outputpath: Path where data and output plots shall be stored.
%             (default current directory)
%
% do_plots: If do_plots==1 plots matrix of z-values of the central time bin
%           of all possible cluster combinations
%
% skip_labeled_spikes: If skip_labeled_spikes==1, all skpikes that
%                         were labeled by the dds algorythms are skipped
%                         
%           
% Outouts:
%
% cross_corr_mat: Matrix (N_cluter X Ncluster X 2 time-bins) 
%          stores cross-correlation counts of all cluster combinations at different time
%          lags
%
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.

%% Input parameters and default values
if ~exist('hist_limit','var')
    hist_limit=20;  % ms
end

if ~exist('bin_width','var')
    bin_width=0.5; % ms
end

if ~exist('outputpath','var')
    outputpath=pwd; 
end

if ~exist('do_plots','var')
    do_plots=1; 
end

if ~exist('session_name','var')
    session_name='temp'; 
end

if ~exist('skip_labeled_spikes','var')
    skip_labeled_spikes=0; 
end


if  mod(hist_limit,bin_width)~=0
    error('The maximal time-lag (hist_limit) must be a multiple of the bin width!')
end

% If there are no input data are provided, load spikes data from current directory
if  ~exist('spikeInfos','var')
    try    
        fprintf('Loading spikes from current folder! \n')
        spike_files=dir('spikeInfos_*');
        if numel(spike_files)==1
            load(spike_files(1).name);
        else 
            error('More than one spikeInfo file in this folder, please specify!');
        end
    catch
        error('No spike data profided or found!');
    end
end

% Check if spikes were already sorted acording to their time stemps 
if ~issorted(spikeInfos.timeStamps)
    [~,idxsort] = sort(spikeInfos.timeStamps);
    spikeInfos=spikeInfos(idxsort,:);
    warning('Input not sorted acording to time stamps!');
    fprintf('Sorting input spikes... \n \n')
end


% exclude all spikes marked by the DDS algo
if skip_labeled_spikes==1
    spikeInfos=spikeInfos(spikeInfos.detectionLabel==1,:);
end

%% assure DER tools are in the path 
code_path=which('der_cal_spike_cross_corr_mat');
addpath(genpath(code_path(1:end-30)));

%% Calculate all cross-correlations in the session

% get list of all cluster in the session
cluster_list=unique([spikeInfos.channelID spikeInfos.clusterID],'rows');
n_cluster=numel(cluster_list(:,1));
central_bin_idx=floor(numel(-hist_limit:bin_width:hist_limit)/2)+1;

% generate the cross-correlation matrix for all cluster in the session
cross_corr_mat=zeros(n_cluster,n_cluster,(2*hist_limit/bin_width)+1);

% save for eac spike it it is in a central time-bin
spikes_in_central_bin=zeros(numel(spikeInfos.timeStamps),1);
percentage_job_done_counter=1;

fprintf('Calculating cross-correlations for all cluster combinations... \n \n')
tic % time calculation of cross correlation

% loop over all spikes in spikeInfos
for idx_lead_spike=1:numel(spikeInfos.timeStamps)
    
    % monitor progress
    if idx_lead_spike/numel(spikeInfos.timeStamps) > 0.01*percentage_job_done_counter
        % calculate the expected computation time till job is done based on first 1% of data
        if percentage_job_done_counter==1
            fprintf('Expected time till all cross_corrs are calculated: %.2f min  \n',100*toc/60);
            fprintf('__________________________________________________________________ \n     \n');
        end
        
        % plot a progressbar in output
        der_progressbar('Calculating cross_corrs',percentage_job_done_counter);
        percentage_job_done_counter=percentage_job_done_counter+1;
    end
    
    % go to next spike after the lead spike
    temp_timeStamp=spikeInfos.timeStamps(idx_lead_spike);
    if idx_lead_spike<numel(spikeInfos.timeStamps) 
        idx_temp_spikes=idx_lead_spike+1;
    else
        continue % finish when the last spike is reached
    end
    
    % get cluster information for the lead spike
    cluster_info_lead_spike   = [spikeInfos.channelID(idx_lead_spike)...
                                 spikeInfos.clusterID(idx_lead_spike)];
    cluster_number_lead_spike = ismember(cluster_list,cluster_info_lead_spike,'rows');

    % check if next spike is within the range of the histogram
    % (e.g. hist_limit+bin_width/2.0)
    while ((spikeInfos.timeStamps(idx_temp_spikes)-temp_timeStamp) < hist_limit+bin_width/2.0) && ...
           idx_temp_spikes < numel(spikeInfos.timeStamps)
        
        % calulate time difference to next spike
        diff_spike_times=spikeInfos.timeStamps(idx_temp_spikes)-temp_timeStamp;
    
        % if the next spike is within less than half a bin_width 
        % mark booth spikes for later as central bin spikes
        if diff_spike_times < bin_width/2.
            spikes_in_central_bin(idx_temp_spikes)=1;
            spikes_in_central_bin(idx_lead_spike)=1;
        end
        
        % information of secound spike
        cluster_info_temp_spike   = [spikeInfos.channelID(idx_temp_spikes)...
                                     spikeInfos.clusterID(idx_temp_spikes)];
        cluster_number_temp_spike  = ismember(cluster_list,cluster_info_temp_spike,'rows');

        % round diff_spike_times to use it as index and increase count in cross-correlation mat
        cross_corr_mat(cluster_number_lead_spike,cluster_number_temp_spike,...
                         central_bin_idx+round(diff_spike_times/bin_width))...
            =   1 +  cross_corr_mat(cluster_number_lead_spike,cluster_number_temp_spike,...
                         central_bin_idx+round(diff_spike_times/bin_width));
        
        % go to next spike
        if idx_temp_spikes<numel(spikeInfos.timeStamps)
            idx_temp_spikes=idx_temp_spikes+1;
        else
            continue % finish if last spike is reached
        end
    end % of while 
end % of for all spikes

% print the time needed to calculate all cross-correlations
corr_calculation_time=toc;
fprintf('\n Time needed to calculate all cross-correlations: %i min \n',corr_calculation_time/60.0);

% transform trianguar matrix to full matrix
A = cross_corr_mat; B = flip(A,3);  C = permute(B ,[2 1 3]); D=A+C;
D(1:end,1:end,central_bin_idx) = A(:,:,central_bin_idx) + C(:,:,central_bin_idx);
cross_corr_mat=D;

cross_corr_bin_width=bin_width;

%% save data
save([outputpath filesep sprintf('cross_corr_mat_%i_ms_width_%i_ms_bins_%s.mat',hist_limit,bin_width,session_name)],'cross_corr_mat','spikes_in_central_bin','cross_corr_bin_width');

%% plot the normalized correlation matrix
    
% z-normalize the cross_corr matrix
norm_cross_corr_mat=zeros(size(cross_corr_mat));
central_z_mat=zeros(size(cross_corr_mat,1),size(cross_corr_mat,2));

for iclus1=1:size(cross_corr_mat,1)
    for iclus2=1:size(cross_corr_mat,2)
        norm_cross_corr_mat(iclus1,iclus2,:) = cross_corr_mat(iclus1,iclus2,:)/sum(cross_corr_mat(iclus1,iclus2,:));
        temp_bincounts= squeeze(cross_corr_mat(iclus1,iclus2,:))';
        
        if std(temp_bincounts(1:end~=central_bin_idx))  == 0 
            central_z_mat(iclus1,iclus2) = 0;
            
        else
            central_z_mat(iclus1,iclus2) = ...
                (temp_bincounts(central_bin_idx)-mean(temp_bincounts(1:end~=central_bin_idx)))...
                    /std(temp_bincounts(1:end~=central_bin_idx));
        end
        
        if isnan(central_z_mat(iclus1,iclus2))
           fprintf('Nan, something went wrong!')        
        end
           
    end
end
    

% save some information for later plotting
clear cluster_labels hemisphere_labels
label_idx=zeros(size(cluster_list,1),1);
label_idx(1)=1;
for ii=1:size(cluster_list,1)
    idx_temp_clus=spikeInfos.channelID == cluster_list(ii,1) & spikeInfos.clusterID == cluster_list(ii,2);

    cluster_labels{ii} = cell2mat(unique(spikeInfos.region(idx_temp_clus)));

    if ii>1 && strcmp(cluster_labels{ii},tmep_clusterlabel)
        cluster_labels{ii} = '';
    else
        tmep_clusterlabel= cell2mat(unique(spikeInfos.region(idx_temp_clus)));
        label_idx(ii)=1;
    end

    %hemisphere_labels{ii} = cluster_labels{ii}(1);
end


if do_plots==1
    % plot z-value of central bins
    ff=figure();
    imagesc(central_z_mat(:,:))
    colormap(jet(100));
    c=colorbar;
    caxis([-5 30]);
    
    ax=gca;
    ax.XTick=1:size(central_z_mat(:,:),1);
    ax.XTickLabel=cluster_labels;
    ax.XTickLabelRotation=90;
    ax.YTick=1:size(central_z_mat(:,:),1);
    ax.YTickLabel=cluster_labels;
    ax.XTickLabelRotation=90;
    ax.YTickLabelRotation=0;
    c.Label.String='z';

    pbaspect([1 1 1])
    ff.Position= [ 0 0 900 800];
    set(findall(ax,'-property','FontSize'),'FontSize',24);
    set(findall(c,'-property','FontSize'),'FontSize',24);

    title('z-values of central bin','FontSize',26)

    % file names 
    file_name=[outputpath filesep session_name '-' ...
        sprintf('z-value-cross-corr-central-bin_with_artifacts_%i_ms_width_%i_mus_bins.png',...
        hist_limit,1000*bin_width)];
    % save figure in outputpath
    print(ff,file_name,'-dpng', '-r300');
end


% save relevant information in an output struct 
cross_corr_mat_info.cross_corr_mat = cross_corr_mat;
cross_corr_mat_info.hist_limit = hist_limit;
cross_corr_mat_info.bin_width = bin_width;
cross_corr_mat_info.spikes_in_central_bin = spikes_in_central_bin;
cross_corr_mat_info.corr_calculation_time = corr_calculation_time;
cross_corr_mat_info.norm_cross_corr_mat = norm_cross_corr_mat;
cross_corr_mat_info.central_z_mat = central_z_mat;
cross_corr_mat_info.cluster_labels=cluster_labels;
cross_corr_mat_info.session_name=session_name;
cross_corr_mat_info.cluster_list=cluster_list;
cross_corr_mat_info.n_cluster=n_cluster;
cross_corr_mat_info.central_bin_idx=central_bin_idx;

save([outputpath filesep sprintf('cross_corr_mat_%i_ms_width_%i_ms_bins_%s.mat',hist_limit,bin_width,session_name)],'cross_corr_mat','spikes_in_central_bin','cross_corr_bin_width','cross_corr_mat_info');


end % of function
