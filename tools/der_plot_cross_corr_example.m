function der_plot_cross_corr_example(spikeInfos,bincounts,bin_width,hist_limit,outputpath,...
    session_name,ch1,ch2,clus_nr1,clus_nr2,idx_clus_1,idx_clus_2,clus_kind_1,clus_kind_2)
%der_plot_cross_corr_example
%   der_plot_cross_corr_example is a funtion to plot examples of cross 
%   correlogramms in der_detect_cross_corr_spikes
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.

fig =figure('visible','off');
hold on;
bb1=bar(-hist_limit-0.5*bin_width:bin_width:hist_limit-0.5*bin_width,bincounts,'histc');
bb1.FaceColor=[0.4 0.4 0.4];
ax1=gca;
xlim([-hist_limit-0.5 hist_limit+0.5]);
xlabel(' t [ms]');
ylabel('N_{spikes}');

% load spikes from cluster 1 & 2 idx_clus_1
a1=axes('Position',[0.15 .6 .25 .25]);
b1=plot(1:64,mean(spikeInfos.SpikeShapes(idx_clus_1,:)),'LineWidth',3,'Color',[0.8 0.2 0.2]);
title(sprintf('%s: %i spikes',clus_kind_1{:},sum(idx_clus_1)));
a1.Visible='off';
a2=axes('Position',[.65 .6 .25 .25]);
b2=plot(1:64,mean(spikeInfos.SpikeShapes(idx_clus_2,:)),'LineWidth',3,'Color',[0.2 0.2 0.8]);
title(sprintf('%s: %i spikes',clus_kind_2{:},sum(idx_clus_2)));
a2.Visible='off';

yscale=[min([a1.YLim,a2.YLim]), max([a1.YLim,a2.YLim])];
a1.YLim = yscale;
a2.YLim = yscale;

set(findall(fig,'-property','FontSize'),'FontSize',20);
%set(findall(fig,'-property','LineWidth'),'LineWidth',2);  
fig.Position = [0 0 600 600];
set(gca, 'FontName', 'Helvetica');


% file names 
file_name=[outputpath filesep session_name '-Cross-corr' ...
    num2str(ch1) '-clus' num2str(clus_nr1) ...
    '-ch' num2str(ch2) '-clus' ...
    num2str(clus_nr2)  '-new.png'];

% save figure in plot_dir
print(fig,file_name,'-dpng', '-r100');
hold off;
close(fig);


end

