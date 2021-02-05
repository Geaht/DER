function []=der_progressbar(input_string,progress,num_of_steps_total,num_of_steps_to_show)
%der_progressbar
%   der_progressbarplots an inline text progressbar for calculations
%   []=dds_progressbar(input_string,progress,num_of_steps_total,num_of_steps_to_show)
%   input_string: text shown in front of the progress bar
%   progress: percent of progress done or if num_of_steps_total is profided
%   number of steps done
%   num_of_steps_to_show: number of symbols in the progress bar
%
%
%   Licence:
%   This source code form is subject to the terms of the Mozilla Public
%   Licence, v. 2.0. if a copy of the MPL was not distributed with this file,
%   you can optain one at http://mozilla.org/MPL/2.0/.

    if ~exist('num_of_steps','var')    
        num_of_steps_to_show=20;
    end
    
    if ~exist('num_of_steps_total','var')    
        progress_percent=progress;
    else
        progress_percent=floor(100*progress/num_of_steps_total);
    end
    
    % character to use in the bar default=9640 is a shaded box
    char_to_use=char(9640);
    
    bar_string_progress=repmat(char_to_use,1,1+floor(progress_percent*num_of_steps_to_show/100));
    bar_string_rest=repmat(' ',1,20-floor(progress_percent*num_of_steps_to_show/100)-1);
    bar_string= ['[ ' bar_string_progress  bar_string_rest ' ]'];
    output_str=sprintf('%s: %03i %%%% %s ',input_string,floor(progress_percent),bar_string);
    
    if floor(progress_percent)==1
        fprintf(output_str)
    else
        refreshbar=repmat('\b',1,numel(output_str)-1);
        fprintf([refreshbar output_str]);
    end
end

