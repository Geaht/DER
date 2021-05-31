function [inspk, inputs] = der_wavedec(spikeShapes)
%der_wavedec
%   der_wavedec creates a wavelet decomposition of spike event shapes using
%   Haar-wavelets
%
%
%   Adapted from Wave_clus 3, https://github.com/csn-le/wave_clus; see
%   Chaure et al: A Novel and fully automatic spike-sorting implementation
%   with variable numer of features, J Neurophysiol. 2018


% number of scales for the wavelet decomposition
scales=5; % Combinato and Wave_clus uses sceales = 4, we use scales = 5; pro

% number of inputs to the clustering
inputs=10;

% Wavelet decomposition
nspk=size(spikeShapes,1);
ls = size(spikeShapes,2);
cc=zeros(nspk,ls);

for spi=1:nspk                                
    [c]=wavedec(spikeShapes(spi,:),scales,'haar');
    cc(spi,1:ls)=c(1:ls);
end


for bin=1:ls                                  % KS test for coefficient selection   
    thr_dist = std(cc(:,bin)) * 3;
    thr_dist_min = mean(cc(:,bin)) - thr_dist;
    thr_dist_max = mean(cc(:,bin)) + thr_dist;
    aux = cc(find(cc(:,bin)>thr_dist_min & cc(:,bin)<thr_dist_max),bin);

    if length(aux) > 10
        [ksstat]=der_test_ks(aux);
        sd(bin)=ksstat;
    else
        sd(bin)=0;
    end
end
[~, ind]=sort(sd);
coeff(1:inputs)=ind(ls:-1:ls-inputs+1);

%CREATES INPUT MATRIX FOR SPC
inspk=zeros(nspk,inputs);
for i=1:nspk
    for j=1:inputs
        inspk(i,j)=cc(i,coeff(j));
    end
end

end