function [ eeg_nleo ] = calc_nleo( seg, Fs)
%calc_nleo Compute non-linear energy operator. 
% INPUT
% ***
% OUTPUT
% ***

seg_p = [zeros(size(seg,1),3) seg zeros(size(seg,1),3)];

eeg_nleo = arrayfun(@(x,x_1,x_2,x_3) (x_1*x_2)-(x*x_3), seg_p(:,4:end-3), seg_p(:,3:end-4),...
    seg_p(:,2:end-5), seg_p(:,1:end-6));

% eeg_nleo_s = simplema(eeg_nleo,(1000/Fs)*120,2);


end

