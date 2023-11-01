function [ thisResult ] = fct_calcPeriodicity( thisChan, thisEEG_BP_smooth, Fs)
%fct_calcPeriodicity 
% INPUT
%   thisChan,   table of events for which to calc periodicity
%   thisEEG_BP_smooth,    full time-series of BP power from channel where
%          events were called
%   Fs, the sampling frequency.
% OUTPUT
%   thisResult,   a table of periodicity calculations based on the events
%       provided in thisChan

% Fs = 200;
Fs = round(Fs);  % approximate
t_lo = 0.33;
t_hi = 3;

% calculate IDI based on the difference between the centroid of
% events. there are several other points within the event that could be
% used alternatively.
if (~isempty(thisChan) & ~isnan(thisChan.Area) & size(thisChan,1)>1)  % not an empty table, a nan table, and has more than 1 entry
    theseIDI = diff(thisChan.WeightedCentroid(:,1))./Fs;
    % figure;histogram(theseIDI,[0:0.25:7])
    
    % fit distribution to log2 transformed data (given distribution around mean
    % is not symmetric d/t ~ 0.5 - 2x)
    try
        thisIDI_log2 = mle(log2(theseIDI));
        thisIDI = 2.^thisIDI_log2;    % phat = mu, sigma
        thisIDI_Hz =1./thisIDI;
        thisIDI_prctNearMode = (numel(theseIDI)-sum(theseIDI>t_hi*thisIDI(1) | theseIDI<t_lo*thisIDI(1)))/numel(theseIDI);
        
    catch ME
        % if not enough data, just use the single value
        if (strcmp(ME.identifier,'stats:ProbDistUnivParam:fit:InsufficientData'))
            thisIDI_Hz = [1./theseIDI(1) NaN];  % mu, sigma
            thisIDI_prctNearMode = NaN;   
        else
        rethrow(ME);
        end
    end
    
    thisIDI_Hz_med = 1./median(theseIDI);
    
    % what % of IDIs are within 0.33-3x of the estimated mode from mle?
    % this is ~ to 'periodicity' from Ruijiter et al. (they used 0.25-1.25 of the median,
    % which I think is not appropriate, as 'IDI' values as not symmetric around the median).

    thisIDI_prctNearMed =  (numel(theseIDI)-sum(theseIDI>t_hi*median(theseIDI) | theseIDI<t_lo*median(theseIDI)))/numel(theseIDI);
    
    thisLOC = zeros(1,size(thisEEG_BP_smooth,2));
    thisLOC(cat(1,thisChan.PixelIdxList{:}))=1;
    %
    % theseLOC = fct_LOCfromEventsTbl(channels,avgEvents)
    % thisLOC = sum(theseLOC,1);
    % figure; plot(thisLOC);
    
    
    %% calculate periodicity based on cross-correlation
    [r,lags]= xcorr(thisLOC,5*Fs,'coeff');
    % figure;plot(lags,r)
    %
    % % find the first peak beyond 0 that is >0.3, report peak height, lag,
    % % prominence, and FWHM --
    % [pks,locs,w,p] = findpeaks(r(5*Fs:end),Fs);
    %
    % % find the next peak closest to zero whose prominence >0
    % indx_prom_gt_0 = round(p,3)>0;
    % indx_prom_gt_0_sub = find(indx_prom_gt_0);
    %
    % [~,thisIndx_minLag] = sort(locs(indx_prom_gt_0));  % includes 0-lag peak
    %
    % thisIndx = indx_prom_gt_0_sub(thisIndx_minLag(2));
    
    % find the first peak beyond 0 that is >0.3, report peak height, lag,
    % prominence, and FWHM --
    [pks,locs,w,p] = findpeaks(r(5*Fs:end),Fs,'MinPeakDistance',0.5);
    
    % find the next peak closest to zero whose prominence >0.1, pk >0.2
    indx_prom_gt_0 = [round(p,3)>0.1 & round(pks,3)>0.2];
    indx_prom_gt_0_sub = find(indx_prom_gt_0);
    
    [~,thisIndx_minLag] = sort(locs(indx_prom_gt_0));  % includes 0-lag peak
    
    % this code should be re-written as a multi-objective optimization to
    % choose peak that minimizes lag distances while max (pk, prominence)
    
    if (any(indx_prom_gt_0))
    thisIndx = indx_prom_gt_0_sub(thisIndx_minLag(1));  % for some reason, no longer includes zero-lag
    
    
    thisIDI_xcorr_pk = pks(thisIndx);
    thisIDI_xcorr_loc = locs(thisIndx);  % in seconds
    thisIDI_xcorrs_w = w(thisIndx);  % width at half prominence
    thisIDI_xcorr_prm = p(thisIndx); % prominence
    thisIDI_xcorr_Hz = 1/locs(thisIndx);
    else
        
    thisIDI_xcorr_pk = nan;
    thisIDI_xcorr_loc = nan;  % in seconds
    thisIDI_xcorrs_w = nan;  % width at half prominence
    thisIDI_xcorr_prm = nan; % prominence
    thisIDI_xcorr_Hz = nan;
        
    end
    
    %% compute cross-correlation for the smoothed BP data (should be less biased
    % by segmentation error)
    % calculate periodicity based on cross-correlation
    [r,lags]= xcorr(thisEEG_BP_smooth,5*Fs,'coeff');
    % figure;plot(lags,r)
    
    % find the first peak beyond 0 that is >0.3, report peak height, lag,
    % prominence, and FWHM --
    [pks,locs,w,p] = findpeaks(r(5*Fs:end),Fs,'MinPeakDistance',0.5);
    
%     figure; plot(lags(5*Fs:end)./Fs,simplema(r(5*Fs:end),10,2)); xlim([0 5]);
    
    % find the next peak closest to zero whose prominence >0.1, pk >0.2
    indx_prom_gt_0 = [round(p,3)>0.1 & round(pks,3)>0.2];
    indx_prom_gt_0_sub = find(indx_prom_gt_0);
    
    [~,thisIndx_minLag] = sort(locs(indx_prom_gt_0));  % includes 0-lag peak
    
    % this code should be re-written as a multi-objective optimization to
    % choose peak that minimizes lag distances while max (pk, prominence)
    if (any(indx_prom_gt_0))
        thisIndx = indx_prom_gt_0_sub(thisIndx_minLag(1));  % for some reason, no longer includes zero-lag
        
        thisIDI_xcorr_BP_pk = pks(thisIndx);
        thisIDI_xcorr_BP_loc = locs(thisIndx);  % in seconds
        thisIDI_xcorr_BP_w = w(thisIndx);  % width at half prominence
        thisIDI_xcorr_BP_prm = p(thisIndx); % prominence
        thisIDI_xcorr_BP_Hz = 1/locs(thisIndx);
    else
        thisIDI_xcorr_BP_pk = nan;
        thisIDI_xcorr_BP_loc = nan;  % in seconds
        thisIDI_xcorr_BP_w = nan;  % width at half prominence
        thisIDI_xcorr_BP_prm =nan; % prominence
        thisIDI_xcorr_BP_Hz = nan;
        
    end
    
    
    %% Calculate power of discharges relative to total power for channel.
    
    power_discharges = sum(cat(2,thisChan.PixelValues{:}));   % will be BP power from events
    power_total = sum(thisEEG_BP_smooth,2);
    power_discharges_rel = power_discharges/power_total;
    

    %% Calculate evolution. 

    indx_IDI = [1:length(theseIDI)]';

    mdl=fitlm(indx_IDI,theseIDI.^-1);
%     mdl=fitlm(indx_IDI,theseIDI.^-1,'y~x1-1','RobustOpts','on');
    evolution_slope = mdl.Coefficients.Estimate(2);   % units is delta(Hz)/delta(step)
    evolution_SE = mdl.Coefficients.SE(2);  % standard error
    evolution_RMSE = mdl.RMSE;

    %% Make a results table
    thisResult = table;
    thisResult.channel = thisChan.channel(1);
    thisResult.channelNum = unique(thisChan.channelNum);
    thisResult.events = size(thisChan,1);
    thisResult.eventRate = size(thisChan,1)/(size(thisEEG_BP_smooth,2)/Fs);
    thisResult.thisIDI_Hz = thisIDI_Hz;
    thisResult.thisIDI_Hz_med = thisIDI_Hz_med;
    thisResult.thisIDI_prctNearMed = thisIDI_prctNearMed;
    thisResult.thisIDI_prctNearMode = thisIDI_prctNearMode;
    thisResult.thisIDI_xcorr_Hz = thisIDI_xcorr_Hz;
    thisResult.thisIDI_xcorr_lag = thisIDI_xcorr_loc;
    thisResult.thisIDI_xcorr_prm = thisIDI_xcorr_prm;
    thisResult.thisIDI_xcorr_BP_Hz = thisIDI_xcorr_BP_Hz;
    thisResult.thisIDI_xcorr_BP_lag = thisIDI_xcorr_BP_loc;
    thisResult.thisIDI_xcorr_BP_prm = thisIDI_xcorr_BP_prm;
    thisResult.power_discharges = power_discharges;
    thisResult.power_total = power_total;
    thisResult.power_discharges_rel = power_discharges_rel;
    thisResult.evolution_slope = evolution_slope;
    thisResult.evolution_SE = evolution_SE;
    thisResult.evolution_RMSE = evolution_RMSE;

else
    
    %oops!
       
    thisResult = table;
    thisResult.channel = thisChan.channel(1);
    thisResult.channelNum = unique(thisChan.channelNum);
    
    if (~isnan(thisChan.Area) & size(thisChan,1)==1) 
    thisResult.events = size(thisChan,1);
    thisResult.eventRate = 1/(size(thisEEG_BP_smooth,2)/Fs);
    else
    thisResult.events = 0;
    thisResult.eventRate = 0;
    end
        
    thisResult.thisIDI_Hz = [nan nan];
    thisResult.thisIDI_Hz_med = nan;
    thisResult.thisIDI_prctNearMed = nan;
    thisResult.thisIDI_prctNearMode = nan;
    thisResult.thisIDI_xcorr_Hz = nan;
    thisResult.thisIDI_xcorr_lag = nan;
    thisResult.thisIDI_xcorr_prm = nan;
    thisResult.thisIDI_xcorr_BP_Hz = nan;
    thisResult.thisIDI_xcorr_BP_lag = nan;
    thisResult.thisIDI_xcorr_BP_prm = nan;
    thisResult.power_discharges = nan;
    thisResult.power_total = nan;
    thisResult.power_discharges_rel = nan;
    thisResult.evolution_slope = nan;
    thisResult.evolution_SE = nan;
    thisResult.evolution_RMSE = nan;


end
end
