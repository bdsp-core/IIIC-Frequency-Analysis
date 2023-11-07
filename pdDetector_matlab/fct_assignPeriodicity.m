function [ theseEvents_new ] = fct_assignPeriodicity( channelLabels_bipolar,theseEvents)
%fct_assignPeriodicity Version 1
%Calculates a periodicity for each channel 

% channels = unique(theseEvents.channelNum);

Fs=200;
LOC = fct_LOCfromEventsTbl(channelLabels_bipolar,theseEvents);
numChan = size(LOC,1);

[autocor,lags] = xcorr(LOC',3*Fs,'coeff');

% for i=1:size(autocor,2)
%     [pks, locs, w, p] = findpeaks(autocor(:,i),'MinPeakDistance',0.5*Fs);
%     pks_c{i} = pks;
%     locs_c{i} = locs;
%     w_c{i} = w;
%     p_c{i} = p;
%     
% end
% 
% p_c_sq = reshape(p_c,18,18);
% pks_c_sq = reshape(pks_c, 18,18);
% locs_c_sq = reshape(locs_c,18,18);
% 
% lags_sq = cellfun(@(x,y) (x(y)./Fs)', repmat({lags},18,18), locs_c_sq,'UniformOutput',0);
% 
% [lags(locs)/Fs_n;p']'

theseCol = 1:numChan:size(autocor,2);
thesePeriods = cell(numChan,1);

theseResults = table();

for i=1:length(theseCol)
    [pks, locs, w, p] = findpeaks(autocor(:,theseCol(i)),'MinPeakDistance',0.5*Fs);
    pks_c{i} = pks;
    locs_c{i} = locs;
    w_c{i} = w;
    p_c{i} = p;
    
    if ~isempty(p)
        offset = ceil(length(p)/2);
        [mostProm, m_l] = max(p(offset+1:end)); % find the most prominent peak after 0-lag
        m_l_true = m_l +offset;
    else
        mostProm = nan;
        m_l_true = nan;
    end
    
    if ~isnan(m_l_true)
        thisPeriod = lags(locs(m_l_true))/Fs;
    else
        thisPeriod = nan;
    end
    
    thesePeriods{i} = thisPeriod;
    thisResult = table();
    
    thisResult.channelNum = i; 
    thisResult.Period1_lag = thisPeriod;
    thisResult.Period1_prom = mostProm;
    
    theseResults = cat(1,theseResults,thisResult);
    
end

% lags_c = cellfun(@(x,y) (x(y)./Fs)', repmat({lags},1,numChan), locs_c,'UniformOutput',0);



end

