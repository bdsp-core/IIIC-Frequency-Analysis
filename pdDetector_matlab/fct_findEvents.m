function [ LOC, theseEvents, eeg_bp_smooth ] = fct_findEvents( seg, Fs, thresh_mag, smoothWin, freq_lo, freq_hi,channelLabels)
%fct_findEvents Naive to specific BP or montage.
% ***
% INPUT seg, nxm double array of EEG in bipolar montage
%       Fs, sampling freq in sec-1
%      ***
% OUTPUT 
%       LOC, nxm logical array corresponding to positions of candidate
%           discharges.
%       theseEvents, nx1 cell array of structs containing stat objects
%       from region props
%
 
        channels = channelLabels;
        LOC = zeros(size(seg));
%         thresh_mag = 15*5;
%         smoothWin = 20;
%                 
        % band-pass filter
        eeg_bp = eegfilt(seg,Fs,freq_lo,freq_hi);   % from EEG Toolbox
                
        eeg_bp_smooth = nan(size(eeg_bp));
        
%         theseEvents = cell(size(eeg_bp,1),1);
        theseEvents = table();
        
        for j=1:size(eeg_bp,1)
           thisEEG_bp_smooth = smooth((eeg_bp(j,:).^2)',smoothWin)';
           eeg_bp_smooth(j,:) = thisEEG_bp_smooth;
           lpds_ = regionprops('table',thisEEG_bp_smooth>thresh_mag,thisEEG_bp_smooth,'Area',...
               'BoundingBox','PixelIdxList','MaxIntensity','MeanIntensity','PixelValues','WeightedCentroid');
           
           lpds_raw = regionprops('table',thisEEG_bp_smooth>thresh_mag,eeg_bp(j,:),...
               'MaxIntensity','MeanIntensity','MinIntensity','PixelValues');
      
           theseVarNames = lpds_raw.Properties.VariableNames;
           
           % append 'raw' to variable names from raw eeg (vs. values from
           % the transformed data)
           for k=1:length(theseVarNames)
              thisVarName = theseVarNames{k};
              theseVarNames(k) = {[thisVarName '_raw']};           
           end
           
           lpds_raw.Properties.VariableNames = theseVarNames;
           
           thisResult = cat(2,lpds_,lpds_raw);
           thisResult.channel = repmat(channels(j),size(lpds_,1),1);
           thisResult.channelNum = repmat(j,size(lpds_,1),1);
           thisResult.channelEventNum = num2cell(1:size(lpds_,1))';
%            theseEvents = cat(1,theseEvents, thisResult);
           
           theseEvents = addToEventsList(thisResult, theseEvents);
        end
        
        theseEvents.eventNum = num2cell(1:size(theseEvents,1))';
        
        LOC = fct_LOCfromEventsTbl(channels, theseEvents); % just include LOC for size
        

end

