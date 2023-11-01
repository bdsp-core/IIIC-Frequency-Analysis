function [ LOC, theseEvents, eeg_bp_smooth ] = fct_findEvents_hilo_dynamic( seg, Fs, ...
    thresh_hi, thresh_lo, smoothWin, freq_lo, freq_hi, boostWin, channelLabels)
%fct_findEvents_hilo_dynamic.  Use a high-low threshold approach to find events. Naive to specific BP or montage.
% Accepts thresh_hi, thresh_lo arguments that reference rms.
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
%         thresh_hi = 15*5;
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
        end
        
        
        % apply a transform (intended for bipolar montages) to boost detection of
        % events that co-occur with other events in other channels.
        eeg_bp_smooth_transf = eeg_bp_smooth.*repmat(smooth(sum(eeg_bp_smooth,1),boostWin)',18,1);
        eeg_bp_smooth = eeg_bp_smooth_transf/50;  % arbitrary scaling.
        
        for j=1:size(eeg_bp_smooth_transf,1)
           
           thisEEG_bp_smooth = eeg_bp_smooth(j,:); 
            
           % high-low threshold
           bw_hi = thisEEG_bp_smooth>thresh_hi;
           bw_lo = thisEEG_bp_smooth>thresh_lo;
           
           lpds_hi = regionprops(bw_hi,'PixelIdxList');
%            lpds_lo = regionprops('table',thisEEG_bp_smooth>thresh_lo,'PixelIdxList');
           
           indx_hi = cat(1, lpds_hi.PixelIdxList);
%            indx_lo = cat(1, lpds_lo.PixelIdxList);
           
           isect = bwselect([bw_lo; bw_lo]', ones(length(indx_hi),1), indx_hi,4);   % must be 2D, order is critical, takes binary image but needs pixel ids for selection
           isect = isect(:,1)';
           
           lpds_ = regionprops('table',isect,thisEEG_bp_smooth,'Area',...
               'BoundingBox','PixelIdxList','MaxIntensity','MeanIntensity','PixelValues','WeightedCentroid');
           
           lpds_raw = regionprops('table',isect,eeg_bp(j,:),...
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
        
        LOC = fct_LOCfromEventsTbl(channels, theseEvents); 
        

end

