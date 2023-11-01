function [ LOC, theseEvents] = fct_findEvents_simple( data_log, thisEEG_bp_smooth, channelLabels)
%fct_findEvents_simple.  Simple version when the input is already logical,
% just want to find the events

% INPUT data_log, 1xm double/logical vector 
%      ***
% OUTPUT 
%       LOC, nxm logical array corresponding to positions of candidate
%           discharges.
%       theseEvents, nx1 cell array of structs containing stat objects
%       from region props
%
        
        channels = channelLabels;
        LOC = zeros(size(data_log));
        totalSamples = size(data_log,2);
        
%         thresh_hi = 15*5;
%         smoothWin = 20;
%                 
        % band-pass filter
%         eeg_bp = eegfilt(seg,Fs,freq_lo,freq_hi);   % from EEG Toolbox
%                 
%         eeg_bp_smooth = nan(size(eeg_bp));
        
%         theseEvents = cell(size(eeg_bp,1),1);
        theseEvents = table();
        
%         for j=1:size(eeg_bp,1)
%            thisEEG_bp_smooth = smooth((eeg_bp(j,:).^2)',smoothWin)';
%            eeg_bp_smooth(j,:) = thisEEG_bp_smooth;
%         end
        
%          eeg_bp_smooth_transf = eeg_bp_smooth;
%         % apply a transform (intended for bipolar montages) to boost detection of
%         % events that co-occur with other events in other channels.      
%         if (length(channels))==18
%         channels_L = [1:4,9:12];
%         channels_R = [5:8, 13:16];
% 
%         eeg_bp_smooth_transfL = eeg_bp_smooth(channels_L,:).*repmat(smooth(sum(eeg_bp_smooth,1),boostWin)',length(channels_L),1);
%         eeg_bp_smooth_transfR = eeg_bp_smooth(channels_R,:).*repmat(smooth(sum(eeg_bp_smooth,1),boostWin)',length(channels_R),1);
%         
%         eeg_bp_smooth_transf(channels_L,:) = eeg_bp_smooth_transfL/50;
%         eeg_bp_smooth_transf(channels_R,:) = eeg_bp_smooth_transfR/50;
%         
% %         eeg_bp_smooth = eeg_bp_smooth_transf/50; % arbitrary scaling    
%         
%         end
%         
        
        for j=1:size(data_log,1)
           
%            thisEEG_bp_smooth = eeg_bp_smooth(j,:); 
%            
%            eeg_prctl_hi = prctile(thisEEG_bp_smooth, thresh_hi);
%            eeg_prctl_lo = prctile(thisEEG_bp_smooth, thresh_lo);
%            
%            % high-low threshold
%            bw_hi = thisEEG_bp_smooth>eeg_prctl_hi;
%            bw_lo1 = diff(eeg_bp(j,:))>0;
%            bw_lo2 = diff(eeg_bp(j,:))<0;
% %            bw_lo = thisEEG_bp_smooth>eeg_prctl_lo;
%            
%            lpds_hi = regionprops(bw_hi,'PixelIdxList');
% %            lpds_lo = regionprops('table',thisEEG_bp_smooth>thresh_lo,'PixelIdxList');
%            
%            indx_hi = cat(1, lpds_hi.PixelIdxList);
% %            indx_lo = cat(1, lpds_lo.PixelIdxList);
%                       
%            isect1 = bwselect([bw_lo1; bw_lo1]', ones(length(indx_hi),1), indx_hi,4);   % must be 2D, order is critical, takes binary image but needs pixel ids for selection
%            isect1 = isect1(:,1)';
%            
%            isect2 = bwselect([bw_lo2; bw_lo2]', ones(length(indx_hi),1), indx_hi,4);   % must be 2D, order is critical, takes binary image but needs pixel ids for selection
%            isect2 = isect2(:,1)';
%            
%            isect = isect1 | isect2;
%            isect = logical([isect 0]);

            isect = data_log;
           
           lpds_ = regionprops('table',isect,thisEEG_bp_smooth,'Area',...
               'BoundingBox','PixelIdxList','MaxIntensity','MeanIntensity','PixelValues','WeightedCentroid');
%            
%            lpds_raw = regionprops('table',isect,seg,...
%                'MaxIntensity','MeanIntensity','MinIntensity','PixelValues');
%       
%            theseVarNames = lpds_raw.Properties.VariableNames;
%            
%            % append 'raw' to variable names from raw eeg (vs. values from
%            % the transformed data)
%            for k=1:length(theseVarNames)
%               thisVarName = theseVarNames{k};
%               theseVarNames(k) = {[thisVarName '_raw']};           
%            end
%            
%            lpds_raw.Properties.VariableNames = theseVarNames;
%            thisResult = cat(2,lpds_,lpds_raw);
           thisResult = lpds_;
           
           if ~isempty(thisResult)
               thisResult.channel = repmat(channels(j),size(lpds_,1),1);
               thisResult.channelNum = repmat(j,size(lpds_,1),1);
               thisResult.channelEventNum = num2cell(1:size(lpds_,1))';
    %            theseEvents = cat(1,theseEvents, thisResult);
               
               theseEvents = addToEventsList(thisResult, theseEvents);
           else
               % empty!
                disp('test');
               
                thisResult = makeNanTable(thisResult);
                thisResult.channel = channels(j);
                thisResult.channelNum = j;
                thisResult.channelEventNum = nan; 
                
                theseEvents = addToEventsList(thisResult, theseEvents);
           end
        end
        
        theseEvents.eventNum = num2cell(1:size(theseEvents,1))';
        
        LOC = fct_LOCfromEventsTbl(channels, theseEvents, totalSamples); 
        

end

