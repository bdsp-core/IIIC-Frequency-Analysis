function [ LOC,artifactEvents ] = fct_findArtifacts ( seg, Fs, thresh_crazy, smoothWin, channelLabels_average)
%fct_findArtifacts Summary of this function goes here
%   Detailed explanation goes here
        
        thresh_high = 30; %demeaned - ma
        thresh_low = 10; %rms

        seg_avg = seg-mean(seg,1);
        channels = channelLabels_average;
        LOC = zeros(size(seg));
        
        artifactEvents = table();
        
        % find events with large amplitude > thresh_crazy
%         thresh_crazy= 3;
        thisRMS = rms(seg_avg);
        
        seg_ma = simplema(seg_avg,smoothWin,2);
%         fct_showEEG_avg(seg_ma);
        
        for j=1:size(seg_avg,1)
            thisEEG = seg_avg(j,:);
            thisEEG_ma = seg_ma(j,:);
%             thisEEG_bp_smooth = smooth((eeg_bp(j,:).^2)',smoothWin)';
%             eeg_bp_smooth(j,:) = thisEEG_bp_smooth;
% 
%             bad_hiamp = regionprops('table',...
%                 thisEEG>(thresh_crazy*thisRMS)| thisEEG<-1*(thresh_crazy*thisRMS),...
%                 'Area','BoundingBox','PixelIdxList');
            
%             bad_hiamp = regionprops('table',...
%                 thisEEG>thresh_crazy | thisEEG<-1*thresh_crazy,...
%                 'Area','BoundingBox','PixelIdxList');
         
            bad_hiamp = regionprops('table',...
                (abs(thisEEG-thisEEG_ma))>thresh_high,...
                'Area','BoundingBox','PixelIdxList');
            
            thisResult = bad_hiamp;
            count = size(bad_hiamp,1);
            thisResult.type = repmat({'high amplitude'},count,1);
            thisResult.typeNum = ones(count,1)*0;
            thisResult.channel = repmat(channels(j),count,1);
            thisResult.channelNum = repmat(j,count,1);
            thisResult.channelEventNum = num2cell(1:count)';
            
            artifactEvents = addToEventsList(thisResult, artifactEvents);
            
            bad_lovolt = regionprops('table',...
                (rms(thisEEG)<thresh_low),...
                'Area','BoundingBox','PixelIdxList');
            
            thisResult = bad_lovolt;
            count = size(bad_lovolt,1);
            thisResult.type = repmat({'lovolt'},count,1);
            thisResult.typeNum = ones(count,1)*1;
            thisResult.channel = repmat(channels(j),count,1);
            thisResult.channelNum = repmat(j,count,1);
            thisResult.channelEventNum = num2cell(1:count)';
 
            artifactEvents = addToEventsList(thisResult, artifactEvents);
            
%             
%             bad_theta = regionprops('table',...
%                 (rms(thisEEG)<thresh_low),...
%                 'Area','BoundingBox','PixelIdxList');
%             
%             thisResult = bad_lovolt;
%             count = size(bad_lovolt,1);
%             thisResult.type = repmat({'lovolt'},count,1);
%             thisResult.typeNum = ones(count,1)*1;
%             thisResult.channel = repmat(channels(j),count,1);
%             thisResult.channelNum = repmat(j,count,1);
%             thisResult.channelEventNum = num2cell(1:count)';
%  
%             artifactEvents = addToEventsList(thisResult, artifactEvents);
            
        end
        
          %areas of low voltage 
        
%   
%         
%         % high-pass filter
%         eeg_bp = eegfilt(seg_avg_bi,Fs,12,0);   % from EEG Toolbox
%                 
%         eeg_bp_smooth = nan(size(eeg_bp));
%         
% %         bipolarEvents = cell(size(eeg_bp,1),1);
%         bipolarEvents = table();
%         
%         for j=1:size(eeg_bp,1)
%            thisEEG_bp_smooth = smooth((eeg_bp(j,:).^2)',smoothWin)';
%            eeg_bp_smooth(j,:) = thisEEG_bp_smooth;
%            lpds_ = regionprops('table',thisEEG_bp_smooth>thresh_mag,'Area','BoundingBox','PixelIdxList');
%            
%            thisResult = lpds_;
%            thisResult.channel = repmat(channels(j),size(lpds_,1),1);
%            thisResult.channelNum = repmat(j,size(lpds_,1),1);
%            thisResult.channelEventNum = num2cell(1:size(lpds_,1))';
%            bipolarEvents = cat(1,bipolarEvents, thisResult);
%         end
%         
%         bipolarEvents.eventNum = num2cell(1:size(bipolarEvents,1))';
%         
%         LOC = fct_LOCfromEventsTbl(channels, bipolarEvents); % just include LOC for size
%         
        
                LOC = fct_LOCfromEventsTbl(channels, artifactEvents); 
        
end

