function [LOC, bipolarEvents]= fct_findBipolarEvents(seg_bi, Fs, thresh_mag, smoothWin, channelLabels_bipolar)
% fct_eventDetector_forLPDs
% INPUT seg_bi, nxm double array of EEG in bipolar montage
%       Fs, sampling freq in sec-1
%      ***
% OUTPUT 
%       LOC, nxm logical array corresponding to positions of candidate
%           discharges.
%       bipolarEvents, nx1 cell array of structs containing stat objects
%       from region props
%
           
        LOC = zeros(size(seg_bi));
%         thresh_mag = 15*5;
%         smoothWin = 20;
%                 
        % band-pass filter
        eeg_bp = eegfilt(seg_bi,Fs,4,20);   % from EEG Toolbox
        eeg_bp_smooth = nan(size(eeg_bp));
        
%         bipolarEvents = cell(size(eeg_bp,1),1);
        bipolarEvents = table();
        
%         for j=1:size(eeg_bp,1)
%            thisEEG_bp_smooth = smooth((eeg_bp(j,:).^2)',smoothWin)';
%            eeg_bp_smooth(j,:) = thisEEG_bp_smooth;
%            lpds_ = regionprops(thisEEG_bp_smooth>thresh_mag,'Area','BoundingBox','PixelIdxList');
%            
%            lpds_Idx = cat(1,lpds_(:).PixelIdxList);
%            LOC(j,lpds_Idx) = 1;
%            bipolarEvents{j} = lpds_;
%         end
%         

        for j=1:size(eeg_bp,1)
           thisEEG_bp_smooth = smooth((eeg_bp(j,:).^2)',smoothWin)';
           eeg_bp_smooth(j,:) = thisEEG_bp_smooth;
           lpds_ = regionprops('table',thisEEG_bp_smooth>thresh_mag,'Area','BoundingBox','PixelIdxList');
           
           thisResult = lpds_;
           thisResult.channel = channels{j}
           
           bipolarEvents = cat(1,bipolarEvents, lpds_);
        end
        
        LOC = fct_LOCfromEvents(LOC, bipolarEvents);
        
        figure(1);
        showEEG(eeg_bp_smooth/5);   
        
        figure(2); 
        
        channels_L = [1:4,9:12];
        channels_R = [5:8, 13:16];

        eeg_bp_smooth_L = sum(eeg_bp_smooth(channels_L,:),1);
        eeg_bp_smooth_R = sum(eeg_bp_smooth(channels_R,:),1);
        
        im_eeg_bp = [eeg_bp_smooth_L; eeg_bp_smooth_R];
        ax1=subplot(2,1,1); plot(1:size(im_eeg_bp,2),im_eeg_bp(1,:)); xlim([0 size(eeg_bp_smooth,2)]);
        ax2=subplot(2,1,2); plot(1:size(im_eeg_bp,2), im_eeg_bp(2,:)); xlim([0 size(eeg_bp_smooth,2)]);
        linkaxes([ax1,ax2],'xy');
        
        eeg_label1 = kmeans([eeg_bp_smooth_L eeg_bp_smooth_R]',2);
        figure(3); imagesc([eeg_label1(1:2800),eeg_label1(2801:end)]');
           
        eeg_label2L = kmeans([eeg_bp_smooth_L]',2);
        eeg_label2R = kmeans([eeg_bp_smooth_R]',2);
        figure(11); imagesc([eeg_label2L,eeg_label2R]');
        
        figure(12); showEEG(seg_bi);
        
        
        coolSide = sum(kurtosis(seg_bi(channels_L,:)')-kurtosis(seg_bi(channels_R,:)'));
        % select group1
        
%         indx_sub1 = find(~logical(eeg_label1-1));
%         
%         eeg_bp_smooth_L_s1 = eeg_bp_smooth_L(indx_sub1);
%         
%         eeg_label2 = kmeans([eeg_bp_smooth_L(find(~(logical(eeg_label1-1))))
%             eeg_bp_smooth_R(find(eq(eeg_label1,1)))]',2);
%         figure(4); imagesc([eeg_label2(1:length(eeg_label2)/2),
%             eeg_label2(length(eeg_label2)/2+1:end)]');
    
    end