function [ theseEvents] = fct_findGroupEvents_avg( avgEvents, seg_avg, eeg_bp_smooth, channels)
%fct_findGroupEvents 

channelNums = 1:length(channels);
channelLabels = channels;
Fs = 200;
totalSamples = size(seg_avg,2);

% for each channel, calculate IDI, periodicity %, cross-correlation,
% discharge power
theseResults = table;

% by individual channel
for j=1:length(channelNums)
    
    thisChan = avgEvents(avgEvents.channelNum==j,:);   % table of m rows, each row = an event
    
    if(isempty(thisChan))
        
        thisChan = makeNanTable(thisChan);
        thisChan.channel = channels(j);
        thisChan.channelNum = channelNums(j);
        
    end
    
    thisEEG_BP_smooth= eeg_bp_smooth(j,:);     % time-series of BP power for this channel
    thisResult = fct_calcPeriodicity(thisChan, thisEEG_BP_smooth,Fs);
    theseResults = cat(1, theseResults,thisResult);
end

% by hemisphere (L,R,C) -- for average montage.
theseGroups_labels = {'Left_t0','Right_t0','Central_t0','Global_t0', ...
    'Left_t2','Right_t2','Central_t1','Global_t4'};
theseGroups_channelNums = {[1:8],[12:19],[9:11],[1:19],[1:8],[12:19],[9:11],[1:19]};
theseGroups_t = {0,0,0,0,2,2,1,4};

for j=1:length(theseGroups_labels)
    
    thisGroup_label = theseGroups_labels(j);
    thisGroups_channelNums = theseGroups_channelNums{j};
    
    theseChan = avgEvents(ismember(avgEvents.channelNum,thisGroups_channelNums),:);
    
    theseLOC = fct_LOCfromEventsTbl(thisGroups_channelNums,theseChan,totalSamples);
    thisLOC = sum(theseLOC,1);
    %     figure;imagesc(thisLOC);
    thisEEG_BP_smooth = sum(eeg_bp_smooth(thisGroups_channelNums,:),1);
    
    [~,thisChan] = fct_findEvents_simple(thisLOC>theseGroups_t{j},thisEEG_BP_smooth,thisGroup_label);
    
    thisResult = fct_calcPeriodicity(thisChan, thisEEG_BP_smooth,Fs);
    theseResults = cat(1, theseResults,thisResult);
end

% by hemisphere (L,R,C) using kmeans to segment
theseGroups_labels = {'Left_kmeans','Right_kmeans','Central_kmeans','Global_kmeans'};
theseGroups_channelNums = {[1:8],[12:19],[9:11],[1:19]};
% theseGroups_t = {0,0,0,0,2,2,1,4};

for j=1:length(theseGroups_labels)
    
    thisGroup_label = theseGroups_labels(j);
    thisGroups_channelNums = theseGroups_channelNums{j};
    
    theseChan = avgEvents(ismember(avgEvents.channelNum,thisGroups_channelNums),:);
    
    theseLOC = fct_LOCfromEventsTbl(thisGroups_channelNums,theseChan,totalSamples);
    thisLOC = sum(theseLOC,1);
    %     figure;imagesc(thisLOC);
    thisEEG_BP_smooth = sum(eeg_bp_smooth(thisGroups_channelNums,:),1);
    thisEEG_kmeans = kmeans_sane(thisEEG_BP_smooth,3);
    
    [~,thisChan] = fct_findEvents_simple(thisEEG_kmeans>0,thisEEG_BP_smooth,thisGroup_label);
    
    thisResult = fct_calcPeriodicity(thisChan, thisEEG_BP_smooth, Fs);
    theseResults = cat(1, theseResults,thisResult);
end

% kmeans to choose globally but then report by hemisphere

theseGroups_labels = {'_global_kmeans'};
theseSubGroup_labels = {'Left','Right'};
theseGroups_channelNums = {[1:19]};
theseSubGroup_channelNums = {[1:8],[12:19]};

for j=1:length(theseGroups_labels)
    
    thisGroup_label = theseGroups_labels(j);
    thisGroups_channelNums = theseGroups_channelNums{j};
    
     eeg_bp_smooth_L = sum(eeg_bp_smooth(theseSubGroup_channelNums{1},:),1);
     eeg_bp_smooth_R = sum(eeg_bp_smooth(theseSubGroup_channelNums{2},:),1);
     
     thisEEG_kmeans = kmeans_sane([eeg_bp_smooth_L; eeg_bp_smooth_R],3);   %prev 2
     
     [~,thisChan] = fct_findEvents_simple(thisEEG_kmeans(1,:)>0,eeg_bp_smooth_L,...
         {[theseSubGroup_labels{1} thisGroup_label{:}]});
     
     thisResult = fct_calcPeriodicity(thisChan, eeg_bp_smooth_L, Fs);
     theseResults = cat(1, theseResults,thisResult);
     
     [~,thisChan] = fct_findEvents_simple(thisEEG_kmeans(2,:)>0,eeg_bp_smooth_R,...
         {[theseSubGroup_labels{2} thisGroup_label{:}]});
     
     thisResult = fct_calcPeriodicity(thisChan, eeg_bp_smooth_R,Fs);
     theseResults = cat(1, theseResults,thisResult);
    
end

theseEvents = theseResults;
% average by hemisphere, and global 


% for each hemisphere, sum events over channels (summed LOC), calculate
% IDI, periodicity%, cross-correlation, discharge power

    % send to 'findEvents' routine?

% for each hemisphere, using BP eeg, calculate IDI, periodicity%,
% cross-correlation, discharge power

    % send to 'findEvents' routine?
    
    
% combine results and return table



end