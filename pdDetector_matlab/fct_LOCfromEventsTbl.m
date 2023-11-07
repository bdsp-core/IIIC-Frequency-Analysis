function [ LOC ] = fct_LOCfromEventsTbl(channels, bipolarEvents, totalSamples)
%fct_LOCfromEventsTbl Takes a table array of events and returns a nxm
%logical/double array indicating the location of all events in an epoch.
%Can process tables with missing channels, etc -- leaving blank spots in
%the LOC image.
% INPUT
%   channels,   cell array of strings of channel names for montage
%   bipolarEvents, *table* array, contains stat objects from regionprops
%       ** events table must contain a column called 'channelNum'
%       indicating the channel numbers.
%   totalSamples,   specify total number of samples per channel
% OUTPUT
%   LOC, nxm logical/double array, indicating location of events (channel, time)
%   in an epoch

channelNums = unique(bipolarEvents.channelNum);
% LOC = blankLOC;
% totalSamples = 2800;
% totalSamples = 3200;

LOC = zeros(length(channels),totalSamples); % cannot rely on channelNums from bipolarEvents to be complete

for j=1:length(channelNums)
    
    lpds_tbl= bipolarEvents(bipolarEvents.channelNum==channelNums(j),:);
    
    thesePixelIdxList = lpds_tbl.PixelIdxList;
    
    tPIL_s = size(thesePixelIdxList,1);
    
    % if table is a single entry table composed only of NaNs, cannot use brace indexing to detect it.
    if (size(lpds_tbl,1)==1 & tPIL_s ==1 & isnan(lpds_tbl.Area))
        % we found a NaN table -- do nothing
        
    else
        
        if (~isempty(lpds_tbl) & ~isnan(lpds_tbl.PixelIdxList{1}))
            
            lpds_Idx = cat(1,lpds_tbl.PixelIdxList{:});
            LOC(channelNums(j),lpds_Idx) = 1;
        else
            % no events in this channel. Just move on to next one
        end
    end
end

end

