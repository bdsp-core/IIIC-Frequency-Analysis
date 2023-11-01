function [ LOC ] = fct_LOCfromEvents( blankLOC, bipolarEvents)
%fct_LOCfromEvents Takes a cell array of events and returns a nxm
%logical/double array indicating the location of all events in an epoch
% INPUT 
%   blankLOC, any nxm array, for size info
%   bipolarEvents, cell array, contains stat objects from regionprops
% OUTPUT
%   LOC, nxm logical/double array, indicating location of events (channel, time) 
%   in an epoch

LOC = zeros(size(blankLOC));

for j=1:size(LOC,1)
    lpds_= bipolarEvents{j};
    
    lpds_Idx = cat(1,lpds_(:).PixelIdxList);
    LOC(j,lpds_Idx) = 1;
end

end

