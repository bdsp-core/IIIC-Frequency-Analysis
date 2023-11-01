function [pdEvents_bi_s] = fct_cullEventsByChan(pdEvents_bi,c_bi)
%fct_cullEventsByChan
%   Select only events from event table that occur in channels found in the list of channels
% INPUT
%   pdEvents_bi,     a table of channel events (rows = events) 
%   c_bi,       a cell array of bipolar channel pairs 
%
% OUTPUT
%   pdEvents_bi_s,  a table of channel events that includes only those involving the specified channel names
%

pdEvents_bi_s = pdEvents_bi(ismember(pdEvents_bi.channel,c_bi),:);

end

