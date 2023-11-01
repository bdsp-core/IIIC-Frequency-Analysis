function [spatialExtent] = fct_calculateSpatialExtent(groupEvents_bi_gt0,c_bi);
%fct_calculateSpatialExtent  
%   Calculate the spatial extent (as percentage of bipolar channels)
% INPUT
%   groupEvents_bi_gt0,     a table of events (rows = channels) with channel statistics
%   c_bi,       a cell array of bipolar channel pairs 
%
% OUTPUT
%   spatialExtent,  a double in the range (0-1) representing % of channels involved
%
% N.b. In theory, bipolar events over-represents frontal/occipital channels -- events
%   from average montage may be less biased. 

channels = groupEvents_bi_gt0.channel;

% leaving open possibilitry to filter channels by some minimum criteria (e.g. # of events or relative power)

spatialExtent = sum(ismember(channels,c_bi))/length(c_bi);


end

