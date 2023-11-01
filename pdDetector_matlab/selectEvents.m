function [ theseEvents_new ] = selectEvents( c_bi, bipolarEvents, criteria, varargin)
%selectEvents Takes event table and returns new table only containing
%events based on criteria. 
% INPUT
%   c_bi, cell array of channels in montage
%   bipolarEvents, nxm table of events (can be any montage)
%   criteria, cell array of strings containing comparative expressions to be evaluated
%   by Matlab using eval() command, referencing columns in bipolarEvents,
%   e.g. '(e.Area>30)' or '~(e.Area>30)'. Use e.*** syntax 
% OUTPUT
%   theseEvents_new, (n-y)xm table containing a subset of the original n
%   rows according to specified criteria

e = bipolarEvents; 
Fs =200;  %used by eval expressions

if size(varargin,1)>0
    
    for k=1:size(varargin,1)
        inVar{k} = varargin{k}; % ignore the rest
    end
    
    if size(inVar,1)==1
        inVar = inVar{:};    % if only one element, unpackage cell
    end
end

theseEvents_new = table;

theseIndx = logical(ones(size(e,1),1));

% combine criteria by Boolean logic
for j=1:length(criteria)
    thisIndx = eval(criteria{j});
    
    if iscell(thisIndx)
        thisIndx = cat(1,thisIndx{:});
    end
    
    theseIndx = [theseIndx & thisIndx];
end

theseEvents_new = e(theseIndx,:);

end

