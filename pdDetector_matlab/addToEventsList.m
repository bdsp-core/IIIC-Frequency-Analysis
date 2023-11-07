function [ artifactEvents_done ] = addToEventsList( thisResult, artifactEvents_old)
%
% This code addresses a matlab bug. need to insulate all PixelIdxList or PixelValues values as
% cells, it should have done it when it created the table from regionprops,
% but the function is new and it didn't.

if any(ismember(tableNames(thisResult),'PixelIdxList'))
if ~iscell(thisResult.PixelIdxList)
    
    newResult = thisResult;
    newResult.PixelIdxList = [];
    newResult.PixelIdxList = cell(size(thisResult,1),1);
    
    for k = 1:size(thisResult,1)
        newResult(k,:).PixelIdxList = {thisResult(k,:).PixelIdxList}; 
    end   
    thisResult = newResult;
end
end

% review the modified result, having been saved back into 'thisResult'
%
% this is actually a different bug -- the values have been entered directly
% as double by regionprops.
if any(ismember(tableNames(thisResult),'PixelValues'))
if ~iscell(thisResult.PixelValues)
    clear newResult;
    newResult = thisResult;
    newResult.PixelValues  = [];
    newResult.PixelValues  = cell(size(thisResult,1),1);
    
    for k = 1:size(thisResult,1)
        newResult(k,:).PixelValues = {thisResult(k,:).PixelValues};
    end
    
    thisResult = newResult;
end    
end

if any(ismember(tableNames(thisResult),'PixelValues_raw'))
if ~iscell(thisResult.PixelValues_raw)
    clear newResult;
    newResult = thisResult;
    newResult.PixelValues_raw= [];
    newResult.PixelValues_raw= cell(size(thisResult,1),1);
    
    for k = 1:size(thisResult,1)
        newResult(k,:).PixelValues_raw= {thisResult(k,:).PixelValues_raw};
    end
    
    thisResult = newResult;
end    
end
artifactEvents_done = cat(1,artifactEvents_old, thisResult);

end

