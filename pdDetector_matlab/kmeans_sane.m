function [ data_km ] = kmeans_sane( rawData, varargin)
%kmeans_sane Wrapper for kmeans function that makes sure that kmeans output labels the data with highest
%number, and background as zero, and returns in same dimensions as provided. 

% run kmeans passing other variables, and reshape output back to original
% dimensions of the rawData. Subtract 1 to make the lowest value be 0
% instead of 1.
% data_km = reshape(kmeans(rawData(:), varargin{:})-1,repmat(1,size(rawData(:))),size(rawData));
data_km = reshape( kmeans(rawData(:), varargin{:})-1, size(rawData));


theseKVals = unique(data_km);

% using heuristic, determine background signal vs noise based simply on
% sum(abs(x)) per pixel -- whichever is greater should have the higher value. 
indx_min = find(data_km==min(theseKVals));  % linear indices, not logical
data_km_min = abs(sum(rawData(indx_min)))/length(indx_min);
indx_max = find(data_km==max(theseKVals));
data_km_max = abs(sum(rawData(indx_max)))/length(indx_max);

if data_km_max<data_km_min
    data_km = data_km*-1+repmat(max(theseKVals),size(data_km));
end


end


