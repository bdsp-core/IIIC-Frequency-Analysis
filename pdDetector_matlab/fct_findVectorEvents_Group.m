function [ eeg_hdv ] = fct_findVectorEvents_Group( seg, Fs)


channels_L = [1:4,9:12];
channels_R = [5:8, 13:16];

thresh_hi = 0;
thresh_lo = 0;

% LOC = fct_LOCfromEventsTbl({1},vectorEvents(strcmp(vectorEvents.type,'vector min'),:));
% figure(37);imagesc(LOC)

thisPatternBlock = seg(:,1:80);
theseEEGChannels = seg;
% 
% thisPatternBlock = seg_bi(channels_L,theseEvents_locs);
% theseEEGChannels = seg_bi(channels_L,:);

eeg_hdv = zeros(size(theseEEGChannels,2),9*size(thisPatternBlock,2));
blockSize = size(thisPatternBlock,2);

vectorEvents_types = { 'vector min', 'vector max', 'vector magnitude'};

vectorEvents_toInclude = {{'PixelValues','PixelValues_min_x','PixelValues_min_y'},...
    {'PixelValues','PixelValues_max_x','PixelValues_max_y'},...
    {'PixelValues','PixelValues_max_x','PixelValues_max_y'}};

for j=1:size(theseEEGChannels,2)-blockSize
    
    
    vectorEvents= fct_findVectorEvents(seg(:,j:j+blockSize-1), Fs, thresh_hi, thresh_lo); %thresh_hi, thresh_lo, smoothWin, c_av)
    
    theseTypes = unique(vectorEvents.type);
    %     thisBlock = theseEEGChannels(:,j:j+blockSize-1);
    
    thisBlock = zeros(1,blockSize*9);
    
    for k= 1:length(theseTypes)
        
        thisType = ismember(vectorEvents_types,theseTypes(k));
        theseData_toInclude = vectorEvents_toInclude{thisType};
        
        for m = 1:length(theseData_toInclude)
            % vectorize all the data from vectorEvents
                                   
            thisData_byType = vectorEvents(ismember(vectorEvents.type,theseTypes(k)),:);
                                    
            thisStart = (k-1)*(80*3) + (m-1)*80+1;
            thisEnd = (k-1)*(80*3) + (m*80);
            thisBlock(1,thisStart:thisEnd) = thisData_byType.(theseData_toInclude{m});
            
        end
    end
    
    
%     eeg_hdv(j,:) = reshape(thisBlock,1,size(thisBlock,1)*size(thisBlock,2));
    eeg_hdv(j,:) = thisBlock;
    
end

end

