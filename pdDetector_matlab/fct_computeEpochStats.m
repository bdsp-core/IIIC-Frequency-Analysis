function [ epochObject ] = fct_computeEpochStats( seg_bi, Fs , channelLabels)
%fct_computeEpochStats 

% channelLabels = c_bi;
channels = channelLabels;
% channelNums = num2cell(1:length(channelLabels));
% channelNums{end+1} = [1:length(channelLabels)];

amp_tbl =table;

%% amplitude
for j=1:size(seg_bi,1)
    
    thisEEG=seg_bi(j,:);
    
    amp_range = max(thisEEG,[],2)-min(thisEEG,[],2);
    amp_rms = rms(thisEEG,2);

    thisResult = table;
    thisResult.amp_range = amp_range;
    thisResult.amp_rms = amp_rms; 
    
%     thisResult.channel = repmat(channels(j),size(thisResult,1),1);
%     thisResult.channelNum = repmat(channelNums(j),size(thisResult,1),1);
    thisResult.channel = repmat(channels(j),size(thisResult,1),1);
    thisResult.channelNum = repmat(j,size(thisResult,1),1);

    amp_tbl = addToEventsList(thisResult, amp_tbl);
end

thisResult =table;
thisResult.amp_range = mean(amp_tbl.amp_range);
thisResult.amp_rms = mean(amp_tbl.amp_rms);
thisResult.channel = repmat('average',size(thisResult,1),1);
thisResult.channelNum = repmat(nan,size(thisResult,1),1);

amp_tbl = addToEventsList(thisResult, amp_tbl);
epochObject.amplitude = amp_tbl;

%% bkgd frequencies

freqLabels = {'wideband','delta','theta','alpha','beta'};
freqRanges = {[0.5 30],[0.5, 4],[4, 8],[8,12],[12, 20]};

power_tbl = table;

for k=1:length(freqLabels)
    
    power_tbl_p = table;
    
    for j=1:size(seg_bi,1)
       
          thisEEG=seg_bi(j,:);
          
          p=bandpower(thisEEG,Fs,freqRanges{k});
        
          thisResult = table;
          
          thisResult.absPower = p;
          if k==1
          thisResult.relPower = 1;
          else
              thisResult.relPower = p/power_tbl.absPower(j);  % relies on data from wideband being first
          end
          thisResult.freqBand = freqLabels(k);
          thisResult.freqBandNum = k;
          thisResult.channel = repmat(channels(j),size(thisResult,1),1);
          thisResult.channelNum = repmat(j,size(thisResult,1),1);

          power_tbl_p = addToEventsList(thisResult, power_tbl_p);
    end
    
    thisResult =table;
    thisResult.absPower = mean(power_tbl_p.absPower);
    thisResult.relPower = mean(power_tbl_p.relPower);
    thisResult.freqBand = freqLabels(k);
    thisResult.freqBandNum = k;
    thisResult.channel = repmat('average',size(thisResult,1),1);
    thisResult.channelNum = repmat(nan,size(thisResult,1),1);
    
    power_tbl_p = addToEventsList(thisResult, power_tbl_p);
    
    power_tbl = addToEventsList(power_tbl_p, power_tbl);
    
    
end

epochObject.power = power_tbl;

%% artifacts


%% cross-correlation 

% within channels

% within hemisphere

% global

%raw signal
%bp signal


%% automatic segmentation based on BP

%% Vector events

%% Discharge segmentation -- BP

%% Discharge segmenation -- NLEO


end

