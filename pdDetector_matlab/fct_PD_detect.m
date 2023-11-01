function  [freq, channels,interpretation,data_obj] = fct_PD_detect( seg )
%fct_PD_detect
% INPUT  eeg, array of double 
% OUTPUT
%   freq, double
%   channels, cell array of strings, specifying the channel name from bipolar montage (eg. "Fp1-F7")
%
% Chris McGraw. 5/1/2019.

% seg = eeg;
seg_bi = fcn_LBipolar(seg);

gap = NaN(1, size(seg , 2));
Fs = 200;
winSize = 5;
stepSize = 20;
smoothWinSize = 120;
thrDur = 40;

debugFlag = false;

channels_L = [1:4,9:12];
channels_R = [5:8, 13:16];

channelLabels_bipolar = {'Fp1-F7';'F7-T3';'T3-T5';'T5-O1'; 'Fp2-F8';'F8-T4';'T4-T6';'T6-O2'; 'Fp1-F3';'F3-C3';'C3-P3';'P3-O1'; 'Fp2-F4';'F4-C4';'C4-P4';'P4-O2'; 'Fz-Cz';'Cz-Pz'};
channelLabels_bipolar_withspace = {'Fp1-F7';'F7-T3';'T3-T5';'T5-O1'; '';'Fp2-F8';'F8-T4';'T4-T6';'T6-O2'; '';'Fp1-F3';'F3-C3';'C3-P3';'P3-O1'; '';'Fp2-F4';'F4-C4';'C4-P4';'P4-O2'; '';'Fz-Cz';'Cz-Pz'};
c_bi = channelLabels_bipolar;

thresh_hi = 80;
thresh_lo = 75;

smoothWin = Fs/10;
boostWin = Fs/10;
meanWin = 80;

% for PD
freqL = 2;
freqH = 14;

% for RDA
% freqL = 1;
% freqH = 3;

totalSamples = size(seg,2);

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, freqL,freqH, boostWin, meanWin, c_bi);

% criteria = {'(e.MaxIntensity)>18','(e.MaxIntensity_raw)-(e.MinIntensity_raw)>25',...
%     '(e.MaxIntensity_raw<100)','(e.MinIntensity_raw>-100)'};

% for PD
% criteria = {'(e.MaxIntensity)>18','(e.MaxIntensity_raw<170)','(e.MinIntensity_raw>-170)',...
%             '(e.MaxIntensity_raw)-(e.MinIntensity_raw)>20'};

% for RDA
% criteria = {'(e.MaxIntensity)>5','(e.MaxIntensity_raw<150)','(e.MinIntensity_raw>-150)',...
%             '(e.MaxIntensity_raw)-(e.MinIntensity_raw)>25'};

% for seizures
criteria = {'(e.MaxIntensity)>18','(e.MaxIntensity_raw)-(e.MinIntensity_raw)>20'};

pdEvents_bi = selectEvents(c_bi, bipolarEvents, criteria);

if debugFlag
    fct_vizEEG_bipolar(eeg_bp_smooth);
% pd_LOC = fct_LOCfromEventsTbl(c_bi, pdEvents_bi, totalSamples);
% fct_vizLOC_bipolar(pd_LOC);
end

groupEvents_bi = fct_findGroupEvents_bi(pdEvents_bi, seg_bi, eeg_bp_smooth, c_bi);

%% logic to select side
% 
% IDI = groupEvents_bi(ismember(groupEvents_bi.channel,'Global_kmeans'),:);
% freq = IDI.thisIDI_Hz_med;
% channels = unique(pdEvents_bi.channel);

% using thresholded values
theseAssays_group = {{'Left_t2';'Right_t2';'Global_t4'},...
                     {'Left_t0';'Right_t0';'Global_t0'}};

events_t4 = groupEvents_bi.events(find(strcmp(groupEvents_bi.channel,'Global_t4')));
events_t0 = groupEvents_bi.events(find(strcmp(groupEvents_bi.channel,'Global_t0')));

if debugFlag
    figure(45); fct_showEEG_wEvents(seg_bi, pdEvents_bi);
end
% determine assay by # of events
if events_t0 == 0
    % no detected events across all channels
    freq = 0;
    channels = cell(0,1);
    interpretation = {'no events'};
    
    %make a null data_obj
    data_obj = NaN; 
    return
elseif events_t0 == 1 && events_t4 ==0 
    % a single event across all channels
    freq = 1/(totalSamples/Fs);
    groupEvents_bi_byChan = groupEvents_bi(1:18,:);
    groupEvents_bi_gt0 = groupEvents_bi_byChan(groupEvents_bi_byChan.events>0,:);
    channels = groupEvents_bi_gt0.channel;
    interpretation = {'one event'};
    %make a null data_obj
    data_obj = NaN; 
    return
elseif events_t4 == 0
    % more than 1 event at theshold=0, but no events at threshold =4 --> use lowest threshold assays
    theseAssays = theseAssays_group{2};
else
    % events even at theshold =4 --> use intermediate threshold assays
    theseAssays = theseAssays_group{1};
end
    
theseGE_b = groupEvents_bi(ismember(groupEvents_bi.channel, theseAssays),:);
[theseGE_b_pdr_max, indx] = max(theseGE_b.power_discharges_rel);

 pow_G =theseGE_b.power_discharges_rel(3);
 pow_R =theseGE_b.power_discharges_rel(2);
 pow_L = theseGE_b.power_discharges_rel(1);
 relPLvPG = abs((pow_L-pow_G)/pow_G);
 relPRvPG = abs((pow_R-pow_G)/pow_G);
%  relPLvPR = abs((pow_L-pow_R)/mean([pow_L,pow_R]));
 relPLvPR = abs(log(pow_L/pow_R));
 
if debugFlag 
 pow_L
 pow_R
 pow_G
 relPLvPG
 relPRvPG 
 relPLvPR 
end

stats_obj.pow_L = pow_L;
stats_obj.pow_R = pow_R;
stats_obj.pow_G = pow_G;
stats_obj.relPLvPG = relPLvPG;
stats_obj.relPRvPG = relPRvPG; 
stats_obj.relPLvPR = relPLvPR; 

% which set of discharges (left, right, global) has the highest "relative power" (ie. greatest percentage of total BP power that is explained by the dischanges)
switch indx
    case 3
        % report as GPD
        % all channels with events >1
        groupEvents_bi_byChan = groupEvents_bi(1:18,:);
        groupEvents_bi_gt0 = groupEvents_bi_byChan(groupEvents_bi_byChan.events>0,:);
        channels = groupEvents_bi_gt0.channel;
        if relPLvPR >0.08
            % bilateral asymmetric
            if pow_L>pow_R
                interpretation = {'GPD-bilateral asym, L>R'};
            else
                interpretation = {'GPD-bilateral asym, R>L'};
            end
        else
            % global
            interpretation = {'GPD'};
        end

        thisLOC = theseGE_b.LOC{3};  % based on indx, ~~ GPD/Global
%         thisLOC = theseGE_b.LOC{3};  % based on indx, ~~ GPD/Global
        thisEvolution_slope = theseGE_b.evolution_slope(3);  
        thisEvolution_RMSE = theseGE_b.evolution_RMSE(3);  

    case 1
        % report as Left LPD
        % all left channels with events >1
      
        % unless the relative difference in relative power b/w L vs G and R vs G is too similar    
%         if abs((pow_L-pow_G)/pow_G) <0.1 & abs((pow_R-pow_G)/pow_G) <0.2
        if (relPLvPG <0.1 && relPLvPR<0.15)
            % then still report as GPD
            groupEvents_bi_byChan = groupEvents_bi(1:18,:);
            groupEvents_bi_gt0 = groupEvents_bi_byChan(groupEvents_bi_byChan.events>0,:);
            channels = groupEvents_bi_gt0.channel;
            interpretation = {'GPD'};
            indx=3; % redefine as GPD
            thisLOC = theseGE_b.LOC{3};  % based on indx, ~~ GPD/Global
            thisEvolution_slope = theseGE_b.evolution_slope(3);  
            thisEvolution_RMSE = theseGE_b.evolution_RMSE(3);  

        else
            % otherwise report as Left LPD
            groupEvents_bi_byChan = groupEvents_bi([channels_L,17,18],:);
            groupEvents_bi_gt0 = groupEvents_bi_byChan(groupEvents_bi_byChan.events>0,:);
            channels = groupEvents_bi_gt0.channel;
            interpretation = {'L LPD'};
            thisLOC = theseGE_b.LOC{1};  % based on indx, ~~ L LPD
            thisEvolution_slope = theseGE_b.evolution_slope(1);  
            thisEvolution_RMSE = theseGE_b.evolution_RMSE(1);  
        
        end
        
    case 2 
        % report as Right LPD
        % all right channels with events >1
        
         % unless the relative difference in relative power b/w L vs G and R vs G is too similar    
%         if abs((pow_R-pow_G)/pow_G) <0.1 & abs((pow_L-pow_G)/pow_G) <0.2 
        if (relPRvPG <0.1 && relPLvPR<0.15)
            % then still report as GPD
            groupEvents_bi_byChan = groupEvents_bi(1:18,:);
            groupEvents_bi_gt0 = groupEvents_bi_byChan(groupEvents_bi_byChan.events>0,:);
            channels = groupEvents_bi_gt0.channel;
            interpretation = {'GPD'};
            indx=3; % redefine as GPD
            thisLOC = theseGE_b.LOC{3};  % based on indx, ~~ GPD/Global
            thisEvolution_slope = theseGE_b.evolution_slope(3);  
            thisEvolution_RMSE = theseGE_b.evolution_RMSE(3);  

        else

            %otherwise report as Right LPD
            groupEvents_bi_byChan = groupEvents_bi([channels_R,17,18],:);
            groupEvents_bi_gt0 = groupEvents_bi_byChan(groupEvents_bi_byChan.events>0,:);
            channels = groupEvents_bi_gt0.channel;
            interpretation = {'R LPD'};
            thisLOC = theseGE_b.LOC{2};  % based on indx, ~~ GPD/Global
            thisEvolution_slope = theseGE_b.evolution_slope(2);  
            thisEvolution_RMSE = theseGE_b.evolution_RMSE(2);  

        end
end
   

% Calculate frequency (Hz) based on the median of the IDI for detected events, unless 
% the % near the median is <75%, in which case use the frequency calculated
% from cross-correlation based method.

freq(1) = theseGE_b(indx,:).thisIDI_xcorr_Hz;
freq(2) = theseGE_b(indx,:).thisIDI_Hz_med;

if debugFlag
disp('xcorr IDI'); 
freq(1) 
disp('med IDI');
freq(2)
disp('event rate');
theseGE_b(indx,:).eventRate
end

if isnan(freq(1)) & isnan(freq(2))
    freq = theseGE_b(indx,:).eventRate;
elseif ~isnan(freq(2)) & ( (theseGE_b(indx,:).thisIDI_prctNearMed>=0.75) | isnan(freq(1)))
    % median of IDI
    freq = freq(2);
%     disp('med IDI');
else
    %cross-correlation    
    freq = freq(1);
%     disp('xcorr IDI');
end

if isnan(freq) 
    freq = 0;  %suppress NaN but should arise.
end

spatialExtent = fct_calculateSpatialExtent(groupEvents_bi_gt0,c_bi);

if debugFlag
%     figure(46); fct_showEEG_wEventsAndStats(seg_bi, pdEvents_bi, thisLOC, stats_obj);
    figure(47); fct_showEEG_wEventsAndStats(seg_bi, fct_cullEventsByChan(pdEvents_bi,unique(groupEvents_bi_gt0.channel)), thisLOC, stats_obj);
end

data_obj = struct;
data_obj.seg_bi = seg_bi;
data_obj.thisLOC = thisLOC;
data_obj.pdEvents_bi = pdEvents_bi;
data_obj.pdEvents_bi_select = fct_cullEventsByChan(pdEvents_bi,unique(groupEvents_bi_gt0.channel));
data_obj.groupEvents_bi = groupEvents_bi; 
data_obj.groupEvents_bi_gt0 = groupEvents_bi_gt0;
data_obj.thisEvolution_slope = thisEvolution_slope;
data_obj.thisEvolution_RMSE = thisEvolution_RMSE; 
data_obj.freq = freq;
data_obj.channels = channels;
data_obj.interpretation = interpretation;
data_obj.spatialExtent = spatialExtent;
data_obj.stats_obj = stats_obj;

end

