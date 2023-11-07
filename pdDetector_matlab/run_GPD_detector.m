%% GPD detector
% Take file, seg. 
dataDirec = 'C:\Users\Chris McGraw\Dropbox (Partners HealthCare)\01 Data\2018\FrequencyOfDischarges_Chris_JJ\code\data\';
subDirec = [dataDirec '\LPD\'];
subDirec = [dataDirec '\GPD\'];
subDirec = [dataDirec '\LRDA\'];
subDirec = [dataDirec '\GRDA\']; 
subDirec = [dataDirec '\Sz\']; 

files = dir([subDirec '*.mat']);

%25 (good lpds, and some noise);
% 8 is ~low voltage
n=31;   
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);
%% Other preliminaries
gap = NaN(1, size(seg , 2));
Fs = 200;
winSize = 5;
stepSize = 20;
smoothWinSize = 120;
thrDur = 40;

%% Generate bipolar events
%  -- bipolarEvents Obj

channels_L = [1:4,9:12];
channels_R = [5:8, 13:16];

channelLabels_bipolar = {'Fp1-F7';'F7-T3';'T3-T5';'T5-O1'; 'Fp2-F8';'F8-T4';'T4-T6';'T6-O2'; 'Fp1-F3';'F3-C3';'C3-P3';'P3-O1'; 'Fp2-F4';'F4-C4';'C4-P4';'P4-O2'; 'Fz-Cz';'Cz-Pz'};
channelLabels_bipolar_withspace = {'Fp1-F7';'F7-T3';'T3-T5';'T5-O1'; '';'Fp2-F8';'F8-T4';'T4-T6';'T6-O2'; '';'Fp1-F3';'F3-C3';'C3-P3';'P3-O1'; '';'Fp2-F4';'F4-C4';'C4-P4';'P4-O2'; '';'Fz-Cz';'Cz-Pz'};
channelLabels_average= {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};
channelLabels_average_withspace = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'';'Fz';'Cz';'Pz';'';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};

c_bi = channelLabels_bipolar;
c_av = channelLabels_average;


% identify epileptiform discharges
thresh_mag = 15*5;
smoothWin = Fs/10;
[LOC, bipolarEvents, eeg_bp_smooth] = fct_findBipolarEvents(seg_bi, Fs, thresh_mag, smoothWin, channelLabels_bipolar);

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents(seg_bi, Fs, thresh_mag, smoothWin, 4, 20, channelLabels_bipolar);

% validate kmeans function -- to determine which label is signal and which
% bkgd. ***

%% Visualize average vs bipolar montage events

n=12;   % 15 = low amplitude, with ?BiPDs?
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

% t_mag_discharge = 15*10;
% smoothWin_discharge = Fs/10;
% [LOC, dischargeEvents, eeg_bp_smooth_discharge] = fct_findEvents(seg_avg, Fs, ...
%     t_mag_discharge, smoothWin_discharge, 4, 20,channelLabels_average);
% fct_vizLOC_average(LOC);
% % fct_showEEG_avg_wEvents(seg_avg, spikeEvents);
% figure(20); fct_showEEG_avg(eeg_bp_smooth_discharge/10);
% figure(21);fct_showEEG_avg_wEvents(seg_avg, dischargeEvents);
% figure(22); fct_showEEG_avg(seg_avg);

thresh_hi = 15*9;
thresh_lo = 15*3.5;
smoothWin = Fs/20;
% [LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents(seg_bi, Fs, thresh_mag, smoothWin, 4, 20, channelLabels_bipolar);
[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4, 20, channelLabels_bipolar);
figure(23); showEEG(seg_bi);
figure(24); fct_showEEG_wEvents(seg_bi, bipolarEvents);
fct_vizLOC_bipolar(LOC);
figure(5); showEEG(eeg_bp_smooth/2);
figure(4); fct_showEEG_wEvents(seg_bi, bipolarEvents);

% select events with specific properties
criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs', 'e.MinIntensity_raw>-100', ...
    'e.MaxIntensity_raw<100'};
% criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs'};
pdEvents = selectEvents(c_bi, bipolarEvents, criteria);
figure(44); fct_showEEG_wEvents(seg_bi, pdEvents);


fct_vizLOC_bipolar(fct_LOCfromEventsTbl(c_bi,pdEvents));
fct_vizEEG_bipolar(eeg_bp_smooth);

%% Visualize bipolar events

figure(12); showEEG(seg_bi);
figure(13); fct_showEEG_avg(seg_avg);
fct_vizLOC_bipolar(LOC);
fct_vizEEG_bipolar(eeg_bp_smooth);
figure(5); showEEG(eeg_bp_smooth/10);

%% Visualize events on raw EEG

figure(4); fct_showEEG_wEvents(seg_bi, bipolarEvents);

%% Flexibly visualize based on which events you would like to see
theseRect = fct_showEEG_wEvents(seg_bi, bipolarEvents(bipolarEvents.Area<Fs/2,:));
fct_vizLOC_bipolar(fct_LOCfromEventsTbl(LOC,bipolarEvents(bipolarEvents.Area<Fs/2,:)));

figure(20); fct_showEEG_wEvents(seg_bi, bipolarEvents(bipolarEvents.Area<0.5*Fs,:));
figure(21); fct_showEEG_wEvents(seg_bi, bipolarEvents((bipolarEvents.Area<0.5*Fs & ...
    bipolarEvents.MaxIntensity>500),:));

figure; histogram(bipolarEvents.MaxIntensity)

%% Using a transformed BP to boost co-occuring signals.

thresh_hi = 15*9;
thresh_lo = 15*3;
smoothWin = Fs/20;
boostWin = Fs/20;  %  prev 3
% [LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents(seg_bi, Fs, thresh_mag, smoothWin, 4, 20, channelLabels_bipolar);
[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4, 20, boostWin, channelLabels_bipolar);
figure(23); showEEG(seg_bi);
figure(24); fct_showEEG_wEvents(seg_bi, bipolarEvents);
fct_vizLOC_bipolar(LOC);
figure(5); showEEG(eeg_bp_smooth/10);
figure(4); fct_showEEG_wEvents(seg_bi, bipolarEvents);

% select events with specific properties
criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs', 'e.MinIntensity_raw>-100', ...
    'e.MaxIntensity_raw<100','e.MaxIntensity_raw-e.MinIntensity_raw>10'};
% criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs'};
pdEvents = selectEvents(c_bi, bipolarEvents, criteria);
figure(44); fct_showEEG_wEvents(seg_bi, pdEvents);


fct_vizLOC_bipolar(fct_LOCfromEventsTbl(c_bi,pdEvents));
fct_vizEEG_bipolar(eeg_bp_smooth);


%% Using a per-channel percentile-based auto threshold 

n=24;   % 15 = low amplitude, with ?BiPDs?
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

thresh_hi = 90;
thresh_lo = 85;
smoothWin = Fs/10;
boostWin = Fs/10;
meanWin = 80;

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto4(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, channelLabels_bipolar);
figure(99); showEEG(seg_bi);
figure(24); fct_showEEG_wEvents(seg_bi, bipolarEvents);
fct_vizLOC_bipolar(LOC);
figure(5); showEEG(eeg_bp_smooth/50);
% figure(4); fct_showEEG_wEvents(seg_bi, bipolarEvents);

% select events with specific properties
% criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs', 'e.MinIntensity_raw>-100', ...
%     'e.MaxIntensity_raw<100'};
% criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs'};
criteria = {'(e.MaxIntensity_raw-e.MinIntensity_raw)>25'};
pdEvents = selectEvents(c_bi, bipolarEvents, criteria);
figure(44); fct_showEEG_wEvents(seg_bi, pdEvents);


% fct_vizLOC_bipolar(fct_LOCfromEventsTbl(c_bi,pdEvents));
fct_vizEEG_bipolar(eeg_bp_smooth);
%% Using a per-channel percentile-based auto threshold with low threshold based on local change in slope

Fs=200;

n=11;   % 15 = low amplitude, with ?BiPDs?; % 20, great LPD_R>L.  28 = artifact
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

thresh_hi = 90;
thresh_lo = 85;
smoothWin = Fs/10;
boostWin = Fs/10;
meanWin = 80;

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, channelLabels_bipolar);
% 
% thresh_hi = 70;
% thresh_lo = 50;
% smoothWin = Fs/2;
% boostWin = Fs/2;
% meanWin = 80;
% 
% [LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_bi, Fs, ...
%     thresh_hi, thresh_lo, smoothWin, 0,4, boostWin, meanWin, channelLabels_bipolar);
% 

figure(99); showEEG(seg_bi);
figure(98); fct_showEEG_avg(seg_avg);
figure(24); fct_showEEG_wEvents(seg_bi, bipolarEvents);
fct_vizLOC_bipolar(LOC);
figure(5); showEEG(eeg_bp_smooth/50);
fct_vizEEG_bipolar(eeg_bp_smooth);

thisFile = files(n);

% fct_makeSuperFigure(seg,seg_bi,LOC,bipolarEvents,eeg_bp_smooth,thisFile);


thresh_hi = 90;
thresh_lo = 85;
smoothWin = Fs/10;
boostWin = Fs/10;
meanWin = 80;

[LOC, avgEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_avg, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_av);

figure(97); fct_showEEG_avg_wEvents(seg_avg, avgEvents);
figure(99); fct_showEEG_avg(seg_avg);
figure(102); fct_showEEG_avg(eeg_bp_smooth/5);
fct_vizEEG_average(eeg_bp_smooth)
fct_vizLOC_average(LOC)


imageAvgEEG(eeg_bp_smooth,seg_bi)
% imageAvgEEG(seg,seg_bi)

% 

% figure(4); fct_showEEG_wEvents(seg_bi, bipolarEvents);

%select events with specific properties
% criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs', 'e.MinIntensity_raw>-100', ...
%     'e.MaxIntensity_raw<100'};
criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs'};
% criteria = {'(e.MaxIntensity_raw-e.MinIntensity_raw)>25'};
criteria = {'(e.MaxIntensity)>100','e.Area < 0.250*Fs'};

pdEvents = selectEvents(c_av, avgEvents, criteria);
figure(44); fct_showEEG_avg_wEvents(seg_avg, pdEvents);


% fct_vizLOC_bipolar(fct_LOCfromEventsTbl(c_bi,pdEvents));
fct_vizEEG_bipolar(eeg_bp_smooth);

%% Group stats 
thresh_hi = 90;
thresh_lo = 85;
smoothWin = Fs/10;
boostWin = Fs/10;
meanWin = 80;


n=41;   % 15 = low amplitude, with ?BiPDs?; % 20, great LPD_R>L. 
%  19 =low amplitude on R with small PLDs, 
% 22= some artifact
%28 = artifact
%39 = noisy on L, R has low amplitude but recurrent blips 
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

% 
% [LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_bi, Fs, ...
%     thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_bi);

[LOC, avgEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_avg, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_av);

figure(96); fct_showEEG_wEvents(seg_bi, bipolarEvents);
figure(97); fct_showEEG_avg_wEvents(seg_avg, avgEvents);
figure(99); fct_showEEG_avg(seg_avg);
figure(100); showEEG(seg_bi);
figure(102); fct_showEEG_avg(eeg_bp_smooth/5);
fct_vizEEG_average(eeg_bp_smooth)
fct_vizLOC_average(LOC)

% criteria = {'(e.MaxIntensity)>50','e.Area < 0.30*Fs'};
criteria = {'(e.MaxIntensity)>40','e.MaxIntensity_raw<100'};

pdEvents = selectEvents(c_av, avgEvents, criteria);
figure(44); fct_showEEG_avg_wEvents(seg_avg, pdEvents);

groupEvents = fct_findGroupEvents(avgEvents, seg_avg, eeg_bp_smooth, c_av);

% groupEvents = fct_findGroupEvents(pdEvents, seg_avg, eeg_bp_smooth, c_av);

%%
thisEvent = fct_calcPeriodicity(pdEvents(pdEvents.channelNum==14,:),eeg_bp_smooth(14,:));
thisEvent = fct_calcPeriodicity(bipolarEvents(bipolarEvents.channelNum==14,:),eeg_bp_smooth(14,:));
thisEvent = fct_calcPeriodicity(bipolarEvents(bipolarEvents.channelNum==14,:),seg_bi(14,:));

thisEvent = fct_calcPeriodicity(bipolarEvents(bipolarEvents.channelNum==14,:),seeg_1Hz_sm(1,:));

figure; plot(lags(5*Fs:end)./Fs,simplema(r(5*Fs:end),10,2)); xlim([0 5]);
figure; plot(eeg_bp_smooth(14,:));

%% Generate synthetic data for cross-correlation;

tlen = size(seg_bi,2);
seeg = zeros(1,tlen);

hits_n= 5;
hits_mean= 130;
hits_sigma = 2;
hits_loc = round(tlen*rand(hits_n,1));
hits_vals = randn(hits_n,1)*hits_sigma + hits_mean;
seeg(1,hits_loc) = hits_vals;

other_n = 50;
other_mean= 35;
other_sigma =2; 
other_loc = round(tlen*rand(other_n,1));
other_vals = randn(other_n,1)*other_sigma + other_mean;
seeg(1,other_loc) = other_vals;

seeg_sm = simplema(simple_blur(seeg,25,2)+(rand(tlen,1)*5)',5,2)
figure;plot(seeg_sm);
figure; plot(eeg_bp_smooth(14,:));

%%
Fs = 200;
tlen = size(seg_bi,2);
seeg_1Hz = zeros(1,tlen);

hits_n= 8;
hits_mean= 130;
hits_sigma = 2;


hits_period = [2*Fs, 1*Fs]; % cycles per second
hits_period_sigma = 15;
hits_period_selector = round(rand(hits_n,1))+1;
hits_init = round(100*rand(1));
hits_loca = round(hits_period_sigma*randn(hits_n,1)+hits_period(hits_period_selector)');
hits_loc = zeros(1,hits_n);

for j=1:length(hits_loc)
    if j==1
    hits_loc(j) = hits_init+hits_loca(j);
    else
        hits_loc(j) = hits_loc(j-1) + hits_loca(j);
    end    
end

hits_vals = randn(hits_n,1)*hits_sigma + hits_mean;
seeg_1Hz(1,hits_loc) = hits_vals';

other_n = 75;
other_mean= 15;
other_sigma =2; 
other_loc = round(tlen*rand(other_n,1));
other_vals = randn(other_n,1)*other_sigma + other_mean;
seeg_1Hz(1,other_loc) = other_vals;

seeg_1Hz_sm = simplema(simple_blur(seeg_1Hz,25,2)+(rand(tlen,1)*5)',5,2)
figure;plot(seeg_1Hz_sm);


thisEvent = fct_calcPeriodicity(bipolarEvents(bipolarEvents.channelNum==14,:),seeg_1Hz_sm(1,:));

%%
figure;
[wt,f]=cwt(eeg_bp_smooth(13,:),Fs);
[wt,f]=cwt(seg_avg(13,:),Fs);


%% Using a per-channel percentile-based auto threshold on NLEO

n=41;   % 15 = low amplitude, with ?BiPDs?
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

thresh_hi = 90;
thresh_lo = 90;
smoothWin = Fs/10;  % typically 120msec = 24 samples
boostWin = Fs/10;
meanWin = 80;

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto4_nleo(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 8,20, boostWin, meanWin, channelLabels_bipolar);
figure(23); showEEG(seg_bi);
figure(24); fct_showEEG_wEvents(seg_bi, bipolarEvents);
fct_vizLOC_bipolar(LOC);
figure(5); showEEG(eeg_bp_smooth);
figure(4); fct_showEEG_wEvents(seg_bi, bipolarEvents);

% select events with specific properties
criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs', 'e.MinIntensity_raw>-100', ...
    'e.MaxIntensity_raw<100'};
% criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs'};
pdEvents = selectEvents(c_bi, bipolarEvents, criteria);
figure(44); fct_showEEG_wEvents(seg_bi, pdEvents);


% fct_vizLOC_bipolar(fct_LOCfromEventsTbl(c_bi,pdEvents));
fct_vizEEG_bipolar(eeg_bp_smooth);
%% Using a per-channel percentile-based auto threshold on NLEO with lower threshold based ...
% on slope == 0. 

n=15;   % 15 = low amplitude, with ?BiPDs?
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

thresh_hi = 90;
thresh_lo = 90;
smoothWin = 30;  % typically 120msec = 24 samples
boostWin = Fs/10;
meanWin = 80;

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6_nleo(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 8,20, boostWin, meanWin, channelLabels_bipolar);
figure(23); showEEG(seg_bi);
figure(24); fct_showEEG_wEvents(seg_bi, bipolarEvents);
fct_vizLOC_bipolar(LOC);
figure(5); showEEG(eeg_bp_smooth);
% figure(4); fct_showEEG_wEvents(seg_bi, bipolarEvents);

% select events with specific properties
criteria = {'e.Area > 0.06*Fs'}
% criteria = {'e.Area > 0.06*Fs', 'e.Area < 0.5*Fs'};
pdEvents = selectEvents(c_bi, bipolarEvents, criteria);
figure(44); fct_showEEG_wEvents(seg_bi, pdEvents);


% fct_vizLOC_bipolar(fct_LOCfromEventsTbl(c_bi,pdEvents));
fct_vizEEG_bipolar(eeg_bp_smooth);
%%
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);   

% select based on line-crossings
criteria_p = {"cellfun(@(x) (length(zci(diff(x)))-1)<9, e.PixelValues_raw, 'UniformOutput', false)"};
criteria_p = {"cellfun(@(x) find(diff(x>0))<3, e.PixelValues_raw, 'UniformOutput', false)"};

% criteria_p = {"cellfun(@(x) sum(eq(diff(x),0)), e.PixelValues_raw, 'UniformOutput', false)"};
criteria_cellfun = join(cat(2,"cat(1,cellfun(@(x) any(x), ", criteria_p{:}, " , 'UniformOutput', false))"));

% pdEvents = selectEvents(c_bi, bipolarEvents, criteria);
% pdEvents = selectEvents(c_bi, bipolarEvents, criteria_p);
pdEvents = selectEvents(c_bi, bipolarEvents, criteria_cellfun);
figure(44); fct_showEEG_wEvents(seg_bi, pdEvents);

zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector


%% Make epoch statistics

epochObject = fct_computeEpochStats(seg_bi, Fs, c_bi); 
% ampltiude
% bkgd
% cross-correlation
%% Grouping

% using BP eeg or NLEO 

%% Generate discharge-specific statistics



bipolarEvents = calc_nleo(c_bi, bipolarEvents)


%% Generate groups of discharges


%% Vector magnitude for naive segmentation
thresh_hi = 80;
thresh_lo = 50;

vectorEvents= fct_findVectorEvents(seg, Fs, thresh_hi, thresh_lo); %thresh_hi, thresh_lo, smoothWin, c_av)

LOC = fct_LOCfromEventsTbl({1},vectorEvents(strcmp(vectorEvents.type,'vector magnitude'),:));
figure(37);imagesc(LOC)

%% Generate events in vector space.  Negative vector

% fct_showEEG_avg(seg_avg);
% imageAvgEEG(seg,seg_bi)
thresh_hi = 80;
thresh_lo = 50;

vectorEvents= fct_findVectorEvents(seg, Fs, thresh_hi, thresh_lo); %thresh_hi, thresh_lo, smoothWin, c_av)

LOC = fct_LOCfromEventsTbl({1},vectorEvents(strcmp(vectorEvents.type,'vector min'),:));
figure(37);imagesc(LOC)

% find x,y pos 
vectorEvents_min = vectorEvents(strcmp(vectorEvents.type,'vector min'),:);
xlen =41;
ylen =41;

LOC_LR = logical(zeros(3,size(seg_bi,2)));
LOC_AP = zeros(3,size(seg_bi,2));
LOC_LRAP = zeros(9, size(seg_bi,2));

theseResults = table;

for k=1:size(vectorEvents_min)
    
   thisRow = vectorEvents_min(k,:);
   
   %calc location w/in AP,LR axes, and report centroid
   
   %hard-coded ranges: x = 1-41 (left to right). y = 1-41 (post-ant; has
   %been flipped at this stage). 
   
%    theseX = trimmean(thisRow.PixelValues_min_x{:},2);
   theseLOCs =thisRow.PixelIdxList{:};

   thisX = thisRow.PixelValues_min_x{:};
%    thisX = mode(thisRow.PixelValues_min_x{:});
    %  L C R
   thisLR = [(thisX<xlen/3*1)', (thisX>=(xlen/3) & thisX<=(xlen/3*2))',(thisX>(xlen/3*2))']';           
    
   thisY = thisRow.PixelValues_min_y{:};
   
   % post mid anterior
   thisAP = [(thisY<ylen/3*1)', (thisY>=(ylen/3) & thisY<=(ylen/3*2))', (thisY>(ylen/3*2))']';
             
   LOC_LR(:,theseLOCs) = thisLR;
   LOC_AP(:,theseLOCs) = thisAP;
   
   ind1 = {[3,2,1],[5,6,4],[9,8,7]};
   
   thisLRAP=zeros(9,size(thisLR,2));
   for j=1:size(thisLR,2)
       thisLRAP(ind1{thisLR(:,j)}(thisAP(:,j)),j) = 1;
   end
   LOC_LRAP(:,theseLOCs)= thisLRAP;
end

LOC_LRAP_min = LOC_LRAP;
% 
% figure(1); imagesc(LOC_LR);
% figure; imagesc(LOC_AP);
% figure(26); imagesc(LOC_LRAP);


%% Generate events in vector space.  Positive vector

fct_showEEG_avg(seg_avg);
% imageAvgEEG(seg,seg_bi)
thresh_hi = 80;
thresh_lo = 50;

vectorEvents= fct_findVectorEvents(seg, Fs, thresh_hi, thresh_lo); %thresh_hi, thresh_lo, smoothWin, c_av)

LOC = fct_LOCfromEventsTbl({1},vectorEvents(strcmp(vectorEvents.type,'vector min'),:));
figure(37);imagesc(LOC)

% find x,y pos 
vectorEvents_min = vectorEvents(strcmp(vectorEvents.type,'vector max'),:);
xlen =41;
ylen =41;

LOC_LR = zeros(3,size(seg_bi,2));
LOC_AP = zeros(3,size(seg_bi,2));
LOC_LRAP = zeros(9, size(seg_bi,2));


theseResults = table;

for k=1:size(vectorEvents_min)
    
   thisRow = vectorEvents_min(k,:);
   
   %calc location w/in AP,LR axes, and report centroid
   
   %hard-coded ranges: x = 1-41 (left to right). y = 1-41 (post-ant; has
   %been flipped at this stage). 
   
%    theseX = trimmean(thisRow.PixelValues_min_x{:},2);
   theseLOCs =thisRow.PixelIdxList{:};

   thisX = thisRow.PixelValues_max_x{:};
%    thisX = mode(thisRow.PixelValues_min_x{:});
    %  L C R
   thisLR = [(thisX<xlen/3*1)', (thisX>=(xlen/3) & thisX<=(xlen/3*2))',(thisX>(xlen/3*2))']';           
    
   thisY = thisRow.PixelValues_max_y{:};
   
   % post mid anterior
   thisAP = [(thisY<ylen/3*1)', (thisY>=(ylen/3) & thisY<=(ylen/3*2))', (thisY>(ylen/3*2))']';
     
   LOC_LR(:,theseLOCs) = thisLR;
   LOC_AP(:,theseLOCs) = thisAP;
   
   ind1 = {[3,2,1],[5,6,4],[9,8,7]};
     
   thisLRAP=zeros(9,size(thisLR,2));
   for j=1:size(thisLR,2)
       thisLRAP(ind1{thisLR(:,j)}(thisAP(:,j)),j) = 1;
   end
   LOC_LRAP(:,theseLOCs)= thisLRAP;
end

LOC_LRAP_max = LOC_LRAP;


% figure(24); imagesc(LOC_LR);
% figure; imagesc(LOC_AP);
% figure(25); imagesc(LOC_LRAP);

LRAP_combo = LOC_LRAP_min+(LOC_LRAP_max*2);
figure(27); imagesc(LOC_LRAP_min+(LOC_LRAP_max*2));

figure(28); plot(sum(LRAP_combo(7:9,:),1)); axis tight;
figure; imagesc(sum(LRAP_combo(7:9,:),1));
figure; imagesc(sum(LRAP_combo(1:3,:),1));
figure; plot(sum(LRAP_combo(1:3,:),1)); axis tight;
%% Do the R^nk approach for naive event clustering in vector space

thresh_hi = 0;
thresh_lo = 0;

LOC = fct_LOCfromEventsTbl({1},vectorEvents(strcmp(vectorEvents.type,'vector min'),:));
figure(37);imagesc(LOC)

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
            thisBlock(j,thisStart:thisEnd) = thisData_byType.(theseData_toInclude{m});
            
        end
    end
    
    
%     eeg_hdv(j,:) = reshape(thisBlock,1,size(thisBlock,1)*size(thisBlock,2));
    eeg_hdv(j,:) = thisBlock;
    
end


%%
eeg_hdv = fct_findVectorEvents_Group( seg , Fs);

tic
Mdl = KDTreeSearcher(eeg_hdv);
toc

idx = knnsearch(Mdl,eeg_hdv(2354,:));

[cIdx,cD] = knnsearch(eeg_hdv,eeg_hdv(2462:2510,:),'K',100,'Distance','euclidean');

theseDistances = pdist2(eeg_hdv(1:80:end,:), eeg_hdv(1:80:end,:),'euclidean');
% figure;imagesc(reshape(theseDistances<0.1,2800,2800)); colormap(cool(7));
figure;imagesc(theseDistances); colormap(cool(7));

figure;plot(-1*theseDistances(1412,:));

[wcoeff,score,latent,tsquared,explained] = pca(eeg_hdv(1:40:end,:));
figure;
plot(score(:,1),score(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')


% convert to data in PC1/2
eeg_hdv(1130:1202,:).*repmat(wcoeff(:,1)',1202-1130+1,1)

%% Generate statistics for groups of discharges


%% Assign technical concerns related to bipolar events from average montage

% identify 'events' of high-noise, high or low amplitude voltage on each channel 

figure(13); fct_showEEG_avg(seg);
channelLabels_average_withspace = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'';'Fz';'Cz';'Pz';'';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};
channelLabels_average = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};

thresh_crazy = 30;
smoothWin = 10;
[LOC_artifact,artifactEvents]= fct_findArtifacts(seg, Fs, thresh_crazy, smoothWin, channelLabels_average);

fct_vizLOC_average(LOC_artifact);

% theseRect = fct_showEEG_avg_wEvents(seg, artifactEvents);
theseRect = fct_showEEG_avg_wEvents(seg, artifactEvents(artifactEvents.Area>Fs/4,:));

% add routine for 'mirroring' artifact in bipolar montage.

% add routine for 'high delta', 'high theta', 'high beta'
t_mag_delta = 15*5*10;
smoothWin_delta = Fs/5;

[LOC, deltaEvents, eeg_bp_smooth] = fct_findEvents(seg_bi, Fs, t_mag_delta, smoothWin_delta, 0, 4,channelLabels_bipolar);
fct_vizLOC_bipolar(LOC);
figure(14); showEEG(eeg_bp_smooth/10);

theseRect = fct_showEEG_wEvents(seg_bi, deltaEvents);

% delta from average 

seg_avg = seg-mean(seg,1);
[LOC, deltaEvents, eeg_bp_smooth] = fct_findEvents(seg_avg, Fs, t_mag_delta, ...
    smoothWin_delta, 0, 4,channelLabels_average);
fct_vizLOC_average(LOC);
figure(14); showEEG(eeg_bp_smooth/10);

figure(19);
theseRect = fct_showEEG_avg_wEvents(seg_avg, deltaEvents);

% theta from bipolar  -- for long stretches of theta

t_mag_theta = 15*2;
smoothWin_theta = Fs;
[LOC, thetaEvents, eeg_bp_smooth] = fct_findEvents(seg_bi, Fs, ...
    t_mag_theta, smoothWin_theta, 4, 9,channelLabels_bipolar);
% fct_vizLOC_bipolar(LOC);
figure; fct_showEEG_wEvents(seg_bi, thetaEvents);
figure(16); showEEG(eeg_bp_smooth);
% figure(17); fct_showEEG_avg(seg_avg);

% theta from average montage.

t_mag_theta = 15*5;
smoothWin_theta = Fs/10;
[LOC, thetaEvents, eeg_bp_smooth] = fct_findEvents(seg_avg, Fs, ...
    t_mag_theta, smoothWin_theta, 4, 8,channelLabels_average);
fct_vizLOC_average(LOC);
theseRect = fct_showEEG_avg_wEvents(seg_avg, thetaEvents);
figure(16); showEEG_avg(eeg_bp_smooth);
figure(17); fct_showEEG_avg(seg-mean(seg,1));

%alpha (8-12)

%beta (15-30)



% high frequency noise (30+Hz) from average.
t_mag_hf = 3;
smoothWin_hf = 4;
[LOC, noiseEvents, eeg_bp_smooth_noise] = fct_findEvents(seg_avg, Fs, ...
    t_mag_hf, smoothWin_hf, 30, 0,channelLabels_average);
fct_vizLOC_average(LOC);
% theseRect = fct_showEEG_avg_wEvents(seg_avg, noiseEvents);
figure(10); fct_showEEG_avg_wEvents(seg_avg, noiseEvents(noiseEvents.Area>0.1*Fs,:));
figure(17); fct_showEEG_avg(eeg_bp_smooth_noise);
figure(16); fct_showEEG_avg(seg_avg);

%*** routine to merge blocks or calculate a probability of noise
%contamination

% just high freq transients
t_mag_spike = 1;
smoothWin_spike = 3;
[LOC, spikeEvents, eeg_bp_smooth_spike] = fct_findEvents(seg_avg, Fs, ...
    t_mag_spike, smoothWin_spike, 14, 30,channelLabels_average);
fct_vizLOC_average(LOC);
% fct_showEEG_avg_wEvents(seg_avg, spikeEvents);
figure(10);fct_showEEG_avg_wEvents(seg_avg, spikeEvents((spikeEvents.MaxIntensity>50),:));
figure(11);fct_showEEG_avg_wEvents(seg_avg, spikeEvents((spikeEvents.MaxIntensity>50 & spikeEvents.Area<0.3*Fs),:));
figure(17); fct_showEEG_avg(eeg_bp_smooth_spike);
figure(16); fct_showEEG_avg(seg_avg);

figure(19); fct_showEEG_avg(eeg_bp_smooth_noise);
figure(18); fct_showEEG_avg((eeg_bp_smooth_spike-eeg_bp_smooth_noise));

% 4-20 on average montage
t_mag_discharge = 15*5;
smoothWin_discharge = Fs/10;
[LOC, dischargeEvents, eeg_bp_smooth_discharge] = fct_findEvents(seg_avg, Fs, ...
    t_mag_discharge, smoothWin_discharge, 4, 20,channelLabels_average);
fct_vizLOC_average(LOC);
% fct_showEEG_avg_wEvents(seg_avg, spikeEvents);
figure(21);fct_showEEG_avg_wEvents(seg_avg, dischargeEvents);


%% Using a 'pseudo-wavelet' based approach for pattern matching. 

thisEvent = pdEvents(3,:);
thisEvent_pixels = thisEvent.PixelValues_raw{:};
thisEvent_locs = thisEvent.PixelIdxList{:};
thisEvent_channel = thisEvent.channelNum;

thisEEGChannel = seg_bi(thisEvent_channel,:);

w = conv(thisEEGChannel, thisEvent_pixels,'same');
figure;plot(w); axis tight;
% figure;plot(thisEEGChannel)
figure; plot(thisEvent_pixels);



% merge event pixels 
theseEvents_locs = unique(cat(1,pdEvents.PixelIdxList{:}));
theseEvents_chans = unique(cat(1,pdEvents.channelNum(:)));

thisPatternBlock = seg_bi(theseEvents_chans,theseEvents_locs);
theseEEGChannels = seg_bi(theseEvents_chans,:);

w = conv2(theseEEGChannels, thisPatternBlock,'same');
figure;plot(w'); axis tight;
figure;plot(prod(w,1)); axis tight;
% figure;plot(thisEEGChannel)
figure; plot(thisEvent_pixels);

figure(59);imagesc(w);

%normalize by rms of eeg x block?
thisDenom = rms(theseEEGChannels,2).*rms(thisPatternBlock,2);

wvlt_events = prod(w./repmat(thisDenom,1,size(w,2)),1).^2;
wvlt_events_smooth = simplema(wvlt_events,Fs/20,2); 

figure; plot(wvlt_events);
figure(50); plot(wvlt_events_smooth); axis tight;

% find other events based on peaks algorithm
[theseOtherEvents_val, theseOtherEvents_locs]= findpeaks(wvlt_events_smooth);
figure(50); hold on; scatter(theseOtherEvents_locs, theseOtherEvents_val);

figure(51); imagesc(wvlt_events_smooth>0.2e11);


%% wavelet-based approach with all channels


% merge event pixels 
theseEvents_locs = unique(cat(1,pdEvents.PixelIdxList{:}));
theseEvents_chans = unique(cat(1,pdEvents.channelNum(:)));

thisPatternBlock = seg_bi(theseEvents_chans,theseEvents_locs);
theseEEGChannels = seg_bi(theseEvents_chans,:);

w = conv2(theseEEGChannels, thisPatternBlock);
figure;plot(w'); axis tight;
figure;plot(prod(w,1)); axis tight;
figure;imagesc(w);
figure;plot(sum(w,1)); axis tight;
% figure;plot(thisEEGChannel)
figure; plot(thisEvent_pixels);


%% iterative 1D convolution for comparison

w_i = zeros(size(theseEEGChannels));
for j=1:size(theseEEGChannels,1)
   
%     thisEEGChannel = theseEEGChannels(j,:);
%     thisPatternBlock_row = thisPatternBlock(j,:);
%     
    thisEEGChannel = theseEEGChannels(j,:)./rms(theseEEGChannels(j,:),2);
    thisPatternBlock_row = thisPatternBlock(j,:)./rms(thisPatternBlock(j,:),2);
    
    w_i(j,:) = conv(thisEEGChannel, thisPatternBlock_row,'same');
    
end

figure(55); imagesc(w_i);
w=w_i;
%%

%normalize by rms of eeg x block?
% thisDenom = max(theseEEGChannels,[],2).*max(thisPatternBlock,[],2);
thisDenom = max(w.^2,[],2);

% wvlt_events = prod(w./repmat(thisDenom,1,size(w,2)),1).^2;
wvlt_events = prod((w.^2)./repmat(thisDenom,1,size(w,2)),1);
wvlt_events_smooth = simplema(wvlt_events,Fs/20,2); 

figure; plot(wvlt_events); axis tight;
figure(50); plot(wvlt_events_smooth); axis tight;
figure(54); imagesc((w.^2)./repmat(thisDenom,1,size(w,2)));


% find other events based on peaks algorithm
[theseOtherEvents_val, theseOtherEvents_locs]= findpeaks(wvlt_events_smooth);
figure(52); hold on; scatter(theseOtherEvents_locs, theseOtherEvents_val);

figure(51); imagesc(wvlt_events_smooth>0.2e11);

%% R^nk vector approach for finding similar events

% criteria = {'(e.MaxIntensity_raw-e.MinIntensity_raw)>25'};
criteria_p = {"cellfun(@(x) abs(diff(x))>7, e.PixelValues_raw, 'UniformOutput', false)"}
criteria_cellfun = join(cat(2,"cat(1,cellfun(@(x) any(x), ", criteria_p{:}, " , 'UniformOutput', false))"));

% pdEvents = selectEvents(c_bi, bipolarEvents, criteria);
pdEvents = selectEvents(c_bi, bipolarEvents, criteria_cellfun);
figure(44); fct_showEEG_wEvents(seg_bi, pdEvents);


thisEvent = pdEvents(13,:);
figure;plot(abs(diff(thisEvent.PixelValues_raw{:})))

theseEvents_locs = unique(cat(1,thisEvent.PixelIdxList{:}));
theseEvents_chans = unique(cat(1,thisEvent.channelNum(:)));

% thisPatternBlock = seg_bi(theseEvents_chans,theseEvents_locs);
% theseEEGChannels = seg_bi(theseEvents_chans,:);

thisPatternBlock = seg_bi(channels_R,theseEvents_locs);
theseEEGChannels = seg_bi(channels_R,:);
% 
% thisPatternBlock = seg_bi(channels_L,theseEvents_locs);
% theseEEGChannels = seg_bi(channels_L,:);

eeg_hdv = zeros(size(theseEEGChannels,2),size(theseEEGChannels,1)*size(thisPatternBlock,2));
blockSize = size(thisPatternBlock,2);

for j=1:size(theseEEGChannels,2)-blockSize
    
    thisBlock = theseEEGChannels(:,j:j+blockSize-1);
    eeg_hdv(j,:) = reshape(thisBlock,1,size(thisBlock,1)*size(thisBlock,2));
    
end

theseDistances = pdist2(eeg_hdv(1:end,:), eeg_hdv(1:end,:),'correlation');
figure;imagesc(theseDistances); colormap(cool(7));

thisSeed = thisEvent.PixelIdxList{:};
thisSeed= thisSeed(1);

figure;imagesc(theseDistances(theseEvents_locs,:)); colormap(cool(7));
figure(90); plot(-1*theseDistances(thisSeed,:)); axis tight;
figure(87); imagesc(theseDistances(thisSeed,:)<0.7); axis tight;

%% Use LOC from R^nk correlation to select events

groupLOC = theseDistances(thisSeed,:)<0.7;
groupLOC_st = regionprops(groupLOC,'PixelIdxList');
groupLOC_locs = cat(1,groupLOC_st(:).PixelIdxList);

bw_groupLOC_d = imdilate(groupLOC,strel('line',5,180));  %% dilate by 1 pixel 
figure(92); ax(1)=subplot(211); imagesc(groupLOC); ax(2)=subplot(212); imagesc(bw_groupLOC_d);
linkaxes(ax);


groupLOC_locs_cell = repmat({groupLOC_locs},size(bipolarEvents,1),1);

% criteria = {'ismember(inVar,e.PixelIdxList)'};
% criteria = {"cat(2,(cellfun(@(x,y) ismember(x,y), inVar, e.PixelIdxList, 'UniformOutput',false)))"};

% criteria = {"ismember(e.channelNum,[1,2,3,4,9,10,11,12])",
%     "cat(1,cellfun(@(x) any(x), (cellfun(@(x,y) ismember(x,y), inVar, e.PixelIdxList, 'UniformOutput',false)), 'UniformOutput', false))"};
% 
% criteria = {"ismember(e.channelNum,[5,6,7,8,13,14,15,16])",
%     "cat(1,cellfun(@(x) any(x), (cellfun(@(x,y) ismember(x,y), inVar, e.PixelIdxList, 'UniformOutput',false)), 'UniformOutput', false))"};


criteria = {"ismember(e.channelNum,inVar{2})",
    "cat(1,cellfun(@(x) any(x), (cellfun(@(x,y) ismember(x,y), inVar{1}, e.PixelIdxList, 'UniformOutput',false)), 'UniformOutput', false))"};

groupEvents = selectEvents(c_bi, bipolarEvents, criteria,{groupLOC_locs_cell,channels_R});

figure; fct_showEEG_wEvents(seg_bi,groupEvents);
figure(87); imagesc(bw_groupLOC_d); axis tight;
LOC = fct_LOCfromEventsTbl({1},groupEvents);
figure(37);imagesc(LOC)


% expand vertically again.

groupEvents_st = regionprops(any(LOC),'PixelIdxList');
groupEvents_locs = unique(cat(1,groupEvents_st(:).PixelIdxList));
groupEvents_locs_cell = repmat({groupEvents_locs},size(bipolarEvents,1),1);

% 
% criteria = {"ismember(e.channelNum,[1,2,3,4,9,10,11,12])",
%     "cat(1,cellfun(@(x) any(x), (cellfun(@(x,y) ismember(x,y), inVar, e.PixelIdxList, 'UniformOutput',false)), 'UniformOutput', false))"};
% % criteria = {"ismember(e.channelNum,[5,6,7,8,13,14,15,16])",
% %     "cat(1,cellfun(@(x) any(x), (cellfun(@(x,y) ismember(x,y), inVar, e.PixelIdxList, 'UniformOutput',false)), 'UniformOutput', false))"};

groupEvents_2 = selectEvents(c_bi, bipolarEvents, criteria,{groupEvents_locs_cell,channels_R});

figure; fct_showEEG_wEvents(seg_bi,groupEvents_2);
% figure(88); imagesc(bw_groupLOC_d); axis tight;
LOC = fct_LOCfromEventsTbl({1},groupEvents_2);
figure(39);imagesc(LOC)
figure; plot(sum(LOC)); axis tight;
fct_vizEEG_bipolar(eeg_bp_smooth);

%% Generate vector events
% -- vectorEvents obj



%% Visualize vector events

%% Classify vector events as "type" of periodic discharges

%% Cross-check bipolar Events by vector events 
% -- inherit classification, assign to VE

%% makeGroups of vector events
% -- use a k-means or pca approach

%% Assess periodicity of vectors group

%% makeGroups of bipolar events


%% compute correlation with vector Events

%% Assess periodicity of bipolar groups

theseEvents = bipolarEvents;

theseEvents_n = fct_assignPeriodicity(channelLabels_bipolar,theseEvents);

% compare results with eeg_bp vs the LOC. 

%%

% junk
channels = unique(theseEvents.channelNum);

for i=1:length(channels)

[autocor,lags] = xcorr(LOC',3*Fs,'coeff');

end



[autocor,lags] = xcorr(eeg_bp_smooth_discharge',3*Fs,'coeff');
toc
figure;
plot(lags/Fs,autocor(:,[1:19:361]));
ylim([0 1]);

figure;
plot(lags/Fs,sum(autocor(:,[1:19:361]),2))



%% Classify epochs
% -- based on superGroup classification of bipolar events
% -- based on group classification of vector events