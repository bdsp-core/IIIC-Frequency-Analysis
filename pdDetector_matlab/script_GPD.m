%% script for running GPD analyzer. 

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


%% Using a per-channel percentile-based auto threshold with lower threshold based ...
% on slope == 0. 

n=6;   % 15 = low amplitude, with ?BiPDs?
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

thresh_hi = 90;
thresh_lo = 90;
smoothWin = 30;  % typically 120msec = 24 samples
boostWin = Fs/10;
meanWin = 80;

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_bi, Fs, ...
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
%% Group stats 
thresh_hi = 90;
thresh_lo = 75;
smoothWin = Fs/10;
boostWin = Fs/10;
meanWin = 80;


n=6;   % 15 = low amplitude, with ?BiPDs?; % 20, great LPD_R>L. 
%  19 =low amplitude on R with small PLDs, 
% 22= some artifact
%28 = artifact
%39 = noisy on L, R has low amplitude but recurrent blips 
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

% 
% [LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_bi, Fs, ...
%     thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_bi);

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_bi);

[LOC, avgEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_avg, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_av);
% 
% [LOC, avgEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_avg, Fs, ...
%     thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_av);

% figure(96); fct_showEEG_wEvents(seg_bi, bipolarEvents);
% figure(97); fct_showEEG_avg_wEvents(seg_avg, avgEvents);
% figure(99); fct_showEEG_avg(seg_avg);
% figure(100); showEEG(seg_bi);
% figure(102); fct_showEEG_avg(eeg_bp_smooth/5);
% fct_vizEEG_average(eeg_bp_smooth)
%  fct_vizLOC_average(LOC)

% criteria = {'(e.MaxIntensity)>50','e.Area < 0.30*Fs'};
criteria = {'(e.MaxIntensity)>18','(e.MaxIntensity_raw)-(e.MinIntensity_raw)>25',...
    '(e.MaxIntensity_raw<100)','(e.MinIntensity_raw>-100)'};

% criteria = {'(e.MaxIntensity)>18','(e.MaxIntensity_raw)-(e.MinIntensity_raw)>25'};

pdEvents_av = selectEvents(c_av, avgEvents, criteria);
figure(44); fct_showEEG_avg_wEvents(seg_avg, pdEvents_av);

pd_av_LOC = fct_LOCfromEventsTbl(c_av, pdEvents_av);
fct_vizLOC_average(pd_av_LOC);

pdEvents_bi = selectEvents(c_bi, bipolarEvents, criteria);
figure(45); fct_showEEG_wEvents(seg_bi, pdEvents_bi);

pd_bi_LOC = fct_LOCfromEventsTbl(c_av, pdEvents_bi);
fct_vizLOC_bipolar(pd_bi_LOC);

% groupEvents = fct_findGroupEvents(avgEvents, seg_avg, eeg_bp_smooth, c_av);

groupEvents_av = fct_findGroupEvents(pdEvents_av, seg_avg, eeg_bp_smooth, c_av);
groupEvents_bi = fct_findGroupEvents(pdEvents_bi, seg_bi, eeg_bp_smooth, c_bi);
%% Results for test run

% Inter-event interval estimate
IDI = groupEvents_av(ismember(groupEvents_av.channel,'Global_kmeans'),:);
IDI.eventRate
% IDI.thisIDI_Hz
IDI.thisIDI_Hz_med
IDI.thisIDI_xcorr_BP_Hz

% Channels involved
theseChannels = groupEvents_av(1:19,:);
theseChannels_involved = theseChannels(theseChannels.events>0,:).channel;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find events in 3D

theseTS = {seg, eeg_bp_smooth, pd_av_LOC};

labels = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};
placeOnHead = [[1 2];[2 2];[3 2];[4 2];[2 1];[3 1];[4 1];[5 2]; [2 3];[3 3];[4 3];[1 4];[2 4];[3 4];[4 4];[2 5];[3 5];[4 5];[5 4]]; 
% (y,x) -- 1st dimension is AP, 2nd dimension is L-R.
nanOnHead = [[1 1]; [1 3]; [1 5];[5 1 ];[5 3];[5 5]]; 
labelMat_head = zeros(5,5,1);

theseImgAvg = cell(2,1);
for k=1:3
    
    thisTS = theseTS{k};

 data_head = zeros(5,5,size(thisTS,2));
 data_head(nanOnHead(:,1),nanOnHead(:,2),:) = nan;
 
 labels_head = cell(5,5,1);
 
 Fs = 200;
%  seg_bp = eegfilt(seg,Fs,4,0);   % high pass above 4
 
 for i = 1:length(labels)
    
    thisIndx = placeOnHead(i,:);
    labels_head(thisIndx(1),thisIndx(2)) = labels(i);
    data_head(thisIndx(1),thisIndx(2),:) = thisTS(i,:);
    
%     data_head(thisIndx(1),thisIndx(2),:) = seg_bp(i,:);
    
    
    labelMat_head(thisIndx(1),thisIndx(2)) = i;
    
 end
 
 % could start with seg_avg instead of seg and avoid this step 
 meanPerUnitTime = mean(reshape(data_head,5*5,size(data_head,3)),'omitnan');
 data_head_n = data_head-repmat(reshape(meanPerUnitTime,1,1,size(meanPerUnitTime,2)),5,5);
 
 img_avg = interpolateImageNan(data_head_n);
 
labelMat_head_io = interpolateImage(labelMat_head);
 
theseImgAvg{k} = img_avg;
end

img_avg_raw = theseImgAvg{1};
img_avg_bp = theseImgAvg{2};
LOC_head = theseImgAvg{3};

 trange = [1:10:2800];
 
 t=310;
 tr = t-Fs/2:t+Fs/2;
 figure(11);imagesc(mean(img_avg(:,:,tr),3),[1 50]);figure(10);imagesc(mean(data_head_n(:,:,tr),3),[1 50]);
%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % regionprops in 3D on interpolated filtered data
 
%   lpds_hi = regionprops(img_avg,'PixelIdxList');
%   CC= bwconncomp(data_head_n>35);
%   img_avg_bp = img_avg;
%   img_avg_raw = img_avg; 
  
  img_avg_cens = img_avg_bp;  % 25x25xt
  corners = {{1:5; 1:5},{21:25;1:5},{21:25,21:25},{1:5,21:25},{1:5,11:15},{21:25,11:15}};
  
  for j=1:length(corners)
  img_avg_cens(corners{j}{1}, corners{j}{2}, 1:2800) = NaN;
  end
  
%   d= img_avg_cens(:,:,10);

  img_avg_cens = img_avg_bp; 
  
  t_bpEEG = 20;
  t_bpEEG_prt = prctile(eeg_bp_smooth(:),75);
  
  t_area = 40;
  t_rawEEG = 40;
  t_rawEEG_range = 50;
  
%   CC= bwconncomp(img_avg_cens>t_bpEEG_prt);
  CC= bwconncomp(img_avg_cens>min(t_bpEEG,t_bpEEG_prt));
%   CC= bwconncomp(img_avg_cens>t_bpEEG_prt);
    
  CC=bwconncomp(LOC_head>0);
%   figure(10); imagesc(mean(LOC_head(:,:,200:300),3),[0 1])
%   
%   events = regionprops('table',CC,img_avg_cens,'PixelIdxList','Centroid','Area','BoundingBox',...
%       'MaxIntensity');
%   

 CC_l = labelmatrix(CC);
 figure; imagesc(CC_l(:,:,200));
  events = regionprops('table',CC,img_avg_raw,'PixelIdxList','Centroid',...
    'Area','BoundingBox','MaxIntensity','MinIntensity');
  
  events.Area_xy = events.BoundingBox(:,4).*events.BoundingBox(:,5);
  events.Area_xy_n = events.Area_xy./25;
  events.IntensityRange = events.MaxIntensity-events.MinIntensity;
  
%   indx = [events.Area_xy>t_area & events.MaxIntensity>t_rawEEG];
  indx = [events.Area_xy>t_area & events.IntensityRange>t_rawEEG_range];
  

  events_t = events.Centroid(indx,3);
  events_s = events.Centroid(indx,1:2);
%   events_l = sub2ind([21,21],events_s(:,1),events_s(:,2));
  
theseEvents = events(indx,:);
%   figure;scatter(events_t,events_l);
  
e_left = theseEvents(theseEvents.Centroid(:,1)<10,:);
e_right= theseEvents(theseEvents.Centroid(:,1)>15,:);
e_center = theseEvents(theseEvents.Centroid(:,1)>10 & theseEvents.Centroid(:,1)<15,:);


  figure(12); histogram2(events_s(:,1),events_s(:,2),1:25,1:25);
  
  figure(11);imagesc(mean(img_avg_cens(:,:,tr),3),[1 50]);figure(10);imagesc(mean(data_head_n(:,:,tr),3),[1 50]);
  figure(13); scatter(e_right.Centroid(:,3), ones(size(e_right,1),1)); xlim([0 2800]);
  figure(14); scatter(e_left.Centroid(:,3), ones(size(e_left,1),1));xlim([0 2800]);
  figure(15); scatter(e_center.Centroid(:,3), ones(size(e_center,1),1));xlim([0 2800]);
  
  figure(13); scatter(e_right.Centroid(:,3), e_right.Area_xy_n); xlim([0 2800]); ylim([0 25]); 
  figure(14); scatter(e_left.Centroid(:,3), e_left.Area_xy_n);xlim([0 2800]); ylim([0 25]);
  figure(15); scatter(e_center.Centroid(:,3), e_center.Area_xy_n);xlim([0 2800]); ylim([0 25]);
  
  loc_center = zeros(2800,1);
  start = floor(e_center.BoundingBox(:,3));
  stop = start+e_center.BoundingBox(:,4);
  
  for k=1:length(start)
  loc_center(start(k):stop(k))=1;
  end
%   loc_center([round(e_center(.BoundingBox(3)),e_center.BoundingBox(6)])=1;
  figure; imagesc(loc_center');
  
  
  %% KNN clustering
  
  % just in spatial domain -- use centroid (x,y), boundingbox (x,y), 
 test = kmeans(cat(2,theseEvents.Centroid(:,1:2),theseEvents.Area_xy_n),3);
 
 for j=1:3
     e_g = theseEvents(test==j,:);
     loc_center = zeros(2800,1);
     start = floor(e_g.BoundingBox(:,3));
     stop = start+e_g.BoundingBox(:,4);
     
     for k=1:length(start)
         loc_center(start(k):stop(k))=1;
     end
     %   loc_center([round(e_center(.BoundingBox(3)),e_center.BoundingBox(6)])=1;
     figure; imagesc(loc_center');
 end
     
 %% For rhythmic delta
 
thresh_hi = 85;
thresh_lo = 50;
smoothWin = Fs;
boostWin = Fs;
meanWin = 80;

n=4;   

[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto4(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 0,4, boostWin, meanWin, channelLabels_bipolar);

[LOC, avgEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto4(seg_avg, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 0,4, boostWin, meanWin, c_av);

criteria = {'(e.MaxIntensity)>18','(e.MaxIntensity_raw)-(e.MinIntensity_raw)>25'};
 
pdEvents_av = selectEvents(c_av, avgEvents, criteria);
figure(44); fct_showEEG_avg_wEvents(seg_avg, pdEvents_av);

pd_LOC = fct_LOCfromEventsTbl(c_av, pdEvents_av);
fct_vizLOC_average(pd_LOC);

pdEvents_bi = selectEvents(c_bi, bipolarEvents, criteria);
figure(45); fct_showEEG_wEvents(seg_bi, pdEvents_bi);

pd_LOC = fct_LOCfromEventsTbl(c_av, pdEvents_bi);
fct_vizLOC_bipolar(pd_LOC);
  
  
  