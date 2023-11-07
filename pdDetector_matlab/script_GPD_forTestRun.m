%% script for running GPD analyzer. 

% Take file, seg. 
% dataDirec = 'C:\Users\Chris McGraw\Dropbox (Partners HealthCare)\01 Data\2018\FrequencyOfDischarges_Chris_JJ\code\data\';
% subDirec = [dataDirec '\LPD\'];
dataDirec = '/data/01 Code/MATLAB/05 PD/02 Code for PD detector/data';
% dataDirec = 'C:\Users\Chris McGraw\Dropbox (Partners HealthCare)\01 Data\2018\FrequencyOfDischarges_Chris_JJ\code\';
% subDirec = [dataDirec 'from JJ\'];
subDirec = [dataDirec '/GPD/'];
% subDirec = [dataDirec '\LRDA\'];
% subDirec = [dataDirec '\GRDA\']; 
% subDirec = [dataDirec '\Sz\']; 

files = dir([subDirec '*.mat']);

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

%% For GPD/LPD detection.
%Using a per-channel percentile-based auto threshold with lower threshold based ...
% on slope == 0. 

n=36;   
[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

thresh_hi = 90;
thresh_lo = 75;
smoothWin = Fs/10;
boostWin = Fs/10;
meanWin = 80;

totalSamples = size(seg,2);

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_bi);
% 
% [LOC, avgEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_avg, Fs, ...
%     thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_av);

% Visualization functions
% figure(96); fct_showEEG_wEvents(seg_bi, bipolarEvents);
% figure(97); fct_showEEG_avg_wEvents(seg_avg, avgEvents);
% figure(99); fct_showEEG_avg(seg_avg);
% figure(100); showEEG(seg_bi);
% figure(102); fct_showEEG_avg(eeg_bp_smooth/5);
% fct_vizEEG_average(eeg_bp_smooth)
%  fct_vizLOC_average(LOC)

fct_vizEEG_bipolar(eeg_bp_smooth)

criteria = {'(e.MaxIntensity)>18','(e.MaxIntensity_raw)-(e.MinIntensity_raw)>25',...
    '(e.MaxIntensity_raw<100)','(e.MinIntensity_raw>-100)'};

% 
% pdEvents_av = selectEvents(c_av, avgEvents, criteria);
% figure(44); fct_showEEG_avg_wEvents(seg_avg, pdEvents_av);
% 
% pd_av_LOC = fct_LOCfromEventsTbl(c_av, pdEvents_av, totalSamples);
% fct_vizLOC_average(pd_av_LOC);

pdEvents_bi = selectEvents(c_bi, bipolarEvents, criteria);
figure(45); fct_showEEG_wEvents(seg_bi, pdEvents_bi);

pd_bi_LOC = fct_LOCfromEventsTbl(c_bi, pdEvents_bi, totalSamples);
fct_vizLOC_bipolar(pd_bi_LOC);

% groupEvents_av = fct_findGroupEvents(pdEvents_av, seg_avg, eeg_bp_smooth, c_av);
groupEvents_bi = fct_findGroupEvents_bi(pdEvents_bi, seg_bi, eeg_bp_smooth, c_bi);

% pow_L = groupEvents_bi(23,:).power_discharges
% pow_R = groupEvents_bi(24,:).power_discharges
% pow_G = groupEvents_bi(26,:).power_discharges
% 
% pow_Lr = groupEvents_bi(23,:).power_discharges_rel
% pow_Rr = groupEvents_bi(24,:).power_discharges_rel
% pow_Gr = groupEvents_bi(26,:).power_discharges_rel
% 
% log(pow_L/pow_R)
% log(pow_Lr/pow_Rr)

[freq, channels, interpretation] = fct_PD_detect(seg);

channels
freq
1/freq


%% Results for test run

% Inter-event interval estimate
IDI = groupEvents_av(ismember(groupEvents_av.channel,'Global_kmeans'),:);
IDI.eventRate
% IDI.thisIDI_Hz
IDI.thisIDI_Hz_med
% IDI.thisIDI_xcorr_BP_Hz

% Channels involved
theseChannels_involved = unique(pdEvents_av.channel);

 %% For rhythmic delta
 
thresh_hi = 85;
thresh_lo = 50;
smoothWin = Fs;
boostWin = Fs;
meanWin = 80;

n=1;   

[seg, seg_bi] = loadFile(files,n);
seg_avg = seg-mean(seg,1);

totalSamples = size(seg,2);

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto4(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 0,4, boostWin, meanWin, channelLabels_bipolar);
% 
% [LOC, avgEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto4(seg_avg, Fs, ...
%     thresh_hi, thresh_lo, smoothWin, 0,4, boostWin, meanWin, c_av);

criteria = {'(e.MaxIntensity)>18','(e.MaxIntensity_raw)-(e.MinIntensity_raw)>25'};
%  
% pdEvents_av = selectEvents(c_av, avgEvents, criteria);
% % figure(44); fct_showEEG_avg_wEvents(seg_avg, pdEvents_av);
% 
% pd_LOC = fct_LOCfromEventsTbl(c_av, pdEvents_av, totalSamples);
% fct_vizLOC_average(pd_LOC);

pdEvents_bi = selectEvents(c_bi, bipolarEvents, criteria);
figure(45); fct_showEEG_wEvents(seg_bi, pdEvents_bi);

pd_LOC = fct_LOCfromEventsTbl(c_bi, pdEvents_bi, totalSamples);
fct_vizLOC_bipolar(pd_LOC);
  
% Channels involved
theseChannels_involved = unique(pdEvents_av.channel);


%% Test of wrapper function. 
%subDirec = [dataDirec '/LPD/'];
%subDirec = [dataDirec '/GPD/'];
%subDirec = [dataDirec '/LRDA/'];
%subDirec = [dataDirec '/GRDA/']; 
subDirec = [dataDirec '/Sz/']; 

files = dir([subDirec '*.mat']);

n=18;   
[seg, seg_bi] = loadFile(files,n);
[freq, channels, interpretation] = fct_PD_detect(seg);


freq
% channels
interpretation

% RDA looks like it would benefit from an alternative approach -- if frequency is all that is required, 
% just take the band pass power from 0.5-4Hz and find the peak. 
% for spatial extent, would need to continue to use my method though...

%% 10-7-22

theseDir = {'/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Hard/', ...
    '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Medium/',...
    '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Easy/'};

for j=1:length(theseDir)
    
    % parentDir = '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Medium/';
    parentDir = theseDir{j};
    
    % subDirec = ['/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Hard/'];
    % subDirec = ['/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Easy/'];
    
    addpath(parentDir);
    
    formatOut = 'yy-mm-dd';
    thisDate = datestr(now,formatOut);
    
    % directory structure
    % >thisDir            a directory containing AVI files
    % >> avi files        to be processed
    % >> analysis_<date>    subDir to place CSV/XLS output from analysis
    % >>> mat               subDir to place MAT files
    % >>> images            subDir to place figures/AVI output
    
    thisAnalysisDir = [parentDir 'analysis_' thisDate '/'];
    thisImageDir = [thisAnalysisDir 'images/'];
    thisMATDir = [thisAnalysisDir 'mat/'];
    make_directories({thisAnalysisDir, thisImageDir, thisMATDir});
    addpath(thisAnalysisDir);
    addpath(thisImageDir);
    
    files = dir([parentDir '*.mat']);
    
    resultTable = table;
    
    tic
    for n=1:length(files)
        
        thisResult=table;
        S =load([files(n).folder '/' files(n).name]);
        % seg= S.data(1:19,:);  % remove EKG
        
        P_model = S.P_model;
        Y_human = S.Y_human;
        
        Y_human_numVotes = sum(Y_human{1,:});
        
        Y_human_p = array2table(Y_human{1,:}./sum(Y_human{1,:}),'VariableNames',Y_human.Properties.VariableNames);  % convert to percentage
        
        % determine which detection function to run based on human rating.
        Y_human_s = Y_human_p(:,Y_human_p{1,:}>=1/3); % all types >1/3
        theseTypes = Y_human_s.Properties.VariableNames;
        theseTypes_p = Y_human_s{1,:};
        
        thisResult.n = n;
        thisResult.file = {files(n).name};
        thisResult.y_human_types = {theseTypes};
        thisResult.y_human_types_prct = {theseTypes_p};
        thisResult.statsObj = {NaN};
        
        Y_human_s_m = Y_human_s(:,ismax(Y_human_s{1,:})); % top type
        thisType = Y_human_s_m.Properties.VariableNames;
        
        if ~isempty(thisType)
            
            switch thisType{:}
                case {'GRDA','LRDA'}
                    disp('grda,lrda');
                    resultTable = cat(1,resultTable,thisResult);
                    continue
                case {'LPD','GPD'}
                    disp('lpd, gpd');
                    
                    %      S =load([files(n).folder '/' files(n).name]);
                    seg= S.data(1:19,:);  % remove EKG
                    
                    [freq, channels, interpretation, data_obj] = fct_PD_detect(seg);
                    
                case 'Seizure'
                    disp('seizure');
                    seg= S.data(1:19,:);  % remove EKG
                    
                    [freq, channels, interpretation, data_obj] = fct_PD_detect(seg);
                    
                case 'Other'
                    disp('other');
                    resultTable = cat(1,resultTable,thisResult);
                    continue
            end
            
            %% Generate figure
            if isstruct(data_obj)
                data_obj.num = n;
                data_obj.filename = files(n).name;
                data_obj.theseTypes = theseTypes;
                data_obj.theseTypes_prctVotes = theseTypes_p;
                data_obj.Y_human = Y_human;
                data_obj.P_model = P_model;

                [f1, thisStatsObj] = fct_showEEG_wEventsAndStats_dataObj(data_obj);
%                 set(f1, 'Visible','on');
                % thisStatsObj contains the stats re: evolution 

                print(f1, '-dpng','-zbuffer','-r600', [thisImageDir files(n).name '.png']);
                close(f1);
                
                save([thisMATDir files(n).name '.mat'],'data_obj','thisStatsObj');
                thisResult.statsObj = {thisStatsObj};
                
            end
            
        end
        %%
        
        resultTable = cat(1,resultTable,thisResult);
        
        
    end
    toc
    save([thisMATDir  'resultTable.mat'],'resultTable');

end

%% Find and move files based on list 

dirForFilenames = '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/Images_Files2Check_2022Oct20/Hard/';
dirWithPNG = '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Hard/analysis_22-10-22/images/';
targetDirToMove = '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Hard/analysis_22-10-22/images/FilesToCheck/';

make_directories({targetDirToMove});

filesForNames = dir([dirForFilenames '*.png']);

filesToMove = dir([dirWithPNG '*.png']);
filesToMove_c = struct2cell(filesToMove)';

filesToMove_ct = cellfun(@(x) x(1:end-8),filesToMove_c(:,1),'UniformOutput',false);

for j=1:length(filesForNames)
    
    [a,b] = ismember(filesForNames(j).name(1:end-4), filesToMove_ct(:,1));
    
    if a
        % found -- move to targetDir
        thisTargetFile = [filesToMove(b).folder '/' filesToMove(b).name];
        disp(['Moving ' thisTargetFile]);
        movefile(thisTargetFile,targetDirToMove);
    end
    
end

%%

dirWithPNG = '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Medium/analysis_22-10-22/images/';
zipFileName= 'images_medium_221013.zip';

system(['zip -r ' zipFileName ' ''' dirWithPNG '''']);
% system(['zip -r images_hard_221023.zip ''' dirWithPNG '''']);

