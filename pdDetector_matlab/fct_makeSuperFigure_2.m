function [ output_args ] = fct_makeSuperFigure_2( seg, file)
%fct_makeSuperFigure 
%INPUT
%   ****
%OUTPUT
% ***

f1= figure('Position',[0 0 2400 2400]);
% [ha, pos] = tight_subplot(8,4,[.01 .03],[.1 .01],[.01 .01]) 
% for ii = 1:16; axes(ha(ii)); plot(randn(10,ii)); end 
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

channels_L = [1:4,9:12];
channels_R = [5:8, 13:16];
channels_C = [17,18];

channelLabels_bipolar = {'Fp1-F7';'F7-T3';'T3-T5';'T5-O1'; 'Fp2-F8';'F8-T4';'T4-T6';'T6-O2'; 'Fp1-F3';'F3-C3';'C3-P3';'P3-O1'; 'Fp2-F4';'F4-C4';'C4-P4';'P4-O2'; 'Fz-Cz';'Cz-Pz'};
channelLabels_bipolar_withspace = {'Fp1-F7';'F7-T3';'T3-T5';'T5-O1'; '';'Fp2-F8';'F8-T4';'T4-T6';'T6-O2'; '';'Fp1-F3';'F3-C3';'C3-P3';'P3-O1'; '';'Fp2-F4';'F4-C4';'C4-P4';'P4-O2'; '';'Fz-Cz';'Cz-Pz'};
channelLabels_average= {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};
channelLabels_average_withspace = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'';'Fz';'Cz';'Pz';'';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};

c_bi = channelLabels_bipolar;
c_av = channelLabels_average;

ax_v = 4;
ax_h = 7;
%% Find bipolar events
thresh_hi = 90;
thresh_lo = 85;
smoothWin = Fs/10;
boostWin = Fs/10;
meanWin = 80;

[LOC, bipolarEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_bi, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_bi);


%% Row 1

% raw bipolar eeg
subplot(ax_v, ax_h, [1 2]); 
showEEG(seg_bi);

% bipolar eeg w events
subplot(ax_v, ax_h, [3 4]);
fct_showEEG_wEvents(seg_bi, bipolarEvents);

% raw avg eeg 
subplot(ax_v, ax_h, [5 6]);
seg_avg = seg-mean(seg,1);
fct_showEEG_avg(seg_avg)


%% Event related figures   

% events as image
subplot(ax_v, ax_h, [10:11]);
im_LOC = LOC;
im_LOC(channels_L,:)=im_LOC(channels_L,:)*3;
im_LOC([17,18],:)=im_LOC([17,18],:)*2;
imagesc(im_LOC([1:4,9:12,5:8,13:18],:));
title('LOC by channel, left v right v center')

% summed events as image, LvRvC
subplot(ax_v, ax_h, [17 18]);
im_LOC2 = [sum(LOC([1:4,9:12],:),1); sum(LOC([5:8,13:16],:),1); ...
    sum(LOC([17,18],:),1)];
imagesc(im_LOC2);
title('summed LOC, L v R v C');

% summed events L vs R
subplot(ax_v, ax_h, [24 25]);
plot(1:size(LOC,2),im_LOC2(1,:)); axis([0 size(LOC,2) 0 8]);
hold on;
plot(1:size(LOC,2), im_LOC2(2,:)); axis([0 size(LOC,2) 0 8]);
title('summed LOC, L v R');

%% BP related figures
im_LOC = eeg_bp_smooth;
imLOC_L = im_LOC(channels_L,:);
imLOC_R = im_LOC(channels_R,:);
imLOC_C = im_LOC([17,18],:);

imLOC_Lk = kmeans_sane(imLOC_L,3);
imLOC_Rk =  kmeans_sane(imLOC_R,3);
imLOC_Ck =  kmeans_sane(imLOC_C,3);
imLOC_LRC = [imLOC_Lk*3; imLOC_Rk; imLOC_Ck*2];

% BP eeg by kmeans as image
subplot(ax_v, ax_h, [8 9]);
imagesc(imLOC_LRC);

% summed BP eeg thresholded by kmeans as image
subplot(ax_v, ax_h, [15 16]);
imLOC_Lk = kmeans_sane(imLOC_L,2);
imLOC_Rk = kmeans_sane(imLOC_R,2);
imLOC_Ck =  kmeans_sane(imLOC_C,2);
imLOC_LRC = [sum(imLOC_Lk,1); sum(imLOC_Rk,1); sum(imLOC_Ck,1)];
imagesc(imLOC_LRC);

%
subplot(ax_v, ax_h, [22 23]);
imLOC_LR = [sum(imLOC_Lk,1); sum(imLOC_Rk,1)];
plot(1:size(imLOC_LR,2),imLOC_LR(1,:)); xlim([0 size(imLOC_LR,2)]);
hold on;
plot(1:size(imLOC_LR,2), imLOC_LR(2,:)); xlim([0 size(imLOC_LR,2)]);
title('summed kmeans BP eeg, L v R ');

%% Find average events 
seg_avg = seg-mean(seg,1);

[LOC, avgEvents, eeg_bp_smooth] = fct_findEvents_hilo_dynamic_auto6(seg_avg, Fs, ...
    thresh_hi, thresh_lo, smoothWin, 4,20, boostWin, meanWin, c_av);

%% Periodicity calculations





%% Save figure;

%       f1 = figure('Visible','off','Position',[0 0 800 800]);
            
        thisDirec = file.folder;
        thisFilename = [file.name '_results.png'];

        suptitle(replaceUnderscore(thisFilename,'-','_'));
        set(f1,'visible','off');
        set(f1, 'PaperPositionMode','auto');
        
        thisSubDir = [thisDirec '\images\'];
        make_directories({thisSubDir});
        print(f1, '-dpng','-zbuffer','-r200', [thisSubDir thisFilename]);
        close(f1);

        
end

