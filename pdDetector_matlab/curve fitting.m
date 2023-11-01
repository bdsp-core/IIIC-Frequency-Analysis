figure; periodogram(eeg_label1(1:2800),[],[0:0.1:3],Fs);

period_est = 500;  % samples
%         X = [1:2800].*(2*pi/period_est);
X = [1:2800];
y = eeg_bp_smooth_L;
X2 = [1:2800]./2800;


%         modelfun = @(b,c,x,phi,H)b + (H/(2*(1-c)))*(cos(x-phi)-c + abs(cos(x-phi)-c));
modelfun_b = @(b,x,x2)b(1) + (b(2)./(2.*(1-b(3)))).*(cos(x.*b(5)-b(4))-b(3) + abs(cos(x.*b(5)-b(4))-b(3)));
%b1 = b; b2 = H, b3= c, b4 =phi, b5= # of cycles
beta0 = [200,2000,0.9,0,(2*pi/period_est)];
opts.RobustWgtFun = 'bisquare';
%         mdl=fitnlm(X,y,modelfun_b,beta0,'ErrorModel', 'constant', 'Options',opts);
mdl=fitnlm(X,y,modelfun_b,beta0,'ErrorModel', 'constant');


figure(2);
plot(X,y);
hold on;
plot(X,mdl.Fitted,'k');
%%
% try with thresholded data
period_est = 500;  % samples
%         X = [1:2800].*(2*pi/period_est);
X = [1:2800];
y = eeg_bp_smooth_L>2000;

%         modelfun = @(b,c,x,phi,H)b + (H/(2*(1-c)))*(cos(x-phi)-c + abs(cos(x-phi)-c));
modelfun_b = @(b,x)b(1) + (b(2)./(2.*(1-b(3)))).*(cos(x.*b(5)-b(4))-b(3) + abs(cos(x.*b(5)-b(4))-b(3)));
%b1 = b; b2 = H, b3= c, b4 =phi, b5= # of cycles
beta0 = [0,1,0.9,0,(2*pi/period_est)];
mdl=fitnlm(X,y,modelfun_b,beta0, 'ErrorModel', 'constant');

figure(2);
plot(X,y);
hold on;
plot(X,mdl.Fitted,'k');

figure(5);
showEEG(seg_bi);

%% stacked matrix

period_est = 500;  % samples
%         X = [1:2800].*(2*pi/period_est);
X = [1:1250];
y = eeg_bp_smooth_L(1:1250)>2000;

%         modelfun = @(b,c,x,phi,H)b + (H/(2*(1-c)))*(cos(x-phi)-c + abs(cos(x-phi)-c));
modelfun_b = @(b,x)b(1) + (b(2)./(2.*(1-b(3)))).*(cos(x.*b(5)-b(4))-b(3) + abs(cos(x.*b(5)-b(4))-b(3)));
%b1 = b; b2 = H, b3= c, b4 =phi, b5= # of cycles
beta0 = [0,1,0.9,0,(2*pi/period_est)];
mdl=fitnlm(X,y,modelfun_b,beta0, 'ErrorModel', 'constant');

figure(6);
plot(X,y);
hold on;
plot(X,mdl.Fitted,'k');
%%
tic
[autocor,lags] = xcorr(eeg_bp_smooth,5*Fs,'coeff');
toc
figure;
plot(lags/Fs,autocor)

%% downsample to 10 hz for cross-correlation

eeg_bp_smooth_lo = downsample(eeg_bp_smooth_discharge',20)';

showEEG(eeg_bp_smooth_lo);

Fs_n = 10;

tic
[autocor,lags] = xcorr(eeg_bp_smooth_lo',5*Fs_n,'coeff');
toc
figure;
plot(lags/Fs_n,autocor)
ylim([0 1]);

for i=1:size(autocor,2)
    [pks, locs, w, p] = findpeaks(autocor(:,i),'MinPeakDistance',0.5*Fs_n);
    pks_c{i} = pks;
    locs_c{i} = locs;
    w_c{i} = w;
    p_c{i} = p;
    
end

p_c_sq = reshape(p_c,18,18);
pks_c_sq = reshape(pks_c, 18,18);
locs_c_sq = reshape(locs_c,18,18);

lags_sq = cellfun(@(x,y) (x(y)./Fs_n)', repmat({lags},18,18), locs_c_sq,'UniformOutput',0);

[lags(locs)/Fs_n;p']'

%% routine to find zero-lag peaks across channels

%for each channel-by-channel cross-correlation, check if zero-lag peak
%exists, and if so, report the prominence. 

zerolag_p = zeros(18,18);

for i=1:18
    
    theseLags = lags_sq(:,i);
    
    for j=1:18
        
        thisLag = theseLags{j};
        
        indx = find(thisLag==0);
        if ~isempty(indx)
            thatProminence = p_c_sq{i,j};
            zerolag_p(i,j)=thatProminence(indx);
        else
            zerolag_p(i,j) = nan;
        end
       
   end
    
end


% figure;
% plot(lags/Fs_n,autocor(:,20))

%% simpler, correlation

[r,p] = corr(eeg_bp_smooth_lo);

[r,p] = corrcoef(eeg_bp_smooth_lo');
c


%% put raw data in x,y,z format

% using 'seg' from freqSpike_GUI, 19x2800 nxm array of double

labels = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};
placeOnHead = [[1 2];[2 2];[3 2];[4 2];[2 1];[3 1];[4 1];[5 2]; [2 3];[3 3];[4 3];[1 4];[2 4];[3 4];[4 4];[2 5];[3 5];[4 5];[5 4]]; 
nanOnHead = [[1 1]; [1 3]; [1 5];[5 1 ];[5 3];[5 5]]; 
labelMat_head = zeros(5,5,1);

 data_head = zeros(5,5,size(seg,2));
 data_head(nanOnHead(:,1),nanOnHead(:,2),:) = nan;
 
 labels_head = cell(5,5,1);
 
 for i = 1:length(labels)
    
    thisIndx = placeOnHead(i,:);
    labels_head(thisIndx(1),thisIndx(2)) = labels(i)
    data_head(thisIndx(1),thisIndx(2),:) = seg(i,:);
    labelMat_head(thisIndx(1),thisIndx(2)) = i;
    
 end
 
 meanPerUnitTime = mean(reshape(data_head,5*5,size(data_head,3)),'omitnan');
 data_head_n = data_head-repmat(reshape(meanPerUnitTime,1,1,size(meanPerUnitTime,2)),5,5);
 
 img_avg = interpolateImageNan(data_head_n);
 
labelMat_head_io = interpolateImage(labelMat_head);
 
 trange = [1:10:2800];
 for l=1:length(trange)
 figure(64); imagesc(img_avg(:,:,trange(l)),[-200 200]); 
 colormap(cool(11)); colorbar;
 figure(61); showEEG(seg_bi);
 gca;
 hold on;
 line([trange(l) trange(l+1)],[0 100],'LineWidth', 2, 'Color','r');
 
 figure(63); showEEG(seg([1:9,11:end],:)-mean(seg([1:9,11:end],:),1));
 gca;
 hold on;
 line([trange(l) trange(l+1)],[0 100],'LineWidth', 2, 'Color','r');
 
 waitforbuttonpress;
 end
 
%%
 figure; periodogram(double(eeg_bp_smooth_L),[],[0:0.1:3],Fs);
 figure; imagesc(eeg_bp_smooth_L>2000);
 
 
 %% interevent interval
  
%%
[Pxx,F] = periodogram(eeg_bp_smooth_lo,[],240,Fs);

plot(F,Pxx)
grid
xlabel('Cycles/Hour')
title('Periodogram of Sleep States')

%%
[seg_bi_psd, f]=cwt(seg_bi(channels_L(4),:)',hanning(Fs),[],[0:0.1:3],Fs);
%         [seg_bi_psd, f]=cwt(eeg_bp_smooth_L',hanning(Fs),[],[0:0.01:1,1:0.1:3],Fs);
figure;
[X,Y]= meshgrid(1:2800,f);
s=surf(X,Y,abs(seg_bi_psd));
s.EdgeColor = 'none';

figure;
plot(f,sum(abs(seg_bi_psd),2));


%% extra removed from freqSpike_GUIv3
    figure; periodogram(eeg_label1(1:2800),[],[0:0.1:3],Fs);
        
        period_est = 500;  % samples
%         X = [1:2800].*(2*pi/period_est);
        X = [1:2800];
        y = eeg_bp_smooth_L;
        X2 = [1:2800]./2800;
   
        
%         modelfun = @(b,c,x,phi,H)b + (H/(2*(1-c)))*(cos(x-phi)-c + abs(cos(x-phi)-c));
        modelfun_b = @(b,x,x2)b(1) + (b(2)./(2.*(1-b(3)))).*(cos(x.*b(5)-b(4))-b(3) + abs(cos(x.*b(5)-b(4))-b(3)));
               %b1 = b; b2 = H, b3= c, b4 =phi, b5= # of cycles
        beta0 = [200,2000,0.9,0,(2*pi/period_est)];
        opts.RobustWgtFun = 'bisquare';
%         mdl=fitnlm(X,y,modelfun_b,beta0,'ErrorModel', 'constant', 'Options',opts);
        mdl=fitnlm(X,y,modelfun_b,beta0,'ErrorModel', 'constant');
                     
        
         figure(2);
         plot(X,y);
         hold on;
         plot(X,mdl.Fitted,'k');
             
         % try with thresholded data    
         period_est = 500;  % samples
%         X = [1:2800].*(2*pi/period_est);
        X = [1:2800];
        y = eeg_bp_smooth_L>2000;
        
%         modelfun = @(b,c,x,phi,H)b + (H/(2*(1-c)))*(cos(x-phi)-c + abs(cos(x-phi)-c));
        modelfun_b = @(b,x)b(1) + (b(2)./(2.*(1-b(3)))).*(cos(x.*b(5)-b(4))-b(3) + abs(cos(x.*b(5)-b(4))-b(3)));
               %b1 = b; b2 = H, b3= c, b4 =phi, b5= # of cycles
        beta0 = [0,1,0.9,0,(2*pi/period_est)];
             mdl=fitnlm(X,y,modelfun_b,beta0, 'ErrorModel', 'constant');
                     
         figure(2);
         plot(X,y);
         hold on;
         plot(X,mdl.Fitted,'k');
         
        figure(5); 
        showEEG(seg_bi);

        [seg_bi_psd, f]=cwt(seg_bi(channels_L(4),:)',hanning(Fs),[],[0:0.1:3],Fs);
%         [seg_bi_psd, f]=cwt(eeg_bp_smooth_L',hanning(Fs),[],[0:0.01:1,1:0.1:3],Fs);
        figure;
        [X,Y]= meshgrid(1:2800,f);
        s=surf(X,Y,abs(seg_bi_psd));
        s.EdgeColor = 'none';
        
        figure;
        plot(f,sum(abs(seg_bi_psd),2));