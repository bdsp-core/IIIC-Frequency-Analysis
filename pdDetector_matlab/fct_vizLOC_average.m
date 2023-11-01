function [ output_args ] = fct_vizLOC_average( LOC_artifact)
%fct_vizLOC_average Show LOCs in average montage

LOC = LOC_artifact;
channelLabels_average = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};

channels_L = [1:8];
channels_R = [12:19];
channels_C = [9:11];

        im_LOC = LOC;
        im_LOC(channels_L,:)=im_LOC(channels_L,:)*3;
        im_LOC(channels_C,:)=im_LOC(channels_C,:)*2;
        
        figure(6);imagesc(im_LOC([channels_L, channels_C, channels_R],:));
        title('LOC by channel, left v right v center')
        
        im_LOC2 = [sum(LOC(channels_L,:),1)>2; sum(LOC(channels_C,:),1)>1;...
            sum(LOC(channels_R,:),1)>2];
        
        figure(8);imagesc(im_LOC2); 
        title('summed LOC >2, L v C v R');
        
        im_LOC2 = [sum(LOC(channels_L,:),1); sum(LOC(channels_C,:),1); ...
            sum(LOC(channels_R,:),1)];
        figure(7);imagesc(im_LOC2);
        title('summed LOC, L v C v R');
        
        figure(9);
        subplot(3,1,1); plot(1:size(LOC,2),im_LOC2(1,:)); axis([0 size(LOC,2) 0 8]);
        subplot(3,1,2); plot(1:size(LOC,2), im_LOC2(3,:)); axis([0 size(LOC,2) 0 8]);
        subplot(3,1,3); plot(1:size(LOC,2), im_LOC2(1,:)-im_LOC2(3,:)); axis([0 size(LOC,2) -8 8]);
        suptitle('summed LOC, L v R');
   
        LOC_diff = im_LOC2(1,:)-im_LOC2(3,:); % difference b/w L vs R. 
        LOC_diffL = sum(LOC_diff(LOC_diff>0))/sum(abs(LOC_diff));
        LOC_diffR = -1*sum(LOC_diff(LOC_diff<0))/sum(abs(LOC_diff));
        LOC_diff_rel = LOC_diffL/LOC_diffR; 
        LOC_diff_abs = LOC_diffL-LOC_diffR;   % pos = L; neg = R.
end

