function [ output_args ] = fct_vizEEG_average( eeg_bp_smooth)
%fct_vizEEG_average Visualize EEG in bipolar, along with events based on
%simple band-pass and smoothing. 
        
        figure(90);
        fct_showEEG_avg(eeg_bp_smooth/5);   
        
        figure(91); 
        
        channels_L = [1:8];
        channels_R = [12:19];
        channels_C = [9:11];

        eeg_bp_smooth_L = sum(eeg_bp_smooth(channels_L,:),1);
        eeg_bp_smooth_R = sum(eeg_bp_smooth(channels_R,:),1);
        
        im_eeg_bp = [eeg_bp_smooth_L; eeg_bp_smooth_R];
        ax1=subplot(2,1,1); plot(1:size(im_eeg_bp,2),im_eeg_bp(1,:)); xlim([0 size(eeg_bp_smooth,2)]);
        ylabel('left');
        ax2=subplot(2,1,2); plot(1:size(im_eeg_bp,2), im_eeg_bp(2,:)); xlim([0 size(eeg_bp_smooth,2)]);
        ylabel('right');
        linkaxes([ax1,ax2],'xy');
        
        suptitle('sum of BP eeg');
        
        eeg_label1 = kmeans([eeg_bp_smooth_L eeg_bp_smooth_R]',3);   %prev 2
        figure(92); imagesc([eeg_label1(1:2800),eeg_label1(2801:end)]'); 
        title('kmeans bp eeg');
           
        eeg_label2L = kmeans([eeg_bp_smooth_L]',3);  %prev 2
        eeg_label2R = kmeans([eeg_bp_smooth_R]',3);  %prev 2
        figure(93); imagesc([eeg_label2L,eeg_label2R]'); 
        title('kmeans bp eeg by channel');
        
        % kmeans on all data
        figure(94);
        im_LOC = reshape(kmeans(eeg_bp_smooth(:),3),size(eeg_bp_smooth));
        im_LOC(channels_L,:)=im_LOC(channels_L,:)*3;
        im_LOC(channels_C,:)=im_LOC(channels_C,:)*2;
        imagesc(im_LOC([channels_L,channels_R,channels_C],:));
   
        figure(95); 
        
%         im_LOC = reshape(kmeans(eeg_bp_smooth(:),3),size(eeg_bp_smooth));
        im_LOC = eeg_bp_smooth;
        imLOC_L = im_LOC(channels_L,:);
        imLOC_R = im_LOC(channels_R,:);
        imLOC_C = im_LOC(channels_C,:);
        
        imLOC_Lk = kmeans_sane(imLOC_L,3);    
        imLOC_Rk =  kmeans_sane(imLOC_R,3);  
        imLOC_Ck =  kmeans_sane(imLOC_C,3);  
               
        imLOC_LRC = [imLOC_Lk*3; imLOC_Rk; imLOC_Ck*2];
        imagesc(imLOC_LRC);
        
        figure(96);        
        imLOC_Lk = kmeans_sane(imLOC_L,2);
        imLOC_Rk = kmeans_sane(imLOC_R,2);
        imLOC_LR = [sum(imLOC_Lk,1); sum(imLOC_Rk,1)];
        ax1=subplot(2,1,1); plot(1:size(imLOC_LR,2),imLOC_LR(1,:)); xlim([0 size(imLOC_LR,2)]);
        ylabel('left');
        ax2=subplot(2,1,2); plot(1:size(imLOC_LR,2), imLOC_LR(2,:)); xlim([0 size(imLOC_LR,2)]);
        ylabel('right');
        linkaxes([ax1,ax2],'xy');
        
        suptitle('sum of LOC from kmeans BP eeg by channel');
        

end

