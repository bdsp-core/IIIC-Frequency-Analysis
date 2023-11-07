function [ output_args ] = fct_vizEEG_bipolar( eeg_bp_smooth)
%fct_vizEEG_bipolar Visualize EEG in bipolar, along with events based on
%simple band-pass and smoothing. 
        
        figure(1);
        showEEG(eeg_bp_smooth/5);   
        
        figure(2); 
        
        channels_L = [1:4,9:12];
        channels_R = [5:8, 13:16];

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
        figure(3); imagesc([eeg_label1(1:2800),eeg_label1(2801:end)]'); 
        title('kmeans bp eeg');
           
        eeg_label2L = kmeans([eeg_bp_smooth_L]',3);  %prev 2
        eeg_label2R = kmeans([eeg_bp_smooth_R]',3);  %prev 2
        figure(11); imagesc([eeg_label2L,eeg_label2R]'); 
        title('kmeans bp eeg by channel');
        
        % kmeans on all data
        figure(12);
        im_LOC = reshape(kmeans(eeg_bp_smooth(:),3),size(eeg_bp_smooth));
        im_LOC(channels_L,:)=im_LOC(channels_L,:)*3;
        im_LOC([17,18],:)=im_LOC([17,18],:)*2;
        imagesc(im_LOC([1:4,9:12,5:8,13:18],:));
   
        figure(90); 
        channels_L = [1:4,9:12];
        channels_R = [5:8, 13:16];
        
%         im_LOC = reshape(kmeans(eeg_bp_smooth(:),3),size(eeg_bp_smooth));
        im_LOC = eeg_bp_smooth;
        imLOC_L = im_LOC(channels_L,:);
        imLOC_R = im_LOC(channels_R,:);
        imLOC_C = im_LOC([17,18],:);
        
%         imLOC_Lk = reshape(kmeans(imLOC_L(:),3),size(imLOC_L));
%         imLOC_Rk = reshape(kmeans(imLOC_R(:),3),size(imLOC_R));
%         imLOC_Ck = reshape(kmeans(imLOC_C(:),3),size(imLOC_C));
%         
        imLOC_Lk = kmeans_sane(imLOC_L,3);    
        imLOC_Rk =  kmeans_sane(imLOC_R,3);  
        imLOC_Ck =  kmeans_sane(imLOC_C,3);  
               
        imLOC_LRC = [imLOC_Lk*3; imLOC_Rk; imLOC_Ck*2];
        imagesc(imLOC_LRC);
        
        figure(91);        
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

