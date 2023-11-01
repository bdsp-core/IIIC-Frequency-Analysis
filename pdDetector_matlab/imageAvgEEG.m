function [ output_args ] = imageAvgEEG( seg_avg, seg_bi)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% put raw data in x,y,z format

% using 'seg' from freqSpike_GUI, 19x2800 nxm array of double

labels = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};
placeOnHead = [[1 2];[2 2];[3 2];[4 2];[2 1];[3 1];[4 1];[5 2]; [2 3];[3 3];[4 3];[1 4];[2 4];[3 4];[4 4];[2 5];[3 5];[4 5];[5 4]]; 
nanOnHead = [[1 1]; [1 3]; [1 5];[5 1 ];[5 3];[5 5]]; 
labelMat_head = zeros(5,5,1);

 data_head = zeros(5,5,size(seg_avg,2));
 data_head(nanOnHead(:,1),nanOnHead(:,2),:) = nan;
 
 labels_head = cell(5,5,1);
 
 Fs = 200;
%  seg_bp = eegfilt(seg,Fs,4,0);   % high pass above 4
 
 for i = 1:length(labels)
    
    thisIndx = placeOnHead(i,:);
    labels_head(thisIndx(1),thisIndx(2)) = labels(i);
    data_head(thisIndx(1),thisIndx(2),:) = seg_avg(i,:);
    
%     data_head(thisIndx(1),thisIndx(2),:) = seg_bp(i,:);
    
    
    labelMat_head(thisIndx(1),thisIndx(2)) = i;
    
 end
 
%  meanPerUnitTime = mean(reshape(data_head,5*5,size(data_head,3)),'omitnan');
%  data_head_n = data_head-repmat(reshape(meanPerUnitTime,1,1,size(meanPerUnitTime,2)),5,5);
%   img_avg = interpolateImageNan(data_head_n);

 
 img_avg = interpolateImageNan(data_head);
 
%  gvf_img = GVF3D(img_avg,0.1,100);
 
 gvf_img = GVF3D(img_avg,0.1,1);
 
 figure(64); plot(squeeze(sum(sum(img_avg(5:8,14:17,:),2),1))); 
 figure(65); plot(squeeze(sum(sum(gvf_img(5:8,14:17,:),2),1))); 
 figure(65); plot(squeeze(sum(sum(gvf_img(14:17,18:21,:),2),1))); 
 
 %%
 indx = 2015;
 figure(66); imagesc(img_avg(:,:,indx)); colormap(cool(11)); colorbar;
 figure(67); imagesc(gvf_img(:,:,indx));colormap(cool(11)); colorbar;
 figure(68); imagesc(gvf_img(:,:,indx)>0.005);colormap(cool(11)); colorbar;
 figure(69); imagesc((squeeze(sum(sum(gvf_img(5:8,18:21,:),2),1))>0.05)'); colormap(cool(11)); colorbar;
 figure(69); imagesc((squeeze(sum(sum(gvf_img(1:21,14:21,:)> (10^-3),2),1)))'); colormap(cool(11)); colorbar;
 figure(70); imagesc((squeeze(sum(sum(gvf_img(1:21,1:8,:)> (10^-3),2),1)))'); colormap(cool(11)); colorbar;
 
 %%
 
 figure(121);
 gvf_img_side = reshape(gvf_img(2:4:end,2:4:end,:),5*5,2800);
 imagesc(gvf_img_side); 
 colormap(cool(11)); colorbar;
 
 xlen = 5;
 ylen = 5;
 zlen = 2800;
 
 img_avg_side = zeros(25,2800);
 gvf_img_side = zeros(25,2800);
 L_side = zeros(25,2800);
 
 cc = bwconncomp(img_avg>10^2)
  L = labelmatrix(cc);
 
 side_indx =1;
 
 for j= 1:ylen
     
     thisY = [(j-1)*4+1:(j)*4];
     
     for m = 1:xlen
         
         thisX = [(m-1)*4+1:(m)*4];
         
         img_block = sum(reshape(img_avg(thisX, thisY,:),length(thisX)*length(thisY),zlen),1);
         gvf_block =  sum(reshape(gvf_img(thisX, thisY,:),length(thisX)*length(thisY),zlen),1);
         L_block = mean(reshape(L(thisX, thisY,:),length(thisX)*length(thisY),zlen),1);
         
         img_avg_side(side_indx,:) = img_block;
         gvf_side(side_indx,:) = gvf_block;
         L_side(side_indx,:) = L_block;
         
         side_indx= side_indx+1;
     end
 end
 
 fake_elec = [1,5,11,15,21,25];
 head_labels = labels_head(:);
 head_labels = head_labels(~ismember(1:25,fake_elec),:);
 
%  theseElecs = ismember(labels,cat(2,labels_head(:)));
 [~,theseElecs ]= ismember(labels,head_labels);
 
 img_avg_side = img_avg_side(~ismember(1:25,fake_elec),:)
 img_avg_side = img_avg_side(theseElecs,:);
 gvf_side = gvf_side(~ismember(1:25,fake_elec),:)
 gvf_side = gvf_side(theseElecs,:);
 
 L_side = L_side(~ismember(1:25,fake_elec),:);
 L_side = L_side(theseElecs,:);
 
 figure(14); imagesc(img_avg_side);
 figure(15); imagesc(gvf_side);
 figure(19); imagesc(img_avg_side>10^2);
 figure(17); imagesc(gvf_side>10^(-1.5));
 figure(38); imagesc(L_side)
 
%  figure; imagesc(kmeans_sane(img_avg_side,2));
%  figure; imagesc(kmeans_sane(seg_avg,2));
 
 %%
%  C= xcorr2(reshape(img_avg,21*21,2800));
 
   m=3; j=3;
  thisY = [(j-1)*4+1:(j)*4];
   thisX = [(m-1)*4+1:(m)*4];
  
 C= xcorr2(reshape(img_avg(thisX,thisY,:),length(thisX)*length(thisY),2800));
 
%  [C1,lags]= xcorr(reshape(img_avg,21*21,2800), reshape(img_avg,21*21,2800),Fs*5);
%  [C1,lags]= xcorr(reshape(img_avg,21*21,2800),Fs*5);
 
C_n = C./max(C(:));

figure;imagesc(C_n');
figure;plot(C_n(length(thisX)*length(thisY),2800:2800+(Fs*5)));
 %%
 
 stats = regionprops('table',L);
 
 
 %%
 stats = regionprops('table',gvf_img>(10^-3),'Area','Centroid','BoundingBox','PixelIdxList');
 
 theseCentroids = cat(1,stats.Centroid);
 theseCentroids_loc = zeros(1,2800);
 theseCentroids_loc(round(theseCentroids(:,3)))=1;
 
 figure;imagesc(theseCentroids_loc);
 figure;scatter(theseCentroids(:,3), theseCentroids(:,2))
 figure;gscatter(theseCentroids(:,3), theseCentroids(:,1),[1:size(theseCentroids,1)])
 figure;gscatter(theseCentroids(:,3), stats.Area,[1:size(theseCentroids,1)])
 xlim([1 2800])
  %%
  cc = bwconncomp(img_avg>10^3.5)
  L = labelmatrix(cc);
  
  stats.GroupNum = num2cell(1:size(stats,1))';
  
  stats_b = stats(logical(kmeans_sane(theseCentroids(:,1),2)),:)
  
  
 %%
 
labelMat_head_io = interpolateImage(labelMat_head);
 
 trange = [1:10:2800];
%  for l=1:length(trange)
%  figure(64); imagesc(img_avg(:,:,trange(l)),[-200 200]); 
%  colormap(cool(11)); colorbar;
%  figure(61); showEEG(seg_bi);
%  gca;
%  hold on;
%  line([trange(l) trange(l+1)],[0 100],'LineWidth', 2, 'Color','r');
%  
%  figure(63); showEEG(seg([1:9,11:end],:)-mean(seg([1:9,11:end],:),1));
%  gca;
%  hold on;
%  line([trange(l) trange(l+1)],[0 100],'LineWidth', 2, 'Color','r');
%  
%  waitforbuttonpress;
%  end

l=1;
while l<length(trange)
    crange = [-100 100];
    
    figure(64); imagesc(img_avg(:,:,trange(l)),crange); colormap(cool(11)); colorbar;
    
    figure(65);
    s=surf(1:41,1:41,img_avg(:,:,trange(l)));
    colormap(cool(11));
    axis([1 41 1 41]);
    zlim(crange);
    caxis(crange);
%     campos(1e3*[-0.1404    0.3099    1.0244]);
    campos([21.0590  367.3533  -31.3796]);
    
%     s.EdgeColor = 'none';
    s.FaceAlpha =0.5;
    colorbar;
%     set(gca,'Xdir','reverse','Ydir','reverse')
    set(gca,'Ydir','reverse');
    
    figure(61); showEEG(seg_bi);
    gca;
    hold on;
    line([trange(l) trange(l+1)],[0 100],'LineWidth', 2, 'Color','r');
    
    f63=figure(63); 
    fct_showEEG_avg(seg_avg);
    gca;
    hold on;
    line([trange(l) trange(l+1)],[0 100],'LineWidth', 2, 'Color','r');
    
    figure(67);
    [px,py] = gradient(img_avg(:,:,trange(l)));  
    contour(1:41,1:41,img_avg(:,:,trange(l)));
    set(gca,'Ydir','reverse');
    colormap(cool(11));
    hold on
    quiver(1:41,1:41,px,py);
    hold off
    
%     L_vec_x = sum(sum(px(:,1:21)));
%     L_vec_y = sum(sum(py(:,1:21)));
%     
%     L_avg_vec = sqrt(L_vec_y.^2 + L_vec_x.^2);
%     L_avg_vec_dir = atan(-1*L_vec_y/L_vec_x)*(180/pi);   % heading wrt x-axis. yaxis is flipped.
%     
%     R_vec_x = sum(sum(px(:,21:end)));
%     R_vec_y = sum(sum(py(:,21:end)));
%     
%     R_avg_vec = sqrt(R_vec_y.^2 + R_vec_x.^2);
%     R_avg_vec_dir = atan(-1*R_vec_y/R_vec_x)*(180/pi);    % heading wrt x-axis. yaxis is flipped.
%     
%     [L_avg_vec R_avg_vec] 
%     [L_avg_vec_dir R_avg_vec_dir]
%     
    
    regions = {{1:41,1:21}; {1:41, 21:41}; {1:21,1:41}; {21:41,1:41}; {1:41,1:41}};  % LEFT, RIGHT, ANT, POST, TOTAL
    
    for i=1:size(regions,1)
       vec_x{i} = sum(sum(px(regions{i}{1},regions{i}{2})));
       vec_y{i} = sum(sum(py(regions{i}{1},regions{i}{2})));
       avg_vec{i} =  sqrt(vec_x{i}.^2 + vec_y{i}.^2);
       avg_vec_dir{i} = atan2(-1*vec_y{i},vec_x{i})*(180/pi);    % heading wrt x-axis. yaxis is flipped.
    end
    
    avg_vec',avg_vec_dir'
    
    
    w = waitforbuttonpress;
    
    switch w
        case 0 %mouse click
         
            [xi,yi] = getpts(f63); 
            l = find(trange>xi(1),1);
            
            x1 = floor(xi(1));
            x2 = floor(xi(2));
        
            hold on;
            thisIndx = x1:1:x2;
            for j=1:length(thisIndx)
                
             [px,py] = gradient(img_avg(:,:,thisIndx(j)));
             
             for i=1:size(regions,1)
                 vec_x{j,i} = sum(sum(px(regions{i}{1},regions{i}{2})));
                 vec_y{j,i} = sum(sum(py(regions{i}{1},regions{i}{2})));
                 
                 avg_vec{j,i} =  sqrt(vec_x{j,i}.^2 + vec_y{j,i}.^2);
                 avg_vec_dir{j,i} = atan2(-1*vec_y{j,i},vec_x{j,i})*(180/pi);    % heading wrt x-axis. yaxis is flipped.
             end
             
%              compass(vec_x{j,5}, vec_y{j,5});
            end
            hold off;
            
            figure(70);
            clear ax;
            ax(1)=subplot(711);
            vecMag = cat(1,avg_vec{1:end,5});
            indx_vecMag_thresh = vecMag>1000;
            plot(1:j,vecMag);
            axis tight;
            ax(2)=subplot(712);
            vecDirec = cat(1,avg_vec_dir{1:end,5});
            
            vecDirec(~indx_vecMag_thresh) = nan;
            plot(1:j,vecDirec);
            
            ylim([-180 180]);
            axis tight;
            
            [xlen,ylen]=size(img_avg(:,:,1));
           
            %max position. 
%             img_max_bin = zeros(size(img_avg(:,:,x1:x2)));
            max_val = zeros(length(thisIndx),1);
            max_loc = max_val;
            for j=1:length(thisIndx)
                thisPage = img_avg(:,:,thisIndx(j));
                [m, m_l] = max(thisPage(:));
                max_val(j) = m;
                max_loc(j) =m_l;
            end
          
            % min position.
            clear m; clear m_l;
%             img_min_bin = nan(size(img_avg(:,:,x1:x2)));
            min_val = zeros(length(thisIndx),1);
            min_loc = min_val;
            for j=1:length(thisIndx)
                thisPage = img_avg(:,:,thisIndx(j));
                [m, m_l] = min(thisPage(:));
                min_val(j) = m;
                min_loc(j) =m_l;
            end
            
           
            [max_loc_y,max_loc_x] = ind2sub(size(img_avg(:,:,1)),max_loc);
                        
            [min_loc_y, min_loc_x] = ind2sub(size(img_avg(:,:,1)),min_loc);
            
            ax(3)=subplot(713);
%             plot(thisIndx, max_val, thisIndx, min_val);   axis tight;
            plot(1:j, max_val, 1:j, min_val);   axis tight;
              
            subplot(714);
%             max_loc_y(~indx_vecMag_thresh)=nan;
            plot(1:j,-1*(max_loc_y-41)); ylim([1 41]);
            
            hold on;
            plot(1:j, ones(j,1)*xlen/3);
            plot(1:j, ones(j,1)*xlen/3*2); axis tight; hold off;
            
            subplot(715);
%             max_loc_x(~indx_vecMag_thresh)=nan;
            plot(1:j,max_loc_x); ylim([1 41]);
            hold on;
            plot(1:j, ones(j,1)*xlen/3);
            plot(1:j, ones(j,1)*xlen/3*2);axis tight; hold off;
            subplot(716);
%             min_loc_y(~indx_vecMag_thresh)=nan;
            plot(1:j,-1*(min_loc_y-41)); ylim([1 41]);
            hold on;
            plot(1:j, ones(j,1)*xlen/3);
            plot(1:j, ones(j,1)*xlen/3*2);axis tight; hold off;
            subplot(717);
%             min_loc_x(~indx_vecMag_thresh)=nan;
            plot(1:j,min_loc_x); ylim([1 41]);
            hold on;
            plot(1:j, ones(j,1)*xlen/3);
            plot(1:j, ones(j,1)*xlen/3*2);axis tight; hold off;
            
            
            figure(72);
            
            ax(4)=subplot(511); 
%             imagesc([max_loc_x';-1*(max_loc_y-41)'],[1 41]); colormap(cool(3));
%             imagesc([-1*(max_loc_y-41)';-1*(min_loc_y-41)'],[1 41]); colormap(cool(3));
            

            max_AP = [(max_loc_y<ylen/3*1)'; (max_loc_y>(ylen/3) & max_loc_y<(ylen/3*2))'; (max_loc_y>(ylen/3*2))'];
            imagesc(max_AP); 

            subplot(512);
            min_AP = [(min_loc_y<ylen/3*1)'; (min_loc_y>(ylen/3) & min_loc_y<(ylen/3*2))'; (min_loc_y>(ylen/3*2))'];
            imagesc(min_AP); 
            
%             max_loc_y_nonan = max_loc_y;
%             max_loc_y_nonan(isnan(max_loc_y)) = -9;
%             min_loc_y_nonan = min_loc_y;
%             min_loc_y_nonan(isnan(min_loc_y)) = -9;
            
%             imagesc([-1*(max_loc_y-41)';-1*(min_loc_y-41)'],[1 41]); colormap(cool(4));
            
            
            ax(6)=subplot(513); 
%             imagesc([min_loc_x'; min_loc_x'],[1 41]); colormap(cool(3));
%             plot(thisIndx, max_loc_y, thisIndx, min_loc_y);   axis tight;
%             plot(1:j, max_loc_y, 1:j, min_loc_y);   axis tight;

            max_LR= [(max_loc_x<xlen/3*1)'; (max_loc_x>(xlen/3) & max_loc_x<(xlen/3*2))'; (max_loc_x>(xlen/3*2))'];
            imagesc(max_LR); 

            subplot(514);
            min_LR = [(min_loc_x<xlen/3*1)'; (min_loc_x>(xlen/3) & min_loc_x<(xlen/3*2))'; (min_loc_x>(xlen/3*2))'];
            imagesc(min_LR); 

            
            ax(7)=subplot(515);
            min2max_heading = atan2(-1*(max_loc_y-min_loc_y),(max_loc_x-min_loc_x))*(180/pi); 
%             plot(thisIndx,min2max_heading);   axis tight;
            plot(1:j,min2max_heading);   axis tight;
            ylim([-180 180]);
            
%             linkaxes(ax,'x');
            
            figure(74);
            % max-min method with threshold/centroids.
%             
            img_avg_max_bw = img_avg>50;
            img_avg_min_bw = img_avg<-50;

%             
            stats_max=cell(length(thisIndx),1);
            stats_min=cell(length(thisIndx),1);
            centroids_max = nan(length(thisIndx),2);
            centroids_min = nan(length(thisIndx),2);
            
            for j=1:length(thisIndx)
                % takes the largest centroid
                thisMax = regionprops(img_avg_max_bw(:,:,thisIndx(j)),'Area','Centroid','PixelIdxList');
                [m,indx] = max(cat(1,thisMax.Area));
                stats_max{j} = thisMax(indx);
                if isempty(thisMax(indx))
                centroids_max(j,1:2) = [nan,nan];
                else
                   centroids_max(j,1:2) = thisMax(indx).Centroid; 
                end
                
                thisMin = regionprops(img_avg_min_bw(:,:,thisIndx(j)),'Area','Centroid','PixelIdxList');
                [m,indx] = max(cat(1,thisMin.Area));
                stats_min{j} = thisMin(indx);
                
                if isempty(thisMin(indx))
                    centroids_min(j,1:2) = [nan,nan];
                else
                    centroids_min(j,1:2) = thisMin(indx).Centroid;
                end
                
            end
%             
%             centroids_max = cat(1,stats_max{:});
%             centroids_max = cat(1,centroids_max.Centroid);
%             
%             centroids_min= cat(1,stats_min{:});
%             centroids_min = cat(1,centroids_min.Centroid);
            
            %max-AP, min-AP
            subplot(412);
            imagesc([-1*(centroids_max(:,2)-41)';-1*(centroids_min(:,2)-41)'],[1 41]); colormap(cool(3));

            %max-LR, min-LR
            subplot(413);
            imagesc([centroids_max(:,1)';centroids_min(:,1)'],[1 41]); colormap(cool(3));
            
            % vector heading
            subplot(414);
            min2max_heading = atan2(-1*(centroids_max(:,2)-centroids_min(:,2)),(centroids_max(:,1)-centroids_min(:,1)))*(180/pi); 
%             plot(1:j,min2max_heading,'LineWidth',2);
            scatter(1:j,min2max_heading);
            ylim([-180 180]);
            xlim([0 length(thisIndx)]);
            
            
%             figure(81);
%             minMax = prctile(reshape(img_avg(:,:,thisIndx),41*41,length(thisIndx)),[1 99]);
%             subplot(411);
%             plot(1:j,minMax(1,:), 1:j, minMax(2,:)); axis tight;
%             
%             img_avg_max_bw = gt(img_avg(:,:,thisIndx),repmat(reshape(minMax(2,:),1,1,length(thisIndx)),41,41));
%             img_avg_min_bw = lt(img_avg(:,:,thisIndx),repmat(reshape(minMax(1,:),1,1,length(thisIndx)),41,41));
%             
%             
            
            figure(71);
%             seg_avg = seg_avg-mean(seg_avg,1);
            showEEG_avg(seg_avg(:,x1:x2));
%             
%             figure(72);
%             seg_bp_avg = seg_bp-mean(seg_bp,1);
%             showEEG_avg(seg_bp_avg(:,x1:x2));
            
            figure(83);
            imagesc(kmeans(smooth(vecMag,5),2)');
            
        case 1 %keyboard
            key = get(gcf,'currentcharacter');
            switch key
                case 27 % escape key
                    break
                case 28
                    l=l-1;
                case 29
                    l=l+1;     
            end
    end
end
end

