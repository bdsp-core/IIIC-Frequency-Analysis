%% take segmentation of EEG, compare events based on "distance" statistics

avg_vec = cell(size(img_avg,3),1);
avg_vec_dir = cell(size(img_avg,3),1);

for j=1:size(img_avg,3)

    [px,py] = gradient(img_avg(:,:,j));
    
    vec_x= sum(sum(px));
    vec_y= sum(sum(py));
    avg_vec{j} =  sqrt(vec_x.^2 + vec_y.^2);
    avg_vec_dir{j} = atan2(-1*vec_y,vec_x)*(180/pi);    % heading wrt x-axis. yaxis is flipped.

end

% max-min value, position
[xlen,ylen]=size(img_avg(:,:,1));

%max position.
%             img_max_bin = zeros(size(img_avg(:,:,x1:x2)));
max_val = zeros(size(img_avg,3),1);
max_loc = max_val;
for j=1:length(max_val)
    thisPage = img_avg(:,:,j);
    [m, m_l] = max(thisPage(:));
    max_val(j) = m;
    max_loc(j) =m_l;
end

% min position.
clear m; clear m_l;
%             img_min_bin = nan(size(img_avg(:,:,x1:x2)));
min_val = zeros(size(img_avg,3),1);
min_loc = min_val;
for j=1:length(thisIndx)
    thisPage = img_avg(:,:,j);
    [m, m_l] = min(thisPage(:));
    min_val(j) = m;
    min_loc(j) =m_l;
end

[max_loc_y, max_loc_x] = ind2sub(size(img_avg(:,:,1)),max_loc);
[min_loc_y, min_loc_x] = ind2sub(size(img_avg(:,:,1)),min_loc);

vecMag = cat(1,avg_vec{:});
thisSeg = kmeans(smooth(vecMag,200),2)';

theseElements = unique(thisSeg);
vecMag_seg = zeros(length(theseElements),1);
for j=1:theseElements
    vecMag_seg(j)= sum(vecMag(thisSeg==theseElements(j)));   
end
[m,m_indx]=max(vecMag_seg);

stats = regionprops((thisSeg==m_indx), 'Area','BoundingBox','PixelIdxList');
gradientArray = cell(length(stats),4);


for k=1:length(stats)
    
    idx_stat = stats(k).PixelIdxList;
    
%     [px,py] = gradient(img_avg(:,:,idx_stat));
%     
%     vec_x= sum(sum(px));
%     vec_y= sum(sum(py));
%     avg_vec{j} =  sqrt(vec_x.^2 + vec_y.^2);
%     avg_vec_dir{j} = atan2(-1*vec_y,vec_x)*(180/pi);    % heading wrt x-axis. yaxis is flipped.
    
    gradientArray{k,1}= vecMag(idx_stat);
    gradientArray{k,2}= [max_val(idx_stat) min_val(idx_stat)];
    gradientArray{k,3}= [max_loc_x(idx_stat), max_loc_y(idx_stat), min_loc_x(idx_stat), min_loc_y(idx_stat)];
end

%%
cov(gradientArray{1,2})

tic
[autocor,lags] = xcorr(thisSeg',5*Fs,'coeff');
toc
figure;
plot(lags/Fs,autocor)
ylim([0 1]);