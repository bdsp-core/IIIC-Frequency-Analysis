function [ thisResult ] = fct_findVectorEvents_hiloprct( thisVal, thresh_hi, thresh_lo, thisLabel, max_xy_min_xy)
%fct_findVectorEvents_hiloprct 
%   ***
           max_xy_min_xy=max_xy_min_xy';
           
           if thresh_hi ==0
               eeg_prctl_hi = 0;
           else
               eeg_prctl_hi = prctile(thisVal, thresh_hi);
           end
           if thresh_lo ==0
               eeg_prctl_lo =0;
           else
               eeg_prctl_lo = prctile(thisVal, thresh_lo);
           end
           
           % high-low threshold
           bw_hi = thisVal>eeg_prctl_hi;
           bw_lo = thisVal>eeg_prctl_lo;
           
           min_vec_hi = regionprops(bw_hi,'PixelIdxList');
%            lpds_lo = regionprops('table',thisEEG_bp_smooth>thresh_lo,'PixelIdxList');
           
           indx_hi = cat(1, min_vec_hi.PixelIdxList);
%            indx_lo = cat(1, lpds_lo.PixelIdxList);
           
           isect = bwselect([bw_lo; bw_lo]', ones(length(indx_hi),1), indx_hi,4);   % must be 2D, order is critical, takes binary image but needs pixel ids for selection
           isect = isect(:,1)';
           
           vecs_ = regionprops('table',isect,thisVal,'Area',...
               'BoundingBox','PixelIdxList','MaxIntensity','MeanIntensity','PixelValues','WeightedCentroid');
           
           vecs_max_x = regionprops('table',isect,max_xy_min_xy(1,:),...
               'MaxIntensity','MeanIntensity','PixelValues','WeightedCentroid');
           vecs_max_x = appendTblNames(vecs_max_x,'_max_x');
           
           vecs_max_y= regionprops('table',isect,max_xy_min_xy(2,:),...
               'MaxIntensity','MeanIntensity','PixelValues','WeightedCentroid');
           vecs_max_y = appendTblNames(vecs_max_y,'_max_y');
           
           vecs_min_x= regionprops('table',isect,max_xy_min_xy(3,:),...
               'MaxIntensity','MeanIntensity','PixelValues','WeightedCentroid');
           vecs_min_x = appendTblNames(vecs_min_x,'_min_x');
           
           vecs_min_y= regionprops('table',isect,max_xy_min_xy(4,:),...
               'MaxIntensity','MeanIntensity','PixelValues','WeightedCentroid');
           vecs_min_y = appendTblNames(vecs_min_y,'_min_y');
           
           thisResult = cat(2, vecs_ , vecs_max_x,vecs_max_y,vecs_min_x,vecs_min_y);
           
           thisResult.type = repmat({thisLabel},size(vecs_ ,1),1);           
           thisResult.channel = repmat({'vector'},size(vecs_,1),1);
           thisResult.channelNum = repmat(1,size(vecs_,1),1);
           thisResult.channelEventNum = num2cell(1:size(vecs_ ,1))';

end

