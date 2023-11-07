function [loc_, lpds_new2] = fct_detectLPDs_v1(seg_bi, LOC, thresh_samples)
      %INPUT: 
      %    seg_bi,  nxm matrix of double, raw EEG data in bipolar montage. 
      %    LOC obj, nxm matrix of logicals, where n= # of channels, m=samples from eeg data
      %    thresh_samples,   a minimum duration threshold for events,
      %    defined in samples. 
      %OUTPUT: 
      %    loc_ obj, 1xm matrix of logical data labeling 'true' LPDs in
      %          sample-space
      %    side,  a double representing sided-ness of LPDs, range [-1 to 1]
      %          L is positive, R is negative; 0 would be balanced, but
      %          these events do not represent LPDs.
        
      loc_ = zeros(1,size(LOC,2));   % no events that meet LPD criteria.
      lpds_new2 = NaN;
        
      
      %sum events for left sided channels
      % in bipolar montage, corresponds to channels 1-4, 9-12
      channels_L = [1:4,9:12];
      loc_L = sum(LOC(channels_L,:),1);
      
      %right channels
      %n.b. this approach ignores the Fz-Cz, Cz-Pz channels; could
      % double-count them with both L and R
      channels_R = [5:8, 13:16];
      loc_R = sum(LOC(channels_R,:),1);
      
      max_loc_L = max(loc_L);
      max_loc_R = max(loc_R);
      min_chan_present = 2;
      
      % must be present in at least 2 channels
      if (max_loc_L >= min_chan_present || max_loc_R >= min_chan_present)
          
          % if 0-1 events in opposite hemisphere, we're done. 
          if lt(max_loc_L, 2)
              loc_ = loc_R;
%               side = log2(max_loc_R/max_loc_L); % channel based
          elseif lt(max_loc_R, 2)
              loc_ = loc_L;
%               side = log2(max_loc_R/max_loc_L);  % channel based
          else  % both sides have at least 2 channels involved.
            
              % compare amplitudes
              % brute force method
              
              % leads with >0 events
%               theseLeads = gt(max(LOC,2),0);
              theseIndices = sum(LOC,1)>=2;
              lpds_ = regionprops(theseIndices, 'Area','PixelIdxList');
              lpds = lpds_(cat(1,lpds_(:).Area)>thresh_samples);
              
              lpds_new = cell(length(lpds),1);
            
              %calculate asymmetry index for each lpd
              for i=1:size(lpds,1)
                  
                  thisLPD = lpds(i);
                  
                  %compute burden as AUC from all channels
                  burden_L_ = sum(trapz(abs(seg_bi(channels_L,thisLPD.PixelIdxList))),2);
                  burden_R_ = sum(trapz(abs(seg_bi(channels_R,thisLPD.PixelIdxList))),2);
                  %               burden = log2(burden_R_/burden_L_);
                  burden_asym = (burden_R_-burden_L_)/(burden_R_+burden_L_);
                  burden_asym2 = log2(burden_R_/burden_L_);
                  
                  thisLPD.burden_asym = burden_asym;
                  thisLPD.burden_asym2 = burden_asym2;
                                    
                  %compute asym index from difference between pairs of R vs L channels,
                  % then take average
                  burden_L_ = abs(seg_bi(channels_L,thisLPD.PixelIdxList));
                  burden_R_ = abs(seg_bi(channels_R,thisLPD.PixelIdxList));
                  burden_asym3 = sum(burden_R_-burden_L_,2)./sum(burden_R_+burden_L_,2);
                  burden_asym4 = max(abs(burden_asym3));
                                    
                  thisLPD.burden_asym3 = burden_asym3;
                  thisLPD.burden_asym4 = burden_asym4;
                  
                  %define LPD as events where asym in at least 1 pair
                  %exceeds 50% (--> abs value exceeding 1/3 wrt asymmetry index) 
%                   thisLPD.isLPD = (sum(burden_asym3>1/3)>0 || sum(burden_asym3<-1/3)>0);
                  thisLPD.isLPD = burden_asym4>1/3;
                                    
                  lpds_new{i} = thisLPD;
              end
              
              % convert to struct array for ease of processing.
              headings = fieldnames(lpds_new{1});
              lpds_new2 = struct();
              for j=1:length(headings)
                  lpds_new2.(headings{j})='';
              end
              for k=1:length(lpds_new)
                  lpds_new2(k)=lpds_new{k};
              end
            
%             isLPD = cat(1,lpds_new2.isLPD);
%             lpd_ids = cat(1,lpds_new2(isLPD).PixelIdxList);
            
            lpd_ids = cat(1,lpds_new2(:).PixelIdxList);
            loc_ = zeros(1,size(seg_bi,2));
            loc_(lpd_ids) = 1;
                             
            cat(2,lpds_new2.burden_asym3)
          end
          
      else
          loc_ = zeros(1,size(LOC,2));   % no events that meet LPD criteria.
          lpds_new2 = NaN;
      end
      
          
    end