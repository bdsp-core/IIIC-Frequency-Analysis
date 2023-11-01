function [ output_args ] = fct_vizLOC_bipolar( LOC )
%fct_vizLOC_bipolar Visualize LOC from bipolar events
% INPUT
%   LOC
% OUTPUT
%   none. Generates figures. 

channels_L = [1:4,9:12];
channels_R = [5:8, 13:16];

        im_LOC = LOC;
        im_LOC(channels_L,:)=im_LOC(channels_L,:)*3;
        im_LOC([17,18],:)=im_LOC([17,18],:)*2;
        
        figure(26);imagesc(im_LOC([1:4,9:12,5:8,13:18],:));
        title('LOC by channel, left v right v center')
        
        im_LOC2 = [sum(LOC([1:4,9:12],:),1)>2; sum(LOC([5:8,13:16],:),1)>2; ...
            sum(LOC([17,18],:),1)>1];
        figure(28);imagesc(im_LOC2); 
        title('summed LOC >2, L v R v C');
        
        im_LOC2 = [sum(LOC([1:4,9:12],:),1); sum(LOC([5:8,13:16],:),1); ...
            sum(LOC([17,18],:),1)];
        figure(27);imagesc(im_LOC2);
        title('summed LOC, L v R v C');
        
        figure(29);
        subplot(3,1,1); plot(1:size(LOC,2),im_LOC2(1,:)); axis([0 size(LOC,2) 0 8]);
        subplot(3,1,2); plot(1:size(LOC,2), im_LOC2(2,:)); axis([0 size(LOC,2) 0 8]);
        subplot(3,1,3); plot(1:size(LOC,2), im_LOC2(1,:)-im_LOC2(2,:)); axis([0 size(LOC,2) -8 8]);
   
        suptitle('summed LOC, L v R');
                


end

