function [ thisImage_interp ] = interpolateImageNan( thisImage_io )
%UNTITLED Like interpolateImage but handles situations where there are NaN
% values at certain 'pixels'

[xlen,ylen] = size(thisImage_io(:,:,1));

% ylen = double(thisWell_bb(4)+3);
% xlen = double(thisWell_bb(3)+3);

preciseFactor = 6;
[X1,Y1] = meshgrid(1:xlen,1:ylen);
[Xq,Yq] = meshgrid(1:1/preciseFactor:xlen,1:1/preciseFactor:ylen);

thisImage_interp = zeros(size(Xq,1),size(Yq,2),size(thisImage_io,3));

%are there nans?

thisImage_io_f = zeros(size(thisImage_io));
% 
% for k=1:size(thisImage_io,3)
%    thisImage_io_f(:,:,k) = inpaint_nans(thisImage_io(:,:,k),1);
% end

for k =1:size(thisImage_io,3)
    thisImage_io_f(:,:,k) = inpaint_nans(thisImage_io(:,:,k),1);
    thisImage_interp(:,:,k) = interp2(X1, Y1,double(thisImage_io_f(:,:,k)), Xq, Yq,'cubic');

end

end

