function [ thisImage_interp ] = interpolateImage( thisImage_io )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% [xlen,ylen] = size(thisImage_io(:,:,1));
[ylen,xlen] = size(thisImage_io(:,:,1));

% ylen = double(thisWell_bb(4)+3);
% xlen = double(thisWell_bb(3)+3);

preciseFactor = 10;
[X1,Y1] = meshgrid(1:xlen,1:ylen);
[Xq,Yq] = meshgrid(1:1/preciseFactor:xlen,1:1/preciseFactor:ylen);

thisImage_interp = zeros(size(Xq,1),size(Yq,2),size(thisImage_io,3));

for k =1:size(thisImage_io,3)
% thisImage_interp(:,:,k) = interp2(X1, Y1,double(thisImage(:,:,k)), Xq, Yq,'cubic');
thisImage_interp(:,:,k) = interp2(X1, Y1,double(thisImage_io(:,:,k)), Xq, Yq,'cubic');

end

end

