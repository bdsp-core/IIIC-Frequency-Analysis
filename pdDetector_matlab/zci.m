function [ indx_zo ] = zci( v )
%zci Indices to zero-crossings 

indx_zo = find(v(:).*circshift(v(:), [-1 0]) <= 0);   

end

