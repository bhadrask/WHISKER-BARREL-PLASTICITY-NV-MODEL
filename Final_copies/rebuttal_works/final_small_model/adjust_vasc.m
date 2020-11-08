function [lissom] = adjust_vasc(lissom,opts)
str_fac=floor(lissom.layers{2}.dim(1)/lissom.layers{3}.dim(1)); %made 1 change here
excess=ceil((lissom.layers{2}.dim(1)-str_fac*lissom.layers{3}.dim(1))/str_fac);
if excess>0
first_indx=str_fac*ceil(excess/2);
last_indx=first_indx+str_fac*(lissom.layers{3}.dim(1)-1);
else
    first_indx=1;
    last_indx=lissom.layers{2}.dim(1);
end
[cx,cy]=meshgrid(first_indx:str_fac:last_indx,first_indx:str_fac:last_indx);
lissom.layers{2}.Adjusted_vasc=zeros(lissom.layers{2}.dim(1));
lissom.layers{2}.Adjusted_vasc(cx(:),cy(:))=1;
 lissom.layers{2}.vessel_adjust_indx.cx=cx;
  lissom.layers{2}.vessel_adjust_indx.cy=cy;
% [a,b]=find(Adjusted_vasc>0);
% lissom.layers{2}.vessel_adjust_indx=sub2ind(size(Adjusted_vasc),a,b);


