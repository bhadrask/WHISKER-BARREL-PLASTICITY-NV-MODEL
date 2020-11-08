function [lissom] = rf_m(lissom,inpsel,ii,opts)

% %% Moving Window with given stride
%-----------------------------------
p=1;stride=lissom.layers{ii}.stride;  inpsel1=inpsel;

%all=0; % make all=0 for using receptive field rf.
  
%% all to all
if opts.all(ii-1)==1
for l = 1 :lissom.layers{ii}.dim(1) 
     
    for k = 1 :lissom.layers{ii}.dim(2) 

 lissom.layers{ii}.X1{l,k} =inpsel;
   
    end
end

elseif opts.all(ii-1)==0
%% zero padding 

% temp1=size(inpsel,1)-(lissom.layers{ii}.dim(1)+lissom.layers{ii}.rf(1)-stride);
% temp2=size(inpsel,2)-(lissom.layers{ii}.dim(2)+lissom.layers{ii}.rf(2)-stride);
temp1=size(inpsel,1)-((lissom.layers{ii}.dim(1)-1)*stride+lissom.layers{ii}.rf(1));
temp2=size(inpsel,2)-((lissom.layers{ii}.dim(2)-1)*stride+lissom.layers{ii}.rf(2));
pad1=0;
pad2=0;
if temp1<0
    pad1=ceil(abs(temp1/2));
end
if temp2<0
    pad2=ceil(abs(temp2/2));
end
inpsel1=padarray(inpsel,[pad1 pad2]);

for l = 1 :lissom.layers{ii}.dim(1) 
        q=1;
    for k = 1 :lissom.layers{ii}.dim(2) 
        lissom.layers{ii}.X1{l,k} =inpsel1(p:p+(lissom.layers{ii}.rf(2)-1),q:q+(lissom.layers{ii}.rf(1)-1));
        q=k*stride+1; 
     end
     p=l*stride+1;
end

%% for feedback from layer 2 BSK
if ii>1 && ii<numel(lissom.layers) && opts.control==1
    stride=lissom.layers{ii+1}.stride;
    p=1;
% temp1=lissom.layers{ii+1}.dim(1)-(lissom.layers{ii}.dim(1)+lissom.layers{ii+1}.rf(1)-stride);
% temp2=lissom.layers{ii+1}.dim(2)-(lissom.layers{ii}.dim(2)+lissom.layers{ii+1}.rf(2)-stride);
temp1=lissom.layers{ii+1}.dim(1)-((lissom.layers{ii}.dim(1)-1)*stride+lissom.layers{ii+1}.rf(1));
temp2=lissom.layers{ii+1}.dim(2)-((lissom.layers{ii}.dim(2)-1)*stride+lissom.layers{ii+1}.rf(2));
pad1=0;
pad2=0;
if temp1<0
    pad1=ceil(abs(temp1/2));
end
if temp2<0
    pad2=ceil(abs(temp2/2));
end
inpsel2=padarray(lissom.layers{ii+1}.Zabove,[pad1 pad2]);


for l = 1 :lissom.layers{ii}.dim(1) 
        q=1;
    for k = 1 :lissom.layers{ii}.dim(2) 
        lissom.layers{ii}.X2{l,k} =inpsel2(p:p+(lissom.layers{ii+1}.rf(2)-1),q:q+(lissom.layers{ii+1}.rf(1)-1));
        q=k*stride+1; 
     end
     p=l*stride+1;
end

end
end
%  figure(3);
%  subplot(2,1,1);imagesc(inpsel);subplot(2,1,2); imagesc(cell2mat(lissom.layers{ii}.X1));pause(0.5);title('RF');
