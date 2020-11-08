function [lissom] = rf_m(lissom,inpsel,ii,opts)
p=1;  stride=lissom.layers{ii}.stride; 
if opts.all(ii-1)==1           %% all to all opts.all(2-1) = opts.all(1) = 0
    for l = 1 :lissom.layers{ii}.dim(1) 
        for k = 1 :lissom.layers{ii}.dim(2) 
            lissom.layers{ii}.X1{l,k} =inpsel;
        end
    end
elseif opts.all(ii-1)==0       %% rf to all opts.all(2-1) = opts.all(1) = 0
    %% zero padding  % n2=(n1-rf)/stride + 1  =>(n2-1)*stride = n1-rf   => (n2-1)*stride - n1 + rf = 0
    temp1=size(inpsel,1)-((lissom.layers{ii}.dim(1)-1)*stride+lissom.layers{ii}.rf(1));  %temp1 = -1 
    temp2=size(inpsel,2)-((lissom.layers{ii}.dim(2)-1)*stride+lissom.layers{ii}.rf(2));  %temp2 = -1 rf=n1-[(n2-1)*stride]
    pad1=0;  %32-[(10-1)*3 + 6]  = n1-[(n2-1)*stride + rf] = n1 - (n2-1)*stride - rf
    pad2=0;
    if temp1<0  %yes
        pad1=ceil(abs(temp1/2)); %-1/2 =abs(-.05 )=ceil(0.5)=1
    end                          %pad1=1
    if temp2<0  %yes
        pad2=ceil(abs(temp2/2)); %-1/2 =abs(-.05 )=ceil(0.5)=1
    end                          %pad2=1
    inpsel1=padarray(inpsel,[pad1 pad2]); %[1 1]
    %%rf selection to X1
    for l = 1 :lissom.layers{ii}.dim(1)     %1 to 10 %[4:9] to [1:6], [4:9] to [4:9], [4:9] to [7:12]....
        q=1;
        for k = 1 :lissom.layers{ii}.dim(2) %1 to 10 %[1:6] to [1:6], [4:9] to [1:6], [7:12] to [1:6]....
            lissom.layers{ii}.X1{l,k} =inpsel1(p:p+(lissom.layers{ii}.rf(2)-1),q:q+(lissom.layers{ii}.rf(1)-1));
            q=k*stride+1;   %q=1,4,7,11...
        end
        p=l*stride+1;       %p=1,4,7,11...
    end
end
% figure(7);imagesc(cell2mat(lissom.layers{2}.X1));title(strcat('rf',num2str(2)));pause(0.5); %colormap(gray);
end